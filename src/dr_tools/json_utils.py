# This file contains utility functions for working with JSON files
import json
import warnings
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Literal, overload

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from ddbj_record.converter.v1_to_v2 import v1_to_v2
from ddbj_record.converter.v2_to_v1 import v2_to_v1
from ddbj_record.converter.v2_to_v3 import v2_to_v3
from ddbj_record.schema import normalize_schema_version
from ddbj_record.schema.v1 import DdbjRecord as DdbjRecordV1
from ddbj_record.schema.v2 import DdbjRecord as DdbjRecordV2
from ddbj_record.schema.v3 import DdbjRecord as DdbjRecordV3
from ddbj_record.schema.v3 import Entry, Feature, Qualifier, Source

# MSS の qualifier の形。値の無い qualifier (/pseudo など) は True で表す
QualifierValues = dict[str, list[str | bool]]

DEFAULT_LOCUS_TAG_PREFIX = "LOCUS"

RecordVersion = Literal["v1", "v2", "v3"]


def _major_version(schema_version: Any) -> str | None:
    """
    schema_version の値から major ("v1" / "v2" / "v3") を決める。分からなければ None
    v3 は minor を持たず "v3" だけを受け付ける。v1 / v2 は ddbj_record の normalize_schema_version で揃える
    """
    if not isinstance(schema_version, str):
        return None
    if schema_version == "v3":
        return "v3"
    normalized = normalize_schema_version(schema_version)
    return normalized.split(".")[0] if normalized else None


@overload
def load_json_to_ddbj_record_instance(json_file: Path, to_record_version: Literal["v1"]) -> DdbjRecordV1: ...


@overload
def load_json_to_ddbj_record_instance(json_file: Path, to_record_version: Literal["v2"]) -> DdbjRecordV2: ...


@overload
def load_json_to_ddbj_record_instance(json_file: Path, to_record_version: Literal["v3"] = "v3") -> DdbjRecordV3: ...


def load_json_to_ddbj_record_instance(
    json_file: Path,
    to_record_version: RecordVersion = "v3",
) -> DdbjRecordV1 | DdbjRecordV2 | DdbjRecordV3:
    """
    JSONファイルを読み込み、to_record_version の DdbjRecord インスタンスを生成して返す

    v1 の JSON は v1 → v2 → v3、v2 の JSON は v2 → v3 と ddbj_record の converter で進める。
    v3 の JSON を v1 / v2 に戻す converter は無いので、その要求は ValueError にする。

    schema_version は ddbj_record の normalize_schema_version でメジャー版に畳む。
    ddbj_record は書き出しの際に schema_version を最新のマイナー版 ("v1.0", "v2.2" 等)
    へ正規化するので、旧来の表記 ("0.1", "v1") だけを受け付けると
    dr_tools が書いた JSON を dr_tools 自身が読めなくなる。
    """
    with open(json_file, encoding="utf-8") as f:
        raw_data = json.load(f)

    schema_version = raw_data.get("schema_version", "")
    major_version = _major_version(schema_version)

    if major_version == "v3":
        if to_record_version != "v3":
            raise ValueError(f"A v3 record cannot be converted to {to_record_version}")
        return DdbjRecordV3.model_validate(raw_data)

    record_v2_instance: DdbjRecordV2
    if major_version == "v1":
        record_v1_instance = DdbjRecordV1.model_validate(raw_data)
        if to_record_version == "v1":
            return record_v1_instance
        record_v2_instance = v1_to_v2(record_v1_instance)
    elif major_version == "v2":
        record_v2_instance = DdbjRecordV2.model_validate(raw_data)
        if to_record_version == "v1":
            return v2_to_v1(record_v2_instance)
    else:
        raise ValueError(f"Unsupported schema_version: {schema_version}")

    if to_record_version == "v2":
        return record_v2_instance
    return v2_to_v3(record_v2_instance)


# === v3 の record を読むための関数 ===


def qualifiers_to_values(qualifiers: Mapping[str, list[Qualifier]] | None) -> QualifierValues:
    """
    v3 の qualifier を MSS の形 ({名前: [値, ...]}) にする
    value を持たない qualifier (/pseudo など) は True にする
    """
    return {
        name: [True if qualifier.value is None else qualifier.value for qualifier in values]
        for name, values in (qualifiers or {}).items()
    }


def values_to_qualifiers(values: Mapping[str, list[str | bool]]) -> dict[str, list[Qualifier]] | None:
    """
    MSS の形の qualifier を v3 の qualifier にする
    True (値の無い qualifier) は value を省く。qualifier が 1 つも無ければ None を返す
    """
    if not values:
        return None
    return {
        name: [Qualifier() if isinstance(value, bool) else Qualifier(value=value) for value in value_list]
        for name, value_list in values.items()
    }


def merge_source_qualifiers(common_source: Source | None, source: Source | None) -> QualifierValues:
    """
    entry の source で common_source を上書きし、MSS の source feature の qualifier にする
    上書きの単位は organism、mol_type、qualifier の名前ごとで、entry の source に無いものは common_source の値を使う
    organism と mol_type は qualifier の先頭に置く
    """
    merged: QualifierValues = {}

    organism = None
    mol_type = None
    for candidate in (common_source, source):
        if candidate is None:
            continue
        if candidate.organism is not None and candidate.organism.name is not None:
            organism = candidate.organism.name
        if candidate.mol_type is not None:
            mol_type = candidate.mol_type
    if organism is not None:
        merged["organism"] = [organism]
    if mol_type is not None:
        merged["mol_type"] = [mol_type]

    for candidate in (common_source, source):
        if candidate is None:
            continue
        for name, value_list in qualifiers_to_values(candidate.qualifiers).items():
            merged[name] = value_list

    return merged


def features_by_entry(record: DdbjRecordV3) -> dict[str, list[Feature]]:
    """
    record の features を entry の alias (sequence_id) ごとにまとめる。並びは record の順を保つ
    sequence_id を持たない feature はどの entry にも置けないので、警告を出して除く
    """
    grouped: dict[str, list[Feature]] = {}
    for feature in record.features or []:
        if feature.sequence_id is None:
            warnings.warn(f"feature {feature.alias!r} has no sequence_id; skipped", UserWarning, stacklevel=2)
            continue
        grouped.setdefault(feature.sequence_id, []).append(feature)
    return grouped


def get_entry_id(entry: Entry) -> str:
    """
    entry の識別子 (MSS の 1 列目と FASTA のヘッダになる名前) を返す
    alias を使い、無ければ name を使う。どちらも無ければ MSS に書けないので ValueError にする
    """
    entry_id = entry.alias or entry.name
    if not entry_id:
        raise ValueError("entry has neither alias nor name")
    return entry_id


def get_locus_tag_prefix(record: DdbjRecordV3) -> str:
    """
    record から locus_tag の prefix を取得する
    projects[].locus_tag_prefix[].prefix の最初のものを使い、無ければ "LOCUS" を返す
    """
    for project in record.projects or []:
        for locus_tag_prefix in project.locus_tag_prefix or []:
            if locus_tag_prefix.prefix:
                return locus_tag_prefix.prefix
    return DEFAULT_LOCUS_TAG_PREFIX


def get_locus_tag(feature: Feature, locus_tag_prefix: str) -> str | None:
    """
    locus_tag_prefix と feature の locus_tag_id を結合した locus_tag を返す
    locus_tag_id を持たない feature では None を返す
    """
    if feature.locus_tag_id is None:
        return None
    return f"{locus_tag_prefix}_{feature.locus_tag_id}"


def feature_qualifier_values(feature: Feature, locus_tag_prefix: str) -> QualifierValues:
    """
    feature の qualifier を MSS の形にし、locus_tag_id があれば locus_tag の qualifier を足す
    """
    qualifiers = qualifiers_to_values(feature.qualifiers)
    locus_tag = get_locus_tag(feature, locus_tag_prefix)
    if locus_tag is not None:
        qualifiers["locus_tag"] = [locus_tag]
    return qualifiers


def get_feature_and_entry(record: DdbjRecordV3, feature_id: str) -> tuple[Feature, Entry]:
    """
    alias が feature_id の feature と、その feature が載る entry を返す
    entry は features_by_entry と同じく get_entry_id (alias、無ければ name) が feature の sequence_id と一致するものを探す
    """
    features = [feature for feature in record.features or [] if feature.alias == feature_id]
    if len(features) == 0:
        raise ValueError(f"Feature with id {feature_id} not found in JSON data")
    if len(features) > 1:
        raise ValueError(f"Multiple features with id {feature_id} found in JSON data")
    feature = features[0]
    if feature.sequence_id is None:
        raise ValueError(f"Feature {feature_id} has no sequence_id")
    entries = [
        entry
        for entry in (record.sequences.entries if record.sequences else None) or []
        if (entry.alias or entry.name) == feature.sequence_id
    ]
    if len(entries) == 0:
        raise ValueError(f"Entry {feature.sequence_id} of feature {feature_id} not found in JSON data")
    return feature, entries[0]


def get_feature_json(json_file: Path, feature_id: str) -> dict[str, Any]:
    """
    JSONデータからfeature_idに対応するfeatureの情報を辞書として取得する
    DFAST web serviceの feature 詳細で表示される情報を取得する

    v3 の feature の項目に加えて、locus_tag の qualifier、nucleotide (feature の塩基配列)、
    translation (CDS のアミノ酸配列) を持つ
    """
    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v3")
    feature, entry = get_feature_and_entry(record, feature_id)
    locus_tag_prefix = get_locus_tag_prefix(record)

    feature_json: dict[str, Any] = feature.model_dump(exclude_none=True)
    locus_tag = get_locus_tag(feature, locus_tag_prefix)
    if locus_tag is not None:
        feature_json.setdefault("qualifiers", {})["locus_tag"] = [{"value": locus_tag}]

    # SeqFeatureオブジェクトを作成 (CDSの場合は、qualifiersにtranslationを追加)
    seq_feature, nucleotide = json_to_seqfeature(feature, entry, locus_tag_prefix)
    feature_json["nucleotide"] = str(nucleotide)
    translate = seq_feature.qualifiers.get("translation", [""])[0]
    feature_json["translation"] = translate
    return feature_json


def json_to_seqfeature(feature: Feature, entry: Entry, locus_tag_prefix: str) -> tuple[SeqFeature, Seq]:
    """
    v3 の feature を BioPython SeqFeature object に変換し、
    そのオブジェクトと、その feature location が指す配列を Seq objectとして返す
    """
    # json2biopython が json_utils を import するので、ここで import して循環を避ける
    from dr_tools.json2biopython import add_translate_qualifier, create_seqfeature  # noqa: PLC0415

    seq = Seq(entry.sequence or "")
    topology = entry.topology or "linear"
    qualifiers = feature_qualifier_values(feature, locus_tag_prefix)

    seq_feature = create_seqfeature(feature.type or "", feature.location or "", qualifiers, len(seq), topology)
    if feature.type == "CDS":
        add_translate_qualifier(seq_feature, seq)
    nucleotide: Seq = seq_feature.extract(seq)  # type: ignore[no-untyped-call]
    return seq_feature, nucleotide
