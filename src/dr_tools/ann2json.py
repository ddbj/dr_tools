import argparse
import logging
import re
import sys
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal, TypeVar

from ddbj_record.converter.v1_to_v2 import v1_to_v2
from ddbj_record.schema.v1 import DdbjRecord as DdBjRecordV1
from ddbj_record.schema.v2 import DdbjRecord as DdBjRecordV2
from ddbj_record.schema.v3 import (
    Address,
    LocusTagPrefix,
    Organism,
    Organization,
    Person,
    Project,
    Provenance,
    Publication,
    Relation,
    RelationTarget,
    Sequences,
    Source,
    SourceFeature,
    StructuredComment,
    Submission,
)
from ddbj_record.schema.v3 import DdbjRecord as DdbjRecordV3
from ddbj_record.schema.v3 import Entry as EntryV3
from ddbj_record.schema.v3 import Feature as FeatureV3
from pydantic import BaseModel

from dr_tools.json_utils import QualifierValues, values_to_qualifiers
from dr_tools.MSS import (
    DBLINK_KEY_TO_DB,
    MSS,
    MSS_DIVISIONS,
    Entry,
    Feature,
    hold_date_to_record,
    reference_status_to_record,
)

DATA_MODEL_VERSION = "0.1"

RecordVersion = Literal["v1", "v2", "v3"]

# record v3 の provenance.source_format に書く値
SOURCE_FORMAT_MSS = "MSS"

# entry の source の qualifier のうち、全ての entry で同じ値でも common_source に移さないもの
ENTRY_SPECIFIC_SOURCE_QUALIFIERS = ("submitter_seqid", "ff_definition")

# Entry.type を決める source の qualifier。先に書いたものを優先する
ENTRY_TYPE_QUALIFIERS = ("plasmid", "chromosome", "organelle", "segment")

_M = TypeVar("_M", bound=BaseModel)

logger = logging.getLogger(__name__)


# === ann の COMMON の読み取り (v1 / v2 用。v3 でも流用する) ===


def get_common_entry(mss: MSS) -> Entry:
    """COMMON entry を返す。無い、または 2 つ以上あれば ValueError"""
    common_entries = [entry for entry in mss.entries if entry.id == "COMMON"]
    if len(common_entries) == 0:
        raise ValueError("COMMON entry not found.")
    if len(common_entries) > 1:
        raise ValueError("Multiple COMMON entries found.")
    return common_entries[0]


def get_dfast_version(entries: list[Entry]) -> str | None:
    """
    DFASTのバージョンを取得する
    COMMENT		line	Annotated by DFAST v.1.3.4 https://dfast.ddbj.nig.ac.jp/
    の1.3.4の部分を取得する
    """
    dfast_version_pat = re.compile(r"DFAST (ver\.|ver|version|v|v\.)?\s?(\d+\.\d+\.\d+)")
    for entry in entries:
        if entry.id == "COMMON":
            for feature in entry.features:
                if feature.type == "COMMENT":
                    for qualifier_value in feature.qualifiers.get("line", []):
                        if isinstance(qualifier_value, bool):
                            continue
                        m = dfast_version_pat.search(qualifier_value)
                        if m:
                            return m.group(2)
    return None


def get_locus_tag_prefix_from_ann(entries: list[Entry]) -> str | None:
    """
    locus_tag_prefix の決定
    entry, featureのqualifierにlocus_tagが含まれる場合、その値の _ より前の部分をlocus_tag_prefixとする
    """
    for entry in entries:
        for feature in entry.features:
            if "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                if locus_tag and isinstance(locus_tag, str) and "_" in locus_tag:
                    return locus_tag.split("_")[0]
    return None


def get_submission_category(mss: MSS) -> str:
    """
    登録区分 (WGS / GNM) の判定
    COMMON entry の DATATYPE feature の type が WGS であれば WGS、それ以外は GNM とする
    """
    for entry in mss.entries:
        if entry.id == "COMMON":
            for feature in entry.features:
                if feature.type == "DATATYPE" and "WGS" in feature.qualifiers.get("type", []):
                    return "WGS"
            break
    return "GNM"


def get_division(common_entry: Entry) -> str | None:
    """
    COMMON entry の DIVISION feature の division。MSS_DIVISIONS に無い値 (TSA、CON など) と、DIVISION が無いときは None
    MSS.from_record が書くのは MSS_DIVISIONS の値だけなので、読むのもその値だけにする。他の値を読むと、
    引数の division より優先され、ann に戻すと消える
    """
    for feature in common_entry.features:
        if feature.type == "DIVISION":
            values = _text_values(feature, "division")
            if values and values[0] in MSS_DIVISIONS:
                return values[0]
    return None


def get_seq_prefix(mss: MSS, submission_category: str) -> str | None:
    """
    seq_prefix の決定 (WGS の場合のみ)
    最初の entry の名前 (sequence01, contig01 等) から末尾の数字を取り除いたもの
    """
    if submission_category != "WGS":
        return None
    for entry in mss.entries:
        if entry.id == "COMMON":
            continue
        return entry.name.rstrip("0123456789") or None
    return None


def common_to_dict(mss: MSS) -> dict[str, dict[str, Any]]:  # Dict[str, V1Common]
    # COMMON entryの情報を取得して辞書で返す
    qualifiers_with_multi_values = ["line", "ab_name", "keyword", "sequence read archive"]
    features_with_multi_values = ["COMMENT", "REFERENCE"]
    common_entry = get_common_entry(mss)
    common: dict[str, dict[str, Any]] = {}
    for feature in common_entry.features:
        inner_dict: dict[str, Any] = {}
        for key, value in feature.qualifiers.items():
            if key in qualifiers_with_multi_values:
                inner_dict[key] = value
            else:
                assert len(value) == 1
                inner_dict[key] = value[0]
        if feature.type in features_with_multi_values:
            common.setdefault(feature.type, []).append(inner_dict)  # type: ignore
        else:
            common[feature.type] = inner_dict
    return {"COMMON": common}


def source_to_dict(mss: MSS) -> dict[str, dict[str, str | list[str]]]:  # Dict[str, V1CommonSource]
    """
    source featureの情報を取得して、共通する項目をCOMMON_SOURCEとして辞書で返す
    COMMON_SOURCEに含まれるqualifierは各entryのsource featureからは削除する
    """
    qualifiers_with_multi_values = ["note"]
    ignore_qualifiers = ["submitter_seqid", "ff_definition"]
    source_features = [feature for entry in mss.entries for feature in entry.features if feature.type == "source"]

    # すべてのentryについて等しい値を持つqualifierを抽出する
    common_qualifiers: dict[str, list[str]] = {}
    for feature in source_features:
        if feature.type != "source":
            continue
        for key, value in feature.qualifiers.items():
            if key not in ignore_qualifiers:
                str_value = "//".join([str(v) for v in value])  # 一時的にリストを文字列に変換
                common_qualifiers.setdefault(key, []).append(str_value)
    common_qualifiers = {
        key: value
        for key, value in common_qualifiers.items()
        if len(value) == len(source_features) and len(set(value)) == 1
    }

    # 各entryのsource featureからCOMMON_SOURCEに含まれるqualifierを削除
    common_qualifiers_keys = set(common_qualifiers.keys())
    for source_feature in source_features:
        for key in common_qualifiers_keys:
            if key in source_feature.qualifiers:
                del source_feature.qualifiers[key]

    # COMMON_SOURCEとして返す辞書を作成
    common_source: dict[str, str | list[str]] = {}
    for c_key, c_value in common_qualifiers.items():
        first_c_value = c_value[0]
        if c_key in qualifiers_with_multi_values:
            common_source[c_key] = first_c_value.split("//")
        else:
            common_source[c_key] = first_c_value
    return {"COMMON_SOURCE": common_source}


SeqInfo = dict[str | int, dict[str, str]]  # {entry_id: {"seq_type": str, "seq_topology": str}}


def infer_meta_from_ann(mss: MSS) -> tuple[dict[str, dict[str, str]], SeqInfo]:
    """
    annファイルに記載された情報から、
    trad_submission_category: GNM; WGS
    seq_type: chromosome; plasmid; unplaced; other
    seq_topology: linear; circular
    を決定する
    common_meta と seq_info の2つの辞書を返す
    注意: DFASTが出力するannファイルにのみ対応
    """
    dfast_version = get_dfast_version(mss.entries)
    # For DFAST, the division is always BCT
    common_meta = {"dfast_version": dfast_version, "division": "BCT"} if dfast_version else {"division": "UNK"}

    trad_submission_category = get_submission_category(mss)
    common_meta["trad_submission_category"] = trad_submission_category

    seq_prefix = get_seq_prefix(mss, trad_submission_category)
    if seq_prefix:
        common_meta["seq_prefix"] = seq_prefix

    locus_tag_prefix = get_locus_tag_prefix_from_ann(mss.entries)
    if locus_tag_prefix:
        common_meta["locus_tag_prefix"] = locus_tag_prefix

    # seq_infoの決定
    # seq_type は ff_definition の文言から読むので、書式が決まっている gnm (complete) でのみ判定する
    # seq_topology は TOPOLOGY 行の有無だけで決まる。wgs (draft) の配列も circular になりうるので、
    # 登録区分によらず読む (ここで落とすと ann -> json -> ann で TOPOLOGY 行が消える)
    seq_info: SeqInfo = {}
    for entry in mss.entries:
        source_or_topology = [feature for feature in entry.features if feature.type in ["source", "TOPOLOGY"]]
        seq_topology, seq_type = "linear", "other"
        for feature in source_or_topology:
            if feature.type == "source" and trad_submission_category == "GNM":
                ff_definition: str = feature.qualifiers.get("ff_definition", [""])[0]  # type: ignore
                if "plasmid" in ff_definition:
                    seq_type = "plasmid"
                elif "unplaced sequence" in ff_definition:
                    seq_type = "unplaced"
                elif "complete genome" in ff_definition:
                    seq_type = "chromosome"

            # seq_topologyの決定 TOPOLOGY featureが記載されていればcircularとする
            if feature.type == "TOPOLOGY" and "circular" in feature.qualifiers:
                seq_topology = "circular"
        seq_info[entry.id] = {"seq_type": seq_type, "seq_topology": seq_topology}

    return {"COMMON_META": common_meta}, seq_info


def make_entry_dict(mss: MSS, seq_info: SeqInfo) -> dict[str, list[dict[str, Any]]]:
    """
    mssオブジェクトを受け取り、各entryの情報を辞書で返す
    """
    entry_dicts = []
    sequence_dict = {sequence.id: sequence.seq for sequence in mss.sequences}
    for entry in mss.entries:
        if entry.id == "COMMON":
            continue
        # note: sequence.id and entry.id is always identical, but entry.name may be different.
        sequence = sequence_dict[entry.id]
        seq_info_dict = seq_info.get(entry.id, {})
        seq_type, seq_topology = seq_info_dict.get("seq_type", "other"), seq_info_dict.get("seq_topology", "linear")
        features = [feature.to_dict() for feature in entry.features if feature.type != "TOPOLOGY"]
        entry_dicts.append(
            {
                "id": entry.id,
                "name": entry.name,
                "type": seq_type,
                "topology": seq_topology,
                "sequence": sequence,
                "features": features,
            }
        )
    return {"ENTRIES": entry_dicts}


# === MSS から record v3 を組み立てる ===


def _or_none(model: _M) -> _M | None:
    """値を 1 つも持たないモデルは省く"""
    return model if model.model_dump(exclude_none=True) else None


def _text_values(feature: Feature, key: str) -> list[str]:
    """qualifier の値のうち、空でない文字列だけを返す"""
    return [value for value in feature.qualifiers.get(key, []) if isinstance(value, str) and value != ""]


def _text_value(feature: Feature, key: str) -> str | None:
    """qualifier の最初の空でない文字列の値を返す。無ければ None"""
    values = _text_values(feature, key)
    return values[0] if values else None


def _normalize_abbr(abbr: str) -> str:
    return re.sub(r"[.\-\s]", "", abbr).lower()


def _abbr_candidates_from_fullname(fullname: str) -> set[str]:
    """
    contact のフルネームから、ab_name として書かれうる略記の候補を作る
    "John A. Doe" なら "Doe,J.A."、"Doe,J."、"Doe," など。姓が先に来る名前 (Yamada Taro) も試す
    """
    name = fullname.strip()
    candidates: set[str] = set()

    def _add_candidate(last: str, initials: str) -> None:
        if initials:
            candidates.add(f"{last},{initials}.")
            candidates.add(f"{last},{initials}")
            candidates.add(f"{last},{initials[0]}.")
            candidates.add(f"{last},{initials[0]}")
            candidates.add(f"{last},{'.'.join(list(initials))}.")
            candidates.add(f"{last},{'.'.join(list(initials))}")
        candidates.add(f"{last},")
        candidates.add(f"{last}")

    if "," in name:  # "Last, First Middle"
        last, rest = [x.strip() for x in name.split(",", 1)]
        rest_parts = [p for p in re.split(r"\s+", rest) if p]
        initials = "".join([p[0] for p in rest_parts if p])
        _add_candidate(last, initials)
    else:  # "First Middle Last"
        parts = [p for p in re.split(r"\s+", name) if p]
        if len(parts) == 1:
            candidates.add(parts[0])
        else:
            first, last = parts[0], parts[-1]
            mids = parts[1:-1]
            initials = first[0] + "".join([m[0] for m in mids if m])
            _add_candidate(last, initials)

            # 姓が先に来る名前 (Yamada Taro -> Last=Yamada, First=Taro)
            last, first = parts[0], parts[-1]
            initials = first[0] + "".join([m[0] for m in mids if m])
            _add_candidate(last, initials)

    return {_normalize_abbr(c) for c in candidates if c}


def find_contact_index(contact: str | None, ab_names: Sequence[str]) -> int | None:
    """contact のフルネームに当たる ab_name の位置を返す。当たるものが無ければ None"""
    if not contact:
        return None
    candidates = _abbr_candidates_from_fullname(contact)
    for i, ab_name in enumerate(ab_names):
        if _normalize_abbr(ab_name) in candidates:
            return i
    return None


def build_submitters(common_entry: Entry) -> list[Person] | None:
    """
    SUBMITTER feature から submitters を作る
    contact に当たる ab_name の人を先頭に置き、name / email / 所属をその人に付ける。先頭が role "contact"、
    それ以外が "submitter"。contact がどの ab_name にも当たらなければ、abbreviation の無い人を先頭に足す
    """
    submitter_features = [feature for feature in common_entry.features if feature.type == "SUBMITTER"]
    if not submitter_features:
        return None
    feature = submitter_features[0]

    ab_names = _text_values(feature, "ab_name")
    contact = _text_value(feature, "contact")
    email = _text_value(feature, "email")

    organizations: list[Organization] = []
    institution = Organization(
        name=_text_value(feature, "institute"),
        url=_text_value(feature, "url"),
        department=_text_value(feature, "department"),
        address=_or_none(
            Address(
                country=_text_value(feature, "country"),
                state=_text_value(feature, "state"),
                city=_text_value(feature, "city"),
                street=_text_value(feature, "street"),
                postal_code=_text_value(feature, "zip"),
            )
        ),
    )
    if _or_none(institution) is not None:
        institution.type = "institution"
        organizations.append(institution)
    consrtm = _text_value(feature, "consrtm")
    if consrtm:
        organizations.append(Organization(name=consrtm, type="consortium"))

    persons = [Person(abbreviation=ab_name) for ab_name in ab_names]
    contact_index = find_contact_index(contact, ab_names)
    if contact_index is not None:
        persons.insert(0, persons.pop(contact_index))
    elif contact:
        warnings.warn(
            f"contact '{contact}' did not match any ab_name; adding as separate Person without abbreviation",
            UserWarning,
            stacklevel=2,
        )
        persons.insert(0, Person())
    elif not persons and (email or organizations):
        persons.append(Person())

    if not persons:
        return None
    persons[0].name = contact
    persons[0].email = email
    persons[0].organizations = organizations or None
    for i, person in enumerate(persons):
        person.role = "contact" if i == 0 else "submitter"
    return persons


def build_publications(common_entry: Entry) -> list[Publication]:
    """REFERENCE feature から publications を作る。year はそのまま date に置く"""
    publications: list[Publication] = []
    for feature in common_entry.features:
        if feature.type != "REFERENCE":
            continue
        status = _text_value(feature, "status")
        publication = _or_none(
            Publication(
                title=_text_value(feature, "title"),
                status=reference_status_to_record(status) if status else None,
                date=_text_value(feature, "year"),
                authors=[Person(abbreviation=ab_name) for ab_name in _text_values(feature, "ab_name")] or None,
                consortiums=_text_values(feature, "consrtm") or None,
            )
        )
        if publication is not None:
            publications.append(publication)
    return publications


def build_relations(common_entry: Entry) -> list[Relation]:
    """DBLINK feature から relations (xref) を作る。DBLINK の並びを保つ"""
    relations: list[Relation] = []
    for feature in common_entry.features:
        if feature.type != "DBLINK":
            continue
        for key, values in feature.qualifiers.items():
            db = DBLINK_KEY_TO_DB.get(key)
            if db is None:
                warnings.warn(f"DBLINK qualifier {key!r} has no place in record; dropped", UserWarning, stacklevel=2)
                continue
            relations.extend(
                Relation(type="xref", target=RelationTarget(db=db, id=value))
                for value in values
                if isinstance(value, str) and value
            )
    return relations


def build_comments(features: list[Feature]) -> list[str]:
    """COMMENT feature を、行を改行でつないだ 1 つの文字列にする。行の無い COMMENT は省く"""
    comments: list[str] = []
    for feature in features:
        if feature.type != "COMMENT":
            continue
        lines = _text_values(feature, "line")
        if lines:
            comments.append("\n".join(lines))
    return comments


def build_structured_comments(common_entry: Entry) -> list[StructuredComment]:
    """ST_COMMENT feature から structured_comments を作る。項目名は ann の書き方のまま fields のキーにする"""
    structured_comments: list[StructuredComment] = []
    for feature in common_entry.features:
        if feature.type != "ST_COMMENT":
            continue
        fields: dict[str, str] = {}
        for key in feature.qualifiers:
            if key == "tagset_id":
                continue
            values = _text_values(feature, key)
            if not values:
                continue
            if len(values) > 1:
                warnings.warn(
                    f"ST_COMMENT field {key!r} has {len(values)} values; only the first is kept",
                    UserWarning,
                    stacklevel=2,
                )
            fields[key] = values[0]
        structured_comment = _or_none(
            StructuredComment(tagset_id=_text_value(feature, "tagset_id"), fields=fields or None)
        )
        if structured_comment is not None:
            structured_comments.append(structured_comment)
    return structured_comments


def build_keywords(common_entry: Entry) -> list[str]:
    """KEYWORD feature の keyword を集める"""
    return [
        keyword
        for feature in common_entry.features
        if feature.type == "KEYWORD"
        for keyword in _text_values(feature, "keyword")
    ]


def get_hold_date(common_entry: Entry) -> str | None:
    """DATE feature の hold_date を record の書き方 (YYYY-MM-DD) で返す"""
    for feature in common_entry.features:
        if feature.type == "DATE":
            hold_date = _text_value(feature, "hold_date")
            if hold_date:
                return hold_date_to_record(hold_date)
    return None


def infer_entry_type(
    source_qualifiers: Mapping[str, Sequence[str | bool]],
    submission_category: str,
    division: str | None = None,
) -> str:
    """
    source の qualifier から Entry.type を決める
    /plasmid、/chromosome、/organelle、/segment の順に見て、あればその名前。どれも無ければ WGS は unplaced、
    division が VRL (ウイルスの配列 1 本のゲノム) は segment、それ以外 (GNM) は chromosome
    """
    for key in ENTRY_TYPE_QUALIFIERS:
        if key in source_qualifiers:
            return key
    if submission_category == "WGS":
        return "unplaced"
    if division == "VRL":
        return "segment"
    return "chromosome"


def infer_completeness(ff_definition: str | None) -> str | None:
    """
    ff_definition の文から Entry.completeness を決める
    "nearly complete" -> nearly complete、"complete genome" / "complete sequence" -> complete、
    "partial genome" / "partial sequence" / "draft genome" / "draft sequence" -> partial。読めなければ None
    """
    if not ff_definition:
        return None
    text = ff_definition.lower()
    if "nearly complete" in text:
        return "nearly complete"
    if "complete genome" in text or "complete sequence" in text:
        return "complete"
    if any(phrase in text for phrase in ("partial genome", "partial sequence", "draft genome", "draft sequence")):
        return "partial"
    return None


def common_source_qualifiers(mss: MSS) -> QualifierValues:
    """
    全ての entry の source feature で同じ値を持つ qualifier を返す (organism / mol_type を含む)
    submitter_seqid と ff_definition は entry 固有なので含めない。MSS は変えない
    """
    source_features = [
        feature
        for entry in mss.entries
        if entry.id != "COMMON"
        for feature in entry.features
        if feature.type == "source"
    ]
    if not source_features:
        return {}
    common: QualifierValues = {}
    first, others = source_features[0], source_features[1:]
    for key, values in first.qualifiers.items():
        if key in ENTRY_SPECIFIC_SOURCE_QUALIFIERS:
            continue
        if all(feature.qualifiers.get(key) == values for feature in others):
            common[key] = list(values)
    return common


def _source_from_qualifiers(qualifiers: Mapping[str, Sequence[str | bool]]) -> Source | None:
    """MSS の source の qualifier から Source を作る。organism と mol_type は qualifiers から分けて持つ"""
    rest: QualifierValues = {}
    organism: str | None = None
    mol_type: str | None = None
    for key, values in qualifiers.items():
        text = next((value for value in values if isinstance(value, str) and value), None)
        if key == "organism":
            organism = text
        elif key == "mol_type":
            mol_type = text
        else:
            rest[key] = list(values)
    return _or_none(
        Source(
            organism=Organism(name=organism) if organism else None,
            mol_type=mol_type,
            qualifiers=values_to_qualifiers(rest),
        )
    )


def _build_feature(feature: Feature, entry_id: str) -> FeatureV3:
    """entry の feature (source / TOPOLOGY / COMMENT 以外) を record の feature にする。locus_tag は locus_tag_id に分ける"""
    qualifiers: QualifierValues = {key: list(values) for key, values in feature.qualifiers.items()}
    locus_tag_id: str | None = None
    locus_tags = qualifiers.pop("locus_tag", [])
    locus_tag = locus_tags[0] if locus_tags else None
    if isinstance(locus_tag, str) and locus_tag:
        if "_" not in locus_tag:
            raise ValueError(f"Invalid locus_tag: {locus_tag}")
        locus_tag_id = locus_tag.split("_", 1)[1]
    return FeatureV3(
        alias=str(feature.id),
        type=feature.type,
        location=feature.location or None,
        sequence_id=entry_id,
        qualifiers=values_to_qualifiers(qualifiers),
        locus_tag_id=locus_tag_id,
    )


def mss_to_record_v3(mss: MSS, division: str | None = None) -> DdbjRecordV3:
    """
    MSS (ann + 配列) から record v3 を組み立てる

    - provenance.source_format は "MSS"、submission_category は DATATYPE から決める (WGS / GNM)
    - DFAST の version が COMMENT にあれば provenance.dfast_version に書く (型に無いキー)
    - division は COMMON の DIVISION が ENV なら ENV を使う。それ以外は引数の division を使い、
      それも省くと、DFAST の出力なら BCT、そうでなければ UNK
    - 全ての entry で同じ値の source の qualifier は sequences.common_source に置き、entry の source には
      entry 固有の qualifier だけを置く
    """
    common_entry = get_common_entry(mss)
    dfast_version = get_dfast_version(mss.entries)
    submission_category = get_submission_category(mss)
    division = get_division(common_entry) or division
    if division is None:
        division = "BCT" if dfast_version else "UNK"

    provenance_data: dict[str, Any] = {
        "source_format": SOURCE_FORMAT_MSS,
        "submission_category": submission_category,
    }
    if dfast_version:
        provenance_data["dfast_version"] = dfast_version

    submission = _or_none(
        Submission(
            submitters=build_submitters(common_entry),
            hold_date=get_hold_date(common_entry),
            comments=build_comments(common_entry.features) or None,
        )
    )

    locus_tag_prefix = get_locus_tag_prefix_from_ann(mss.entries)
    project = _or_none(
        Project(
            publications=build_publications(common_entry) or None,
            locus_tag_prefix=[LocusTagPrefix(prefix=locus_tag_prefix)] if locus_tag_prefix else None,
        )
    )

    common_qualifiers = common_source_qualifiers(mss)
    sequence_dict = {sequence.id: sequence.seq for sequence in mss.sequences}
    entries: list[EntryV3] = []
    features: list[FeatureV3] = []
    for entry in mss.entries:
        if entry.id == "COMMON":
            continue
        if entry.id not in sequence_dict:
            raise ValueError(f"Sequence for entry {entry.id} not found in the sequence file")

        source_features: list[SourceFeature] = []
        entry_type: str | None = None
        completeness: str | None = None
        topology = "linear"
        for feature in entry.features:
            if feature.type == "source":
                entry_specific = {
                    key: values
                    for key, values in feature.qualifiers.items()
                    if key != "ff_definition" and key not in common_qualifiers
                }
                source_features.append(
                    SourceFeature(
                        alias=str(feature.id),
                        location=feature.location or None,
                        source=_source_from_qualifiers(entry_specific),
                        definition=_text_values(feature, "ff_definition") or None,
                    )
                )
                if entry_type is None:
                    entry_type = infer_entry_type(feature.qualifiers, submission_category, division)
                if completeness is None:
                    completeness = infer_completeness(_text_value(feature, "ff_definition"))
            elif feature.type == "TOPOLOGY":
                if "circular" in feature.qualifiers:
                    topology = "circular"
            elif feature.type == "COMMENT":
                continue  # build_comments でまとめて読む
            else:
                features.append(_build_feature(feature, entry.id))

        entries.append(
            EntryV3(
                alias=entry.id,
                name=entry.name,
                type=entry_type if entry_type is not None else infer_entry_type({}, submission_category, division),
                topology=topology,
                division=division,
                completeness=completeness,
                sequence=sequence_dict[entry.id],
                comments=build_comments(entry.features) or None,
                source_features=source_features or None,
            )
        )

    sequences = Sequences(
        seq_prefix=get_seq_prefix(mss, submission_category),
        common_source=_source_from_qualifiers(common_qualifiers),
        entries=entries or None,
        structured_comments=build_structured_comments(common_entry) or None,
        keywords=build_keywords(common_entry) or None,
    )

    return DdbjRecordV3(
        schema_version="v3",
        provenance=Provenance.model_validate(provenance_data),
        submission=submission,
        projects=[project] if project is not None else None,
        sequences=sequences,
        features=features or None,
        relations=build_relations(common_entry) or None,
    )


def ann2json_for_dfast(
    ann_file: Path,
    seq_file: Path,
    out_json_file: Path | None = None,
    division: str | None = None,
    record_version: RecordVersion = "v2",
) -> None:
    """
    DFAST が生成する ann ファイルと seq ファイルを読み込み、json形式のファイルに変換する

    record_version が v3 のときは MSS から v3 を直接組み立てる。v1 / v2 のときは v1 を組み立て、
    v2 は ddbj_record の converter で v1 から変換する

    todo: CDS featureのtranslationが含まれていないので, jsonにtranslationを含めるか、jsonの情報を元にtranslationを生成する機能が必要
    """
    # MSSインスタンスを作成してパース
    mss = MSS.parse(ann_file, seq_file)

    ddbj_record: DdBjRecordV1 | DdBjRecordV2 | DdbjRecordV3
    if record_version == "v3":
        ddbj_record = mss_to_record_v3(mss, division=division)
    else:
        dat: dict[str, Any] = {"schema_version": DATA_MODEL_VERSION}

        common_dict = common_to_dict(mss)
        common_source = source_to_dict(mss)
        common_meta, seq_info = infer_meta_from_ann(mss)
        trad_submission_category = common_meta["COMMON_META"].pop("trad_submission_category")
        common_dict["COMMON"]["trad_submission_category"] = trad_submission_category
        if division:
            # divisionが指定されていれば変更 (デフォルトでは UNK)
            common_meta["COMMON_META"]["division"] = division

        dict_entry = make_entry_dict(mss, seq_info)

        dat.update(common_dict)
        dat.update(common_source)
        dat.update(common_meta)
        dat.update(dict_entry)

        if record_version == "v1":
            ddbj_record = DdBjRecordV1.model_validate(dat)
        else:
            v1_record = DdBjRecordV1.model_validate(dat)
            ddbj_record = v1_to_v2(v1_record)

    json_str = ddbj_record.model_dump_json(indent=4, exclude_none=True, by_alias=True)
    if out_json_file is None:
        print(json_str)
    else:
        with open(out_json_file, "w", encoding="utf-8") as f:
            f.write(json_str)


def main() -> None:
    logging.basicConfig(level=logging.WARNING)

    # 引数のパース
    argparser = argparse.ArgumentParser(description="Convert MSS annotation and sequence files to DFAST json file")
    argparser.add_argument("ann_file", type=str, help="MSS Annotation file")
    argparser.add_argument("seq_file", type=str, help="MSS Sequence file")
    argparser.add_argument(
        "-o",
        "--out_json_file",
        type=str,
        help="Output json file. If not specified, output will written to stdout.",
        default=None,
    )
    argparser.add_argument(
        "-d", "--division", type=str, help="Division. Default is BCT for DFAST result or UNK", default=None
    )
    argparser.add_argument(
        "-r",
        "--record_version",
        type=str,
        choices=["v1", "v2", "v3"],
        help="DDBJ Record version. Default is v2",
        default="v2",
    )

    if len(sys.argv) == 1:
        argparser.print_help()
        sys.exit(1)

    args = argparser.parse_args()

    ann_file = args.ann_file
    seq_file = args.seq_file
    out_json_file = args.out_json_file

    logger.info("ann_file: %s, seq_file: %s, out_json_file: %s", ann_file, seq_file, out_json_file)
    ann2json_for_dfast(ann_file, seq_file, out_json_file, division=args.division, record_version=args.record_version)


if __name__ == "__main__":
    main()
