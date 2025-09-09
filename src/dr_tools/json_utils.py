# This file contains utility functions for working with JSON files
import json
from pathlib import Path
from typing import Any, Dict, Literal, Tuple

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from ddbj_record.converter.v1_to_v2 import v1_to_v2
from ddbj_record.converter.v2_to_v1 import v2_to_v1
from ddbj_record.schema.v1 import DdbjRecord as DdbjRecordV1
from ddbj_record.schema.v2 import DdbjRecord as DdbjRecordV2


def load_json_to_ddbj_record_instance(
    json_file: Path,
    to_record_version: Literal["v1", "v2"] = "v1"
) -> DdbjRecordV1 | DdbjRecordV2:
    """
    JSONファイルを読み込み、DdbjRecordインスタンスを生成して返す
    """
    with open(json_file, encoding="utf-8") as f:
        raw_data = json.load(f)

    if raw_data["schema_version"] in ("0.1", "v1"):
        record_v1_instance = DdbjRecordV1.model_validate(raw_data)
        if to_record_version == "v1":
            return record_v1_instance
        else:
            record_v2_instance = v1_to_v2(record_v1_instance)
            return record_v2_instance
    elif raw_data["schema_version"] in ("0.2", "v2"):
        record_v2_instance = DdbjRecordV2.model_validate(raw_data)
        if to_record_version == "v2":
            return record_v2_instance
        else:
            record_v1_instance = v2_to_v1(record_v2_instance)
            return record_v1_instance
    else:
        raise ValueError(f"Unsupported schema_version: {raw_data['schema_version']}")


def get_feature_and_entry_json(json_dat: Dict[str, Any], feature_id: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    feature_idを指定しmJSONデータからfeature_idに対応するfeatureとentryのjsonデータを取得する
    """
    features = [(feature, entry) for entry in json_dat.get("ENTRIES", [])
                for feature in entry.get("features", []) if feature.get("id", "unknown") == feature_id]
    if len(features) == 0:
        raise ValueError(f"Feature with id {feature_id} not found in JSON data")
    elif len(features) > 1:
        raise ValueError(f"Multiple features with id {feature_id} found in JSON data")
    feature_json, entry_json = features[0]
    return feature_json, entry_json


def get_feature_json(json_file: Path, feature_id: str) -> Dict[str, Any]:
    """
    JSONデータからfeature_idに対応するfeatureの情報を辞書として取得する
    DFAST web serviceの feature 詳細で表示される情報を取得する
    """
    # === ddbj record v2 対応 ===
    ddbj_record_instance = load_json_to_ddbj_record_instance(json_file, to_record_version="v1")
    json_dat = ddbj_record_instance.model_dump(exclude_none=True, by_alias=True)  # dict形式に変換
    feature_json, entry_json = get_feature_and_entry_json(json_dat, feature_id)
    locus_tag_prefix = get_locus_tag_prefix(json_dat)
    set_locus_tag(feature_json, locus_tag_prefix)

    # SeqFeatureオブジェクトを作成 (CDSの場合は、qualifiersにtranslationを追加)
    seq_feature, nucleotide = json_to_seqfeature(feature_json, entry_json)
    feature_json["nucleotide"] = str(nucleotide)
    translate = seq_feature.qualifiers.get("translation", [""])[0]
    feature_json["translation"] = translate
    return feature_json


def get_locus_tag_prefix(json_dat: Dict[str, Any]) -> str:
    """
    JSONデータからlocus_tagのprefixを取得する
    JSONに定義されていない場合は、デフォルト値"LOCUS"を返す
    """
    return json_dat.get("COMMON_META", {}).get("locus_tag_prefix", "LOCUS")  # type: ignore


def set_locus_tag(feature_json: Dict[str, Any], locus_tag_prefix: str) -> None:
    """
    feature_jsonにlocus_tagを設定する
    locus_tag_prefixとfeatureのlocus_tag_idを結合した文字列をlocus_tagとして設定する
    locus_tag_id を持たない場合は、locus_tagを設定しない
    """
    if "locus_tag_id" in feature_json:
        locus_tag = locus_tag_prefix + "_" + feature_json["locus_tag_id"]
        feature_json["qualifiers"]["locus_tag"] = [locus_tag]


def json_to_seqfeature(feature_json: Dict[str, Any], entry_json: Dict[str, Any]) -> SeqFeature:
    """
    JSON で書かれた featureを BioPython SeqFeature object に変換し、
    そのオブジェクトと、その feature location が指す配列を Seq objectとして返す
    """
    from dr_tools.json2biopython import (add_translate_qualifier,
                                         create_seqfeature)

    seq = Seq(entry_json.get("sequence", ""))
    topology = entry_json.get("topology", "linear")
    seq_length = entry_json.get("length") or len(entry_json.get("sequence", ""))

    feature_type = feature_json["type"]
    feature_location_str = feature_json.get("location", "")
    qualifiers = feature_json.get("qualifiers", {})

    feature = create_seqfeature(feature_type, feature_location_str, qualifiers, seq_length, topology)
    if feature_type == "CDS":
        add_translate_qualifier(feature, seq)
    nucleotide = feature.extract(seq)
    return feature, nucleotide
