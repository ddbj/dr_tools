# This file contains utility functions for working with JSON files
import json
from typing import Dict, Tuple
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from mss_tools.json2biopython import create_seqfeature, add_translate_qualifier

def test():
    json_file = "examples/ecoli_example.json"
    feature_id = "feature_9"
    feature_json = get_feature_json(json_file, feature_id)
    print(feature_json)

def get_feature_and_entry_json(json_dat: Dict, feature_id: str) -> Tuple[Dict, Dict]:
    """
    feature_idを指定しｍJSONデータからfeature_idに対応するfeatureとentryのjsonデータを取得する
    """
    features = [(feature, entry) for entry in json_dat.get("ENTRIES", []) for feature in entry.get("features", [])  if feature.get("id", "unknown") == feature_id]
    if len(features) == 0:
        raise ValueError(f"Feature with id {feature_id} not found in JSON data")
    elif len(features) > 1:
        raise ValueError(f"Multiple features with id {feature_id} found in JSON data")
    feature_json, entry_json = features[0]
    return feature_json, entry_json

def get_feature_json(json_file: Path, feature_id: str) -> Dict:
    """
    JSONデータからfeature_idに対応するfeatureの情報を辞書として取得する
    DFAST web serviceの feature 詳細で表示される情報を取得する
    """
    json_dat = json.load(open(json_file))
    feature_json, entry_json = get_feature_and_entry_json(json_dat, feature_id)
    locus_tag_prefix = get_locus_tag_prefix(json_dat)
    set_locus_tag(feature_json, locus_tag_prefix)

    # SeqFeatureオブジェクトを作成 (CDSの場合は、qualifiersにtranslationを追加)
    seq_feature, nucleotide = json_to_seqfeature(feature_json, entry_json)
    feature_json["nucleotide"] = str(nucleotide)
    translate = seq_feature.qualifiers.get("translation", [""])[0]
    feature_json["translation"] = translate
    return feature_json

def get_locus_tag_prefix(json_dat: Dict) -> str:
    """
    JSONデータからlocus_tagのprefixを取得する
    JSONに定義されていない場合は、デフォルト値"LOCUS"を返す
    """
    return json_dat.get("COMMON_META", {}).get("locus_tag_prefix", "LOCUS")

def set_locus_tag(feature_json: Dict, locus_tag_prefix: str) -> None:
    """
    feature_jsonにlocus_tagを設定する
    locus_tag_prefixとfeatureのlocus_tag_idを結合した文字列をlocus_tagとして設定する
    locus_tag_id を持たない場合は、locus_tagを設定しない
    """
    if "locus_tag_id" in feature_json:
        locus_tag = locus_tag_prefix + "_" + feature_json["locus_tag_id"]
        feature_json["qualifiers"]["locus_tag"] = [locus_tag]


def json_to_seqfeature(feature_json: Dict, entry_json: Dict) -> SeqFeature:
    """
    JSON で書かれた featureを BioPython SeqFeature object に変換し、
    そのオブジェクトと、その feature location が指す配列を Seq objectとして返す
    """
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

if __name__ == "__main__":
    test()