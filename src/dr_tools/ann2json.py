#!/usr/bin/env python3

import json
import re
from pathlib import Path
from dr_tools.MSS import MSS

DATA_MODEL_VERSION = "0.1"

def common_to_dict(mss: MSS) -> dict:
    # COMMON entryの情報を取得して辞書で返す
    qualifiers_with_multi_values = ["line", "ab_name", "keyword", "sequence read archive"]
    features_with_multi_values = ["COMMENT", "REFERENCE"]
    common_entry = [entry for entry in mss.entries if entry.id == "COMMON"]
    len_common_entry = len(common_entry)
    if len_common_entry == 0:
        raise ValueError("COMMON entry not found.")
    elif len_common_entry > 1:
        raise ValueError("Multiple COMMON entries found.")
    common_entry = common_entry[0]
    D = {}
    for feature in common_entry.features:
        inner_dict = {}
        for key, value in feature.qualifiers.items():
            if key in qualifiers_with_multi_values:
                inner_dict[key] = value
            else:
                assert len(value) == 1
                inner_dict[key] = value[0]
        if feature.type in features_with_multi_values:
            D.setdefault(feature.type, []).append(inner_dict)
        else:
            D[feature.type] = inner_dict
    return {"COMMON": D}

def source_to_dict(mss):
    """
    source featureの情報を取得して、共通する項目をCOMMON_SOURCEとして辞書で返す
    COMMON_SOURCEに含まれるqualifierは各entryのsource featureからは削除する
    """
    qualifiers_with_multi_values = ["note"]
    ignore_qualifiers = ["submitter_seqid", "ff_definition"]
    source_features = [feature for entry in mss.entries for feature in entry.features if feature.type == "source"]

    # すべてのenrtyについて等しい値を持つqualifierを抽出する
    common_qualifiers = {}
    for feature in source_features:
        if feature.type != "source":
            continue
        for key, value in feature.qualifiers.items():
            if key not in ignore_qualifiers:
                value = "//".join(value)  # 一時的にリストを文字列に変換
                common_qualifiers.setdefault(key, []).append(value)
    common_qualifiers = {key: value for key, value in common_qualifiers.items() if len(value) == len(source_features) and len(set(value))==1}

    # 各entryのsource featureからCOMMON_SOURCEに含まれるqualifierを削除
    common_qualifiers_keys = set(common_qualifiers.keys())
    for source_feature in source_features:
        for key in common_qualifiers_keys:
            if key in source_feature.qualifiers:
                del source_feature.qualifiers[key]

    # COMMON_SOURCEとして返す辞書を作成
    D = {}
    for key, value in common_qualifiers.items():
        value = value[0]
        if key in qualifiers_with_multi_values:
            D[key] = value.split("//")
        else:
            D[key] = value
    return {"COMMON_SOURCE": D}

def infer_meta_from_ann(mss):
    """
    annファイルに記載された情報から、
    trad_submission_category: GNM; WGS
    seq_type: chromosome; plasmid; unplaced; other
    seq_topology: linear; circular
    を決定する
    common_meta と seq_info の2つの辞書を返す
    注意: DFASTが出力するannファイルにのみ対応
    """

    def _get_locus_tag_prefix(entries) -> str|None:
        """
        locus_tag_prefix の決定
        entry, featureのqualifierにlocus_tagが含まれる場合、その値の _ より前の部分をlocus_tag_prefixとする
        """
        for entry in entries:
            for feature in entry.features:
                if "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    if locus_tag:
                        return locus_tag.split("_")[0]
        return None

    def _get_dfast_version(entries) -> str|None:
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
                            m = dfast_version_pat.search(qualifier_value)
                            if m:
                                dfast_version = m.group(2)
                                # print(f"dfast_version: {dfast_version}, match1 '{m.group(1)}', match2 '{m.group(2)}'")  # debug
                                return dfast_version
        return None

    dfast_version = _get_dfast_version(mss.entries)
    if dfast_version:
        common_meta = {"dfast_version": dfast_version, "division": "BCT"} # For DFAST, the division is always BCT
    else:
        common_meta = {"division": "UNK"}  

    # trad_submission_category (WGS, GNM) かの判定
    # COMMON entryのDATATYPE featureのtypeがWGSであればWGSとする
    # それ以外はGNMとする 
    trad_submission_category = "GNM"
    for entry in mss.entries:
        if entry.id == "COMMON":
            for feature in entry.features:
                if feature.type == "DATATYPE" and "WGS" in feature.qualifiers.get("type", []):
                    trad_submission_category = "WGS"
                    break
            break
    common_meta["trad_submission_category"] = trad_submission_category
    # print(f"_trad_submission_category: {trad_submission_category}")  # debug

    # seq_prefixの決定 (wgsの場合のみ)
    # sequence01 contig01 等から数字の部分を取り除いたもの。
    seq_prefix = ""
    if trad_submission_category == "WGS":
        for entry in mss.entries:
            if entry.id == "COMMON":
                continue
            seq_name = entry.name
            seq_prefix = seq_name.rstrip("0123456789")
            break
    if seq_prefix:
        common_meta["seq_prefix"] = seq_prefix
    # print(f"seq_prefix: {seq_prefix}")  # debug

    # locus_tag_prefix の決定
    locus_tag_prefix = _get_locus_tag_prefix(mss.entries)
    if locus_tag_prefix:
        common_meta["locus_tag_prefix"] = locus_tag_prefix

    # seq_infoの決定
    # gnm (complete) の場合、 seq_type, seq_topology
    # wgs (draft) の場合、不要
    seq_info = {}
    if trad_submission_category == "GNM":
        for entry in mss.entries:
            source_or_tpoology = [feature for feature in entry.features if feature.type in ["source", "TOPOLOGY"]]
            seq_topology, seq_type = "linear", "other"
            for feature in source_or_tpoology:
                if feature.type == "source":
                    ff_definition = feature.qualifiers.get("ff_definition", [""])[0]
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

def make_entry_dict(mss, seq_info):
    """
    mssオブジェクトを受け取り、各entryの情報を辞書で返す
    """
    # D = {}
    L = []
    sequence_dict = {sequence.id: sequence.seq for sequence in mss.sequences}
    for entry in mss.entries:
        if entry.id == "COMMON":
            continue
        sequence = sequence_dict[entry.id]  # note: sequence.id and entry.id is always identical, but entry.name may be different.
        # sequence = sequence[:10] + "..."  # for debug
        seq_info_dict = seq_info.get(entry.id, {})
        seq_type, seq_topology = seq_info_dict.get("seq_type", "other"), seq_info_dict.get("seq_topology", "linear")
        features = [feature.to_dict() for feature in entry.features if feature.type != "TOPOLOGY"]
        # features = {feature.id: feature.to_dict() for feature in entry.features if feature.type != "TOPOLOGY"}
        # D[entry.id] = {"name": entry.name, "_type": seq_type, "_topology": seq_topology, "sequence": sequence, "features": features}
        L.append({"id": entry.id, "name": entry.name, "type": seq_type, "topology": seq_topology, "sequence": sequence, "features": features})
    return {"ENTRIES": L}


def ann2json_for_dfast(ann_file: Path, seq_file: Path, out_json_file: Path=None, division: str|None=None) -> None:
    """
    DFAST が生成する ann ファイルと seq ファイルを読み込み、json形式のファイルに変換する

    todo: CDS featureのtranslationが含まれていないので, jsonにtranslationを含めるか、jsonの情報を元にtranslationを生成する機能が必要
    """
    dat = {"schema_version": DATA_MODEL_VERSION}

    # MSSインスタンスを作成してパース
    mss = MSS.parse(ann_file, seq_file)
    

    common_dict = (common_to_dict(mss))

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

    if out_json_file is None:
        print(json.dumps(dat, indent=4))
    else:
        with open(out_json_file, "w") as f:
            json.dump(dat, f, indent=4)

def main():
    import sys
    import logging
    import argparse
    logging.basicConfig(level=logging.WARNING)

    # 引数のパース
    argparser = argparse.ArgumentParser(description='Convert MSS annotation and sequence files to DFAST json file')
    argparser.add_argument('ann_file', type=str, help='MSS Annotation file')
    argparser.add_argument('seq_file', type=str, help='MSS Sequence file')
    argparser.add_argument('-o', '--out_json_file', type=str, help='Output json file. If not specified, output will written to stdout.', default=None)
    argparser.add_argument('-d', '--division', type=str, help='Division. Default is BCT for DFAST result or UNK', default=None)

    if len(sys.argv) == 1:
        argparser.print_help()
        sys.exit(1)

    args = argparser.parse_args()

    ann_file = args.ann_file
    seq_file = args.seq_file
    out_json_file = args.out_json_file

    logging.info(f"ann_file: {ann_file}, seq_file: {seq_file}, out_json_file: {out_json_file}")
    ann2json_for_dfast(ann_file, seq_file, out_json_file, division=args.division)

if __name__ == "__main__":
    main()