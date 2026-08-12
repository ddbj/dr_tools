import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

from ddbj_record.schema.v1 import Common, CommonMeta, CommonSource
from ddbj_record.schema.v2 import DdbjRecord as V2DdbjRecord

from dr_tools.ann2json import (ann2json_for_dfast, common_to_dict,
                               infer_meta_from_ann, make_entry_dict,
                               source_to_dict)
from dr_tools.json2ann import json2ann
from dr_tools.MSS import Entry, Feature, MSS


def test_common_to_dict_complete(mss_complete: MSS) -> None:
    common_dict = common_to_dict(mss_complete)
    assert "COMMON" in common_dict
    common_obj = common_dict["COMMON"]
    common_obj["trad_submission_category"] = "GNM"  # for test
    Common.model_validate(common_obj)


def test_common_to_dict_vrl(mss_vrl: MSS) -> None:
    common_dict = common_to_dict(mss_vrl)
    assert "COMMON" in common_dict
    common_obj = common_dict["COMMON"]
    common_obj["trad_submission_category"] = "WGS"  # for test
    Common.model_validate(common_obj)


def test_source_to_dict_complete(mss_complete: MSS) -> None:
    source_dict = source_to_dict(mss_complete)
    assert "COMMON_SOURCE" in source_dict
    source_obj = source_dict["COMMON_SOURCE"]
    CommonSource.model_validate(source_obj)


def test_source_to_dict_vrl(mss_vrl: MSS) -> None:
    source_dict = source_to_dict(mss_vrl)
    assert "COMMON_SOURCE" in source_dict
    source_obj = source_dict["COMMON_SOURCE"]
    CommonSource.model_validate(source_obj)


def test_infer_meta_from_ann_complete(mss_complete: MSS) -> None:
    common_meta, _seq_info = infer_meta_from_ann(mss_complete)
    assert "COMMON_META" in common_meta
    common_meta_obj = common_meta["COMMON_META"]
    common_meta_obj.pop("trad_submission_category", None)  # not in ANN
    CommonMeta.model_validate(common_meta_obj)


def test_ann2json_for_dfast_complete(eg_complete: Dict[str, Path]) -> None:
    ann2json_for_dfast(
        ann_file=eg_complete["ann"],
        seq_file=eg_complete["seq"],
        out_json_file=Path("/dev/null"),
        record_version="v1",
    )
    ann2json_for_dfast(
        ann_file=eg_complete["ann"],
        seq_file=eg_complete["seq"],
        out_json_file=Path("/dev/null"),
        record_version="v2",
    )


def test_ann2json_for_dfast_vrl(eg_vrl: Dict[str, Path]) -> None:
    ann2json_for_dfast(
        ann_file=eg_vrl["ann"],
        seq_file=eg_vrl["seq"],
        out_json_file=Path("/dev/null"),
        record_version="v1",
    )
    ann2json_for_dfast(
        ann_file=eg_vrl["ann"],
        seq_file=eg_vrl["seq"],
        out_json_file=Path("/dev/null"),
        record_version="v2",
    )


def test_ann2_json_with_entry_comment(eg_complete: Dict[str, Path]) -> None:
    ann_file = eg_complete["ann"].parent.joinpath("complete_genome_with_comment.ann")
    seq_file = eg_complete["seq"]

    with tempfile.NamedTemporaryFile() as tmp:
        tmp_path = Path(tmp.name)

        ann2json_for_dfast(
            ann_file=ann_file,
            seq_file=seq_file,
            out_json_file=tmp_path,
            record_version="v2",
        )

        json_str = tmp_path.read_text(encoding="utf-8")

    record = V2DdbjRecord.model_validate_json(json_str)
    has_entry_comment = False
    for entry in record.sequences.entries:
        if entry.comments:
            has_entry_comment = True
            break

    assert has_entry_comment


def _count_value_lines(ann_text: str) -> int:
    """5列TSVの5列目 (qualifierの値) が非空の行数を数える"""
    count = 0
    for line in ann_text.splitlines():
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) >= 5 and cols[4] != "":
            count += 1
    return count


def _wgs_mss_with_ff_definition(ff_definition: str) -> MSS:
    """
    WGS の COMMON entry と、ff_definition のみを持つ source feature の entry からなる
    最小の MSS を組み立てる。seq_type の ff_definition 判定が GNM に限定されることを
    確認するために使う (examples の WGS 相当データには plasmid/complete genome 等の
    文言を含む ff_definition が無いため、専用に組み立てる)。
    """
    datatype_qualifiers: Dict[str, List[str | bool]] = {"type": ["WGS"]}
    common_entry = Entry(
        id="COMMON",
        name="COMMON",
        features=[Feature(type="DATATYPE", id="feature_1", qualifiers=datatype_qualifiers)],
    )
    source_qualifiers: Dict[str, List[str | bool]] = {"ff_definition": [ff_definition]}
    contig_entry = Entry(
        id="contig01",
        name="contig01",
        features=[Feature(type="source", id="feature_2", location="1..30", qualifiers=source_qualifiers)],
    )
    return MSS(ann_file=None, seq_file=None, entries=[common_entry, contig_entry], sequences=[])


def test_infer_meta_from_ann_wgs_entry_with_topology_feature_sets_seq_topology_circular(
    wgs_mss_files: Tuple[Path, Path],
) -> None:
    ann_file, seq_file = wgs_mss_files
    mss = MSS.parse(ann_file, seq_file)
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig01"]["seq_topology"] == "circular"


def test_infer_meta_from_ann_wgs_entry_without_topology_feature_sets_seq_topology_linear(
    wgs_mss_files: Tuple[Path, Path],
) -> None:
    ann_file, seq_file = wgs_mss_files
    mss = MSS.parse(ann_file, seq_file)
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig02"]["seq_topology"] == "linear"


def test_infer_meta_from_ann_wgs_ff_definition_plasmid_keyword_seq_type_stays_other() -> None:
    mss = _wgs_mss_with_ff_definition("Escherichia coli DNA, plasmid pXYZ")
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig01"]["seq_type"] == "other"


def test_infer_meta_from_ann_wgs_ff_definition_complete_genome_keyword_seq_type_stays_other() -> None:
    mss = _wgs_mss_with_ff_definition("Escherichia coli DNA, complete genome")
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig01"]["seq_type"] == "other"


def test_infer_meta_from_ann_gnm_ff_definition_plasmid_keyword_sets_seq_type_plasmid(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    assert seq_info["pPLH-1"]["seq_type"] == "plasmid"


def test_infer_meta_from_ann_gnm_ff_definition_complete_genome_keyword_sets_seq_type_chromosome(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    assert seq_info["chromosome"]["seq_type"] == "chromosome"


def test_infer_meta_from_ann_common_entry_included_in_seq_info(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    assert "COMMON" in seq_info


def test_make_entry_dict_excludes_common_entry_from_output(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    entry_dict = make_entry_dict(mss_complete, seq_info)
    entry_ids = [entry["id"] for entry in entry_dict["ENTRIES"]]
    assert "COMMON" not in entry_ids


def test_ann2json_json2ann_roundtrip_wgs_topology_line_count_matches_circular_entries(
    wgs_mss_files: Tuple[Path, Path], tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("wgs.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    topology_lines = [line for line in out_ann_text.splitlines() if "\tTOPOLOGY\t" in line]
    # wgs_mss_files では contig01 のみ circular (TOPOLOGY 行を持つ)
    assert len(topology_lines) == 1


def test_ann2json_json2ann_roundtrip_wgs_locus_tag_lines_preserved(
    wgs_mss_files: Tuple[Path, Path], tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("wgs.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    assert "\t\t\tlocus_tag\tTEST_00010" in out_ann_text
    assert "\t\t\tlocus_tag\tTEST_00020" in out_ann_text


def test_ann2json_json2ann_roundtrip_wgs_value_line_count_not_decreased(
    wgs_mss_files: Tuple[Path, Path], tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("wgs.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    orig_count = _count_value_lines(ann_file.read_text(encoding="utf-8"))
    new_count = _count_value_lines(tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8"))
    assert new_count >= orig_count


def test_ann2json_json2ann_roundtrip_vrl_empty_qualifier_value_line_dropped(
    eg_vrl: Dict[str, Path], tmp_path: Path,
) -> None:
    json_file = tmp_path.joinpath("vrl.json")
    ann2json_for_dfast(eg_vrl["ann"], eg_vrl["seq"], json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    # examples/vrl_result.ann の geo_loc_name は値を持たない (bool qualifierでもない) 唯一の qualifier
    assert "geo_loc_name" not in out_ann_text


def test_ann2json_json2ann_roundtrip_vrl_value_line_count_not_decreased(
    eg_vrl: Dict[str, Path], tmp_path: Path,
) -> None:
    json_file = tmp_path.joinpath("vrl.json")
    ann2json_for_dfast(eg_vrl["ann"], eg_vrl["seq"], json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    orig_count = _count_value_lines(eg_vrl["ann"].read_text(encoding="utf-8"))
    new_count = _count_value_lines(tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8"))
    assert new_count >= orig_count
