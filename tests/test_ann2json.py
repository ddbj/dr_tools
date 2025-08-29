from pathlib import Path
from typing import Dict

from ddbj_record.schema.v1 import Common, CommonMeta, CommonSource

from dr_tools.ann2json import (ann2json_for_dfast, common_to_dict,
                               infer_meta_from_ann, source_to_dict)
from dr_tools.MSS import MSS


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
