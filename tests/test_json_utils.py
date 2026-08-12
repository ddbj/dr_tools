import json
from pathlib import Path
from typing import Any, Dict

import pytest
from ddbj_record.schema.v1 import DdbjRecord as DdbjRecordV1
from ddbj_record.schema.v2 import DdbjRecord as DdbjRecordV2

from dr_tools.json_utils import (get_feature_and_entry_json,
                                 get_locus_tag_prefix,
                                 load_json_to_ddbj_record_instance,
                                 set_locus_tag)


def _write_json_with_schema_version(
    src: Path, schema_version: Any, tmp_path: Path, filename: str
) -> Path:
    """src の JSON を読み込み、schema_version だけ差し替えて tmp_path に書き出す"""
    raw: Dict[str, Any] = json.loads(src.read_text(encoding="utf-8"))
    raw["schema_version"] = schema_version
    json_file = tmp_path / filename
    json_file.write_text(json.dumps(raw), encoding="utf-8")
    return json_file


# --- load_json_to_ddbj_record_instance: 旧来の表記 ---

def test_load_json_to_ddbj_record_instance_v1_legacy_schema_version_to_v1_returns_v1_instance(
    eg_complete: Dict[str, Path],
) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v1_json"], to_record_version="v1")
    assert isinstance(record, DdbjRecordV1)


def test_load_json_to_ddbj_record_instance_v1_legacy_schema_version_to_v2_returns_v2_instance(
    eg_complete: Dict[str, Path],
) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v1_json"], to_record_version="v2")
    assert isinstance(record, DdbjRecordV2)


def test_load_json_to_ddbj_record_instance_v2_legacy_schema_version_to_v2_returns_v2_instance(
    eg_complete: Dict[str, Path],
) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v2_json"], to_record_version="v2")
    assert isinstance(record, DdbjRecordV2)


def test_load_json_to_ddbj_record_instance_v2_legacy_schema_version_to_v1_returns_v1_instance(
    eg_complete: Dict[str, Path],
) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v2_json"], to_record_version="v1")
    assert isinstance(record, DdbjRecordV1)


# --- load_json_to_ddbj_record_instance: 正規化済みの表記 ---

def test_load_json_to_ddbj_record_instance_v1_normalized_schema_version_to_v1_returns_v1_instance(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v1_json"], "v1.0", tmp_path, "v1_0.json")
    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v1")
    assert isinstance(record, DdbjRecordV1)


def test_load_json_to_ddbj_record_instance_v2_normalized_intermediate_minor_version_to_v2_returns_v2_instance(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v2_json"], "v2.1", tmp_path, "v2_1.json")
    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v2")
    assert isinstance(record, DdbjRecordV2)


def test_load_json_to_ddbj_record_instance_v2_normalized_latest_minor_version_to_v1_returns_v1_instance(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v2_json"], "v2.2", tmp_path, "v2_2.json")
    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v1")
    assert isinstance(record, DdbjRecordV1)


# --- load_json_to_ddbj_record_instance: 異常系 ---

def test_load_json_to_ddbj_record_instance_unknown_schema_version_raises_value_error_with_original_value(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v1_json"], "v3", tmp_path, "unknown.json")
    with pytest.raises(ValueError, match=r"^Unsupported schema_version: v3$"):
        load_json_to_ddbj_record_instance(json_file, to_record_version="v1")


def test_load_json_to_ddbj_record_instance_empty_schema_version_raises_value_error(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v1_json"], "", tmp_path, "empty.json")
    with pytest.raises(ValueError, match=r"^Unsupported schema_version: $"):
        load_json_to_ddbj_record_instance(json_file, to_record_version="v1")


def test_load_json_to_ddbj_record_instance_missing_schema_version_key_raises_value_error(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    raw: Dict[str, Any] = json.loads(eg_complete["v1_json"].read_text(encoding="utf-8"))
    del raw["schema_version"]
    json_file = tmp_path / "missing.json"
    json_file.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match=r"^Unsupported schema_version: $"):
        load_json_to_ddbj_record_instance(json_file, to_record_version="v1")


def test_load_json_to_ddbj_record_instance_non_string_schema_version_raises_value_error(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v1_json"], 2, tmp_path, "non_string.json")
    with pytest.raises(ValueError, match=r"^Unsupported schema_version: 2$"):
        load_json_to_ddbj_record_instance(json_file, to_record_version="v1")


# --- get_locus_tag_prefix ---

def test_get_locus_tag_prefix_common_meta_has_prefix_returns_prefix() -> None:
    json_dat: Dict[str, Any] = {"COMMON_META": {"locus_tag_prefix": "ABCD"}}
    assert get_locus_tag_prefix(json_dat) == "ABCD"


def test_get_locus_tag_prefix_missing_common_meta_key_returns_default() -> None:
    json_dat: Dict[str, Any] = {}
    assert get_locus_tag_prefix(json_dat) == "LOCUS"


def test_get_locus_tag_prefix_common_meta_without_prefix_key_returns_default() -> None:
    json_dat: Dict[str, Any] = {"COMMON_META": {}}
    assert get_locus_tag_prefix(json_dat) == "LOCUS"


# --- set_locus_tag ---

def test_set_locus_tag_feature_has_locus_tag_id_sets_qualifier() -> None:
    feature_json: Dict[str, Any] = {"locus_tag_id": "00010", "qualifiers": {}}
    set_locus_tag(feature_json, "PREFIX")
    assert feature_json["qualifiers"]["locus_tag"] == ["PREFIX_00010"]


def test_set_locus_tag_feature_without_locus_tag_id_does_not_add_qualifier() -> None:
    feature_json: Dict[str, Any] = {"qualifiers": {"product": ["hypothetical protein"]}}
    set_locus_tag(feature_json, "PREFIX")
    assert "locus_tag" not in feature_json["qualifiers"]
    assert feature_json["qualifiers"] == {"product": ["hypothetical protein"]}


# --- get_feature_and_entry_json ---

def test_get_feature_and_entry_json_feature_id_found_returns_feature_and_entry() -> None:
    entry_a: Dict[str, Any] = {"id": "entry_a", "features": [{"id": "feat_1", "type": "CDS"}]}
    entry_b: Dict[str, Any] = {"id": "entry_b", "features": [{"id": "feat_2", "type": "gene"}]}
    json_dat: Dict[str, Any] = {"ENTRIES": [entry_a, entry_b]}

    feature_json, entry_json = get_feature_and_entry_json(json_dat, "feat_2")

    assert feature_json == {"id": "feat_2", "type": "gene"}
    assert entry_json is entry_b


def test_get_feature_and_entry_json_feature_id_not_found_raises_value_error() -> None:
    json_dat: Dict[str, Any] = {"ENTRIES": [{"id": "entry_a", "features": [{"id": "feat_1"}]}]}
    with pytest.raises(ValueError, match="unknown_id"):
        get_feature_and_entry_json(json_dat, "unknown_id")


def test_get_feature_and_entry_json_duplicate_feature_id_raises_value_error() -> None:
    json_dat: Dict[str, Any] = {
        "ENTRIES": [
            {"id": "entry_a", "features": [{"id": "dup"}]},
            {"id": "entry_b", "features": [{"id": "dup"}]},
        ]
    }
    with pytest.raises(ValueError, match="dup"):
        get_feature_and_entry_json(json_dat, "dup")
