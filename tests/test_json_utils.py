import json
import warnings
from pathlib import Path
from typing import Any, Dict, List

import pytest
from ddbj_record.schema.v1 import DdbjRecord as DdbjRecordV1
from ddbj_record.schema.v2 import DdbjRecord as DdbjRecordV2
from ddbj_record.schema.v3 import DdbjRecord as DdbjRecordV3
from ddbj_record.schema.v3 import Feature, Organism, Qualifier, Source
from hypothesis import given
from hypothesis import strategies as st

from dr_tools.json_utils import (features_by_entry, get_entry_id,
                                 get_feature_and_entry, get_feature_json,
                                 get_locus_tag, get_locus_tag_prefix,
                                 load_json_to_ddbj_record_instance,
                                 merge_source_qualifiers,
                                 qualifiers_to_values, values_to_qualifiers)


def _write_json_with_schema_version(
    src: Path, schema_version: Any, tmp_path: Path, filename: str
) -> Path:
    """src の JSON を読み込み、schema_version だけ差し替えて tmp_path に書き出す"""
    raw: Dict[str, Any] = json.loads(src.read_text(encoding="utf-8"))
    raw["schema_version"] = schema_version
    json_file = tmp_path / filename
    json_file.write_text(json.dumps(raw), encoding="utf-8")
    return json_file


def _v3_record(**kwargs: Any) -> DdbjRecordV3:
    return DdbjRecordV3.model_validate({"schema_version": "v3", **kwargs})


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


# --- load_json_to_ddbj_record_instance: v3 ---

def test_load_json_to_ddbj_record_instance_default_to_record_version_is_v3(eg_complete: Dict[str, Path]) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v1_json"])
    assert isinstance(record, DdbjRecordV3)
    assert record.schema_version == "v3"


def test_load_json_to_ddbj_record_instance_v1_to_v3_keeps_entries_and_features(eg_complete: Dict[str, Path]) -> None:
    record_v1 = load_json_to_ddbj_record_instance(eg_complete["v1_json"], to_record_version="v1")
    record_v3 = load_json_to_ddbj_record_instance(eg_complete["v1_json"], to_record_version="v3")
    assert record_v3.sequences is not None and record_v3.sequences.entries is not None
    assert [entry.alias for entry in record_v3.sequences.entries] == [entry.id for entry in record_v1.ENTRIES]
    n_v1_features = sum(1 for entry in record_v1.ENTRIES for feature in entry.features if feature.type != "source")
    assert len(record_v3.features or []) == n_v1_features


def test_load_json_to_ddbj_record_instance_v2_to_v3_returns_v3_instance(eg_complete: Dict[str, Path]) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v2_json"], to_record_version="v3")
    assert isinstance(record, DdbjRecordV3)
    assert record.sequences is not None and record.sequences.common_source is not None
    assert record.sequences.common_source.organism is not None
    assert record.sequences.common_source.organism.name == "Paucilactobacillus hokkaidonensis"


def test_load_json_to_ddbj_record_instance_v3_to_v3_returns_v3_instance_without_conversion(
    eg_complete: Dict[str, Path],
) -> None:
    record = load_json_to_ddbj_record_instance(eg_complete["v3_json"], to_record_version="v3")
    assert isinstance(record, DdbjRecordV3)
    assert record.provenance is not None
    assert record.provenance.source_format == "MSS"


def test_load_json_to_ddbj_record_instance_v3_to_v1_raises_value_error(eg_complete: Dict[str, Path]) -> None:
    with pytest.raises(ValueError, match="v3 record cannot be converted to v1"):
        load_json_to_ddbj_record_instance(eg_complete["v3_json"], to_record_version="v1")


def test_load_json_to_ddbj_record_instance_v3_to_v2_raises_value_error(eg_complete: Dict[str, Path]) -> None:
    with pytest.raises(ValueError, match="v3 record cannot be converted to v2"):
        load_json_to_ddbj_record_instance(eg_complete["v3_json"], to_record_version="v2")


def test_load_json_to_ddbj_record_instance_v3_with_minor_version_is_unsupported(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v3_json"], "v3.0", tmp_path, "v3_0.json")
    with pytest.raises(ValueError, match=r"^Unsupported schema_version: v3\.0$"):
        load_json_to_ddbj_record_instance(json_file)


# --- load_json_to_ddbj_record_instance: 異常系 ---

def test_load_json_to_ddbj_record_instance_unknown_schema_version_raises_value_error_with_original_value(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json_file = _write_json_with_schema_version(eg_complete["v1_json"], "v9", tmp_path, "unknown.json")
    with pytest.raises(ValueError, match=r"^Unsupported schema_version: v9$"):
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


# --- get_locus_tag_prefix / get_locus_tag ---

def test_get_locus_tag_prefix_project_has_prefix_returns_prefix() -> None:
    record = _v3_record(projects=[{"locus_tag_prefix": [{"prefix": "ABCD"}]}])
    assert get_locus_tag_prefix(record) == "ABCD"


def test_get_locus_tag_prefix_no_projects_returns_default() -> None:
    assert get_locus_tag_prefix(_v3_record()) == "LOCUS"


def test_get_locus_tag_prefix_project_without_prefix_returns_default() -> None:
    record = _v3_record(projects=[{"publications": [{"title": "T"}]}])
    assert get_locus_tag_prefix(record) == "LOCUS"


def test_get_locus_tag_prefix_first_project_without_prefix_uses_later_project() -> None:
    record = _v3_record(projects=[{"accession": "PRJDB1"}, {"locus_tag_prefix": [{"prefix": "XYZ"}]}])
    assert get_locus_tag_prefix(record) == "XYZ"


def test_get_locus_tag_prefix_empty_prefix_is_skipped() -> None:
    record = _v3_record(projects=[{"locus_tag_prefix": [{"prefix": ""}, {"prefix": "SECOND"}]}])
    assert get_locus_tag_prefix(record) == "SECOND"


def test_get_locus_tag_feature_has_locus_tag_id_returns_prefix_and_id() -> None:
    assert get_locus_tag(Feature(locus_tag_id="00010"), "PREFIX") == "PREFIX_00010"


def test_get_locus_tag_feature_without_locus_tag_id_returns_none() -> None:
    assert get_locus_tag(Feature(qualifiers={"product": [Qualifier(value="x")]}), "PREFIX") is None


# --- qualifiers_to_values / values_to_qualifiers ---

def test_qualifiers_to_values_qualifier_without_value_becomes_true() -> None:
    values = qualifiers_to_values({"pseudo": [Qualifier()], "product": [Qualifier(value="p")]})
    assert values == {"pseudo": [True], "product": ["p"]}


def test_qualifiers_to_values_none_returns_empty_dict() -> None:
    assert qualifiers_to_values(None) == {}


def test_values_to_qualifiers_true_becomes_qualifier_without_value() -> None:
    qualifiers = values_to_qualifiers({"pseudo": [True]})
    assert qualifiers == {"pseudo": [Qualifier()]}
    assert qualifiers["pseudo"][0].model_dump(exclude_none=True) == {}


def test_values_to_qualifiers_empty_dict_returns_none() -> None:
    assert values_to_qualifiers({}) is None


def test_values_to_qualifiers_empty_string_is_kept_as_value() -> None:
    assert values_to_qualifiers({"note": [""]}) == {"note": [Qualifier(value="")]}


_qualifier_value = st.one_of(st.just(True), st.text(max_size=8))


@given(values=st.dictionaries(st.text(min_size=1, max_size=5), st.lists(_qualifier_value, min_size=1, max_size=3), max_size=4))
def test_qualifier_values_roundtrip_through_v3_qualifiers_is_identity(values: Dict[str, List[str | bool]]) -> None:
    qualifiers = values_to_qualifiers(values)
    assert qualifiers_to_values(qualifiers) == values


# --- merge_source_qualifiers ---

def test_merge_source_qualifiers_entry_qualifier_replaces_common_qualifier_of_same_name() -> None:
    common = Source(qualifiers={"strain": [Qualifier(value="A")], "isolation_source": [Qualifier(value="soil")]})
    entry = Source(qualifiers={"strain": [Qualifier(value="B")]})
    merged = merge_source_qualifiers(common, entry)
    assert merged == {"strain": ["B"], "isolation_source": ["soil"]}


def test_merge_source_qualifiers_entry_only_qualifier_is_added_after_common() -> None:
    common = Source(organism=Organism(name="E. coli"), mol_type="genomic DNA", qualifiers={"strain": [Qualifier(value="A")]})
    entry = Source(qualifiers={"plasmid": [Qualifier(value="pX")]})
    merged = merge_source_qualifiers(common, entry)
    assert list(merged.items()) == [
        ("organism", ["E. coli"]), ("mol_type", ["genomic DNA"]), ("strain", ["A"]), ("plasmid", ["pX"]),
    ]


def test_merge_source_qualifiers_entry_organism_and_mol_type_override_common() -> None:
    common = Source(organism=Organism(name="E. coli"), mol_type="genomic DNA")
    entry = Source(organism=Organism(name="Vibrio cholerae"), mol_type="genomic RNA")
    assert merge_source_qualifiers(common, entry) == {"organism": ["Vibrio cholerae"], "mol_type": ["genomic RNA"]}


def test_merge_source_qualifiers_entry_organism_without_name_does_not_override_common() -> None:
    common = Source(organism=Organism(name="E. coli"))
    entry = Source(organism=Organism(taxonomy_id="562"))
    assert merge_source_qualifiers(common, entry) == {"organism": ["E. coli"]}


def test_merge_source_qualifiers_both_none_returns_empty_dict() -> None:
    assert merge_source_qualifiers(None, None) == {}


def test_merge_source_qualifiers_common_none_uses_entry_source_only() -> None:
    entry = Source(mol_type="genomic DNA", qualifiers={"note": [Qualifier(value="n1"), Qualifier(value="n2")]})
    assert merge_source_qualifiers(None, entry) == {"mol_type": ["genomic DNA"], "note": ["n1", "n2"]}


# --- features_by_entry / get_entry_id ---

def test_features_by_entry_groups_by_sequence_id_keeping_record_order() -> None:
    record = _v3_record(features=[
        {"alias": "f1", "type": "CDS", "sequence_id": "a"},
        {"alias": "f2", "type": "CDS", "sequence_id": "b"},
        {"alias": "f3", "type": "gene", "sequence_id": "a"},
    ])
    grouped = features_by_entry(record)
    assert {key: [feature.alias for feature in features] for key, features in grouped.items()} == {
        "a": ["f1", "f3"], "b": ["f2"],
    }


def test_features_by_entry_feature_without_sequence_id_is_skipped_with_warning() -> None:
    record = _v3_record(features=[{"alias": "orphan", "type": "CDS"}])
    with pytest.warns(UserWarning, match="orphan"):
        grouped = features_by_entry(record)
    assert grouped == {}


def test_features_by_entry_no_features_returns_empty_dict() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert features_by_entry(_v3_record()) == {}


def test_get_entry_id_uses_alias_before_name() -> None:
    record = _v3_record(sequences={"entries": [{"alias": "a1", "name": "n1"}, {"name": "n2"}]})
    assert record.sequences is not None and record.sequences.entries is not None
    assert [get_entry_id(entry) for entry in record.sequences.entries] == ["a1", "n2"]


def test_get_entry_id_entry_without_alias_and_name_raises_value_error() -> None:
    record = _v3_record(sequences={"entries": [{"type": "chromosome"}]})
    assert record.sequences is not None and record.sequences.entries is not None
    with pytest.raises(ValueError, match="neither alias nor name"):
        get_entry_id(record.sequences.entries[0])


# --- get_feature_and_entry ---

def test_get_feature_and_entry_feature_id_found_returns_feature_and_its_entry() -> None:
    record = _v3_record(
        sequences={"entries": [{"alias": "entry_a"}, {"alias": "entry_b"}]},
        features=[{"alias": "feat_1", "type": "CDS", "sequence_id": "entry_a"},
                  {"alias": "feat_2", "type": "gene", "sequence_id": "entry_b"}],
    )
    feature, entry = get_feature_and_entry(record, "feat_2")
    assert feature.type == "gene"
    assert entry.alias == "entry_b"


def test_get_feature_and_entry_feature_id_not_found_raises_value_error() -> None:
    record = _v3_record(features=[{"alias": "feat_1", "type": "CDS", "sequence_id": "a"}])
    with pytest.raises(ValueError, match="unknown_id"):
        get_feature_and_entry(record, "unknown_id")


def test_get_feature_and_entry_duplicate_feature_id_raises_value_error() -> None:
    record = _v3_record(features=[{"alias": "dup", "type": "CDS", "sequence_id": "a"},
                                  {"alias": "dup", "type": "CDS", "sequence_id": "b"}])
    with pytest.raises(ValueError, match="dup"):
        get_feature_and_entry(record, "dup")


def test_get_feature_and_entry_entry_of_feature_missing_raises_value_error() -> None:
    record = _v3_record(features=[{"alias": "feat_1", "type": "CDS", "sequence_id": "nowhere"}])
    with pytest.raises(ValueError, match="nowhere"):
        get_feature_and_entry(record, "feat_1")


# --- get_feature_json ---

def test_get_feature_json_cds_has_locus_tag_nucleotide_and_translation(full_record: Dict[str, Path]) -> None:
    data = get_feature_json(full_record["v3"], "feature_11")
    assert data["alias"] == "feature_11"
    assert data["type"] == "CDS"
    assert data["sequence_id"] == "chromosome"
    assert data["qualifiers"]["locus_tag"] == [{"value": "PLH_00020"}]
    assert data["qualifiers"]["gene"] == [{"value": "dnaN"}]
    assert data["nucleotide"].upper().startswith("ATG")
    assert data["translation"].startswith("M")


def test_get_feature_json_same_result_from_v1_and_v3_json(full_record: Dict[str, Path]) -> None:
    assert get_feature_json(full_record["v1"], "feature_11") == get_feature_json(full_record["v3"], "feature_11")
