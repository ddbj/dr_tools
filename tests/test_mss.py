import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest
from hypothesis import given
from hypothesis import strategies as st

from dr_tools.MSS import MSS, Entry, Feature, Sequence

# === helpers for MSS.from_json tests ===
# from_json は ddbj_record スキーマでバリデーションされた JSON を読むので、
# 契約を確かめるのに必要な最小限のフィールドだけを埋めた record を組み立てる。


def _minimal_record(entries: List[Dict[str, Any]], locus_tag_prefix: Optional[str] = None) -> Dict[str, Any]:
    common_meta: Dict[str, Any] = {"division": "BCT"}
    if locus_tag_prefix is not None:
        common_meta["locus_tag_prefix"] = locus_tag_prefix
    return {
        "schema_version": "v1",
        "COMMON": {
            "SUBMITTER": {
                "contact": "Test Contact",
                "email": "test@example.com",
                "institute": "Test Institute",
                "country": "Japan",
                "city": "Tokyo",
                "street": "1-2-3",
                "zip": "100-0001",
            },
            "trad_submission_category": "GNM",
        },
        "COMMON_SOURCE": {
            "organism": "Escherichia coli",
            "mol_type": "genomic DNA",
        },
        "COMMON_META": common_meta,
        "ENTRIES": entries,
    }


def _write_record(tmp_path: Path, record: Dict[str, Any]) -> Path:
    json_file = tmp_path.joinpath("record.json")
    json_file.write_text(json.dumps(record), encoding="utf-8")
    return json_file


# === Sequence tests ===


def test_sequence_to_fasta_default() -> None:
    seq = Sequence(id="seq1", seq="ATGC" * 20)
    fasta = seq.to_fasta()
    lines = fasta.strip().split("\n")
    assert lines[0] == ">seq1"
    assert all(len(line) <= 60 for line in lines[1:])


def test_sequence_to_fasta_custom_width() -> None:
    seq = Sequence(id="seq2", seq="ATGC" * 10)
    fasta = seq.to_fasta(width=4)
    lines = fasta.strip().split("\n")
    assert lines[0] == ">seq2"
    assert all(len(line) == 4 for line in lines[1:])


def test_sequence_to_fasta_separator() -> None:
    seq = Sequence(id="seq3", seq="ATGC" * 5)
    fasta = seq.to_fasta(separator=True)
    assert fasta.endswith("//\n")


# === Feature.to_dict tests ===


def test_feature_to_dict_valid_locus_tag() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"locus_tag": ["abc_123"], "note": ["test"]})
    d = f.to_dict()
    assert d["locus_tag_id"] == "123"
    assert "locus_tag" not in d["qualifiers"]


def test_feature_to_dict_invalid_locus_tag() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"locus_tag": ["abc123"]})
    with pytest.raises(ValueError):
        f.to_dict()


# === Feature.to_tsv tests ===


def test_feature_to_tsv_empty_string_value_omits_row() -> None:
    f = Feature(type="CDS", id=1, location="1..10",
                qualifiers={"note": [""], "product": ["hypothetical protein"]})
    tsv = f.to_tsv()
    assert ["", "", "", "note", ""] not in tsv
    assert [row for row in tsv if row[3] == "note"] == []
    assert [row for row in tsv if row[3] == "product"] == [["", "CDS", "1..10", "product", "hypothetical protein"]]


def test_feature_to_tsv_bool_qualifier_yields_row_with_empty_value_column() -> None:
    f = Feature(type="TOPOLOGY", id=1, qualifiers={"circular": [True]})
    tsv = f.to_tsv()
    assert tsv == [["", "TOPOLOGY", "", "circular", ""]]


def test_feature_to_tsv_same_key_mixed_empty_and_nonempty_keeps_only_nonempty() -> None:
    f = Feature(type="gene", id=1, location="1..5", qualifiers={"note": ["", "real note", ""]})
    tsv = f.to_tsv()
    assert tsv == [["", "gene", "1..5", "note", "real note"]]


def test_feature_to_tsv_no_qualifiers_yields_single_placeholder_row() -> None:
    f = Feature(type="UTR", id=1, location="1..50", qualifiers={})
    tsv = f.to_tsv()
    assert tsv == [["", "UTR", "1..50", "", ""]]


def test_feature_to_tsv_only_empty_string_qualifiers_yields_placeholder_row_without_indexerror() -> None:
    # bool でない qualifier が全て空文字列だと行が1つも積まれないので、
    # placeholder を足す分岐を通らないと tsv[0] で IndexError になる
    f = Feature(type="gene", id=1, location="1..5", qualifiers={"note": [""], "product": [""]})
    tsv = f.to_tsv()
    assert tsv == [["", "gene", "1..5", "", ""]]


def test_feature_to_tsv_first_row_second_column_is_always_type() -> None:
    f = Feature(type="mRNA", id=1, location="10..20", qualifiers={"product": ["x"], "note": ["y"]})
    tsv = f.to_tsv()
    assert tsv[0][1] == "mRNA"
    assert tsv[0] == ["", "mRNA", "10..20", "product", "x"]
    assert tsv[1] == ["", "", "", "note", "y"]


def test_feature_to_tsv_location_set_on_first_row_when_present() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"note": ["test"]})
    tsv = f.to_tsv()
    assert tsv[0] == ["", "gene", "1..100", "note", "test"]


def test_feature_to_tsv_no_location_leaves_third_column_empty() -> None:
    f = Feature(type="DBLINK", id=1, location="", qualifiers={"project": ["PRJDB1"]})
    tsv = f.to_tsv()
    assert tsv[0] == ["", "DBLINK", "", "project", "PRJDB1"]


_qualifier_value_strategy = st.one_of(
    st.just(""),
    st.just(" "),
    st.just("　"),  # 全角スペース
    st.just("日本語"),
    st.text(min_size=1, max_size=10),
)


@given(values=st.lists(_qualifier_value_strategy, min_size=0, max_size=6))
def test_feature_to_tsv_property_rows_match_nonempty_values_in_order(values: List[str]) -> None:
    qualifiers: Dict[str, List[str | bool]] = {"note": list(values)} if values else {}
    f = Feature(type="gene", id=1, location="1..5", qualifiers=qualifiers)
    tsv = f.to_tsv()
    expected = [v for v in values if v != ""]
    if expected:
        assert tsv == [["", "gene", "1..5", "note", v] if i == 0 else ["", "", "", "note", v]
                       for i, v in enumerate(expected)]
    else:
        assert tsv == [["", "gene", "1..5", "", ""]]


# === Entry.to_tsv tests ===


def test_entry_to_tsv_first_column_uses_name_not_id() -> None:
    f1 = Feature(type="gene", id=1, location="1..100", qualifiers={"note": ["test"]})
    entry = Entry(id="internal-id-1", name="user-visible-name", features=[f1])
    tsv = entry.to_tsv()
    assert tsv[0][0] == "user-visible-name"
    assert tsv[0][0] != "internal-id-1"


def test_entry_to_tsv_empty_features_returns_empty_list() -> None:
    entry = Entry(id="entry2", name="entry2", features=[])
    assert entry.to_tsv() == []


def test_entry_requires_name_argument() -> None:
    with pytest.raises(TypeError):
        Entry(id="only-id")  # type: ignore[call-arg]


# === MSS.parse / round-trip smoke tests ===


def test_mss_parse_with_example_files(eg_complete: Dict[str, Path], eg_vrl: Dict[str, Path]) -> None:
    MSS.parse(eg_complete["ann"], eg_complete["seq"])
    MSS.parse(eg_vrl["ann"], eg_vrl["seq"])


def test_mss_to_tsv(mss_complete: MSS, mss_vrl: MSS) -> None:
    mss_complete.to_tsv()
    mss_vrl.to_tsv()


def test_mss_to_fasta(mss_complete: MSS, mss_vrl: MSS) -> None:
    mss_complete.to_fasta()
    mss_vrl.to_fasta()


# === MSS.from_json tests ===


def test_from_json_smoke_examples(eg_complete: Dict[str, Path], eg_vrl: Dict[str, Path]) -> None:
    MSS.from_json(eg_complete["v1_json"])
    MSS.from_json(eg_complete["v2_json"])
    MSS.from_json(eg_vrl["v1_json"])
    MSS.from_json(eg_vrl["v2_json"])


def test_from_json_circular_topology_adds_topology_feature_with_bool_circular_qualifier(tmp_path: Path) -> None:
    entries = [{
        "id": "e1", "name": "e1", "type": "chromosome", "topology": "circular",
        "sequence": "ATGCATGCATGC",
        "features": [
            {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
        ],
    }]
    json_file = _write_record(tmp_path, _minimal_record(entries))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    topology_features = [f for f in entry.features if f.type == "TOPOLOGY"]
    assert len(topology_features) == 1
    assert topology_features[0].qualifiers == {"circular": [True]}


def test_from_json_linear_topology_does_not_add_topology_feature(tmp_path: Path) -> None:
    entries = [{
        "id": "e1", "name": "e1", "type": "chromosome", "topology": "linear",
        "sequence": "ATGCATGCATGC",
        "features": [
            {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
        ],
    }]
    json_file = _write_record(tmp_path, _minimal_record(entries))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    assert [f.type for f in entry.features if f.type == "TOPOLOGY"] == []


def test_from_json_sequence_id_uses_entry_name_not_id(tmp_path: Path) -> None:
    entries = [{
        "id": "internal-entry-id", "name": "displayed-name", "type": "chromosome", "topology": "linear",
        "sequence": "ATGCATGCATGC",
        "features": [
            {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
        ],
    }]
    json_file = _write_record(tmp_path, _minimal_record(entries))
    mss = MSS.from_json(json_file)

    seq_ids = [s.id for s in mss.sequences]
    assert seq_ids == ["displayed-name"]
    assert "internal-entry-id" not in seq_ids


def test_from_json_builds_locus_tag_from_prefix_and_locus_tag_id(tmp_path: Path) -> None:
    entries = [{
        "id": "e1", "name": "e1", "type": "chromosome", "topology": "linear",
        "sequence": "ATGCATGCATGC",
        "features": [
            {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
            {"id": "f2", "type": "CDS", "location": "1..12", "qualifiers": {"product": ["hypothetical protein"]},
             "locus_tag_id": "00010"},
        ],
    }]
    json_file = _write_record(tmp_path, _minimal_record(entries, locus_tag_prefix="PLH"))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    cds = next(f for f in entry.features if f.type == "CDS")
    assert cds.qualifiers["locus_tag"] == ["PLH_00010"]


def test_from_json_default_locus_tag_prefix_is_locus_when_missing_from_common_meta(tmp_path: Path) -> None:
    entries = [{
        "id": "e1", "name": "e1", "type": "chromosome", "topology": "linear",
        "sequence": "ATGCATGCATGC",
        "features": [
            {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
            {"id": "f2", "type": "CDS", "location": "1..12", "qualifiers": {"product": ["hypothetical protein"]},
             "locus_tag_id": "00020"},
        ],
    }]
    json_file = _write_record(tmp_path, _minimal_record(entries, locus_tag_prefix=None))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    cds = next(f for f in entry.features if f.type == "CDS")
    assert cds.qualifiers["locus_tag"] == ["LOCUS_00020"]
