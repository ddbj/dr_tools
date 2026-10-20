import json
import tempfile
import warnings
from pathlib import Path
from typing import Any

import pytest
from ddbj_record.schema.v3 import DdbjRecord as DdbjRecordV3
from hypothesis import given
from hypothesis import strategies as st

from dr_tools.MSS import (
    MSS,
    Entry,
    Feature,
    Sequence,
    hold_date_to_mss,
    hold_date_to_record,
    reference_status_to_mss,
    reference_status_to_record,
)

# === helpers for MSS.from_json tests ===
# from_json は ddbj_record スキーマでバリデーションされた JSON を読むので、
# 契約を確かめるのに必要な最小限のフィールドだけを埋めた record を組み立てる。


def _minimal_record(entries: list[dict[str, Any]], locus_tag_prefix: str | None = None) -> dict[str, Any]:
    common_meta: dict[str, Any] = {"division": "BCT"}
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


def _write_record(tmp_path: Path, record: dict[str, Any]) -> Path:
    json_file = tmp_path.joinpath("record.json")
    json_file.write_text(json.dumps(record), encoding="utf-8")
    return json_file


def _v3_record(**kwargs: Any) -> DdbjRecordV3:
    """v3 の最小の record。sequences を省くと entry を 1 つ置く"""
    data: dict[str, Any] = {"schema_version": "v3"}
    data.update(kwargs)
    if "sequences" not in data:
        data["sequences"] = {"entries": [{"alias": "e1", "name": "e1", "sequence": "acgt"}]}
    return DdbjRecordV3.model_validate(data)


def _features_of(mss: MSS, entry_id: str) -> list[Feature]:
    return next(entry for entry in mss.entries if entry.id == entry_id).features


def _common_feature(mss: MSS, feature_type: str) -> Feature | None:
    features = [feature for feature in _features_of(mss, "COMMON") if feature.type == feature_type]
    return features[0] if features else None


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
    with pytest.raises(ValueError, match="Invalid locus_tag"):
        f.to_dict()


def test_feature_to_dict_does_not_change_feature_qualifiers() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"locus_tag": ["abc_123"], "note": ["test"]})
    first = f.to_dict()
    second = f.to_dict()
    assert first == second
    assert second["locus_tag_id"] == "123"
    assert f.qualifiers == {"locus_tag": ["abc_123"], "note": ["test"]}


# === Feature.to_tsv tests ===


def test_feature_to_tsv_empty_string_value_omits_row() -> None:
    f = Feature(type="CDS", id=1, location="1..10", qualifiers={"note": [""], "product": ["hypothetical protein"]})
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
def test_feature_to_tsv_property_rows_match_nonempty_values_in_order(values: list[str]) -> None:
    qualifiers: dict[str, list[str | bool]] = {"note": list(values)} if values else {}
    f = Feature(type="gene", id=1, location="1..5", qualifiers=qualifiers)
    tsv = f.to_tsv()
    expected = [v for v in values if v != ""]
    if expected:
        assert tsv == [
            ["", "gene", "1..5", "note", v] if i == 0 else ["", "", "", "note", v] for i, v in enumerate(expected)
        ]
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


def test_mss_parse_with_example_files(eg_complete: dict[str, Path], eg_vrl: dict[str, Path]) -> None:
    MSS.parse(eg_complete["ann"], eg_complete["seq"])
    MSS.parse(eg_vrl["ann"], eg_vrl["seq"])


_PARSE_ANN = "COMMON\tSUBMITTER\t\tab_name\tMishima,H.\nseq1\tsource\t1..4\torganism\tEscherichia coli\n"


def _parse(tmp_path: Path, ann: str, seq: str) -> MSS:
    ann_file = tmp_path.joinpath("test.ann")
    seq_file = tmp_path.joinpath("test.fa")
    ann_file.write_text(ann, encoding="utf-8")
    seq_file.write_text(seq, encoding="utf-8")
    return MSS.parse(ann_file, seq_file)


@pytest.mark.parametrize(
    "seq",
    [
        ">seq1\nACGT\n//\n>seq2\nGG\nCC\n//\n",
        ">seq1\nACGT\n//\n>seq2\nGG\nCC\n//",
        ">seq1\nACGT\n>seq2\nGG\nCC\n",
        ">seq1\nACGT\n>seq2\nGG\nCC",
        ">seq1\r\nACGT\r\n//\r\n>seq2\r\nGG\r\nCC\r\n//\r\n",
        ">seq1\nACGT\n\n//\n\n>seq2\nGG\nCC\n//\n\n",
    ],
)
def test_mss_parse_seq_separator_and_line_endings_do_not_reach_sequence(tmp_path: Path, seq: str) -> None:
    mss = _parse(tmp_path, _PARSE_ANN, seq)
    assert [(sequence.id, sequence.seq) for sequence in mss.sequences] == [("seq1", "ACGT"), ("seq2", "GGCC")]


@given(
    st.lists(st.text(alphabet="ACGTN", min_size=1, max_size=130), min_size=1, max_size=3), st.booleans(), st.booleans()
)
def test_mss_parse_seq_reads_back_what_to_fasta_writes(
    seqs: list[str], separator: bool, trailing_newline: bool
) -> None:
    mss = MSS(ann_file=None, seq_file=None, sequences=[Sequence(id=f"seq{i}", seq=seq) for i, seq in enumerate(seqs)])
    fasta = mss.to_fasta(separator=separator)
    if not trailing_newline:
        fasta = fasta.rstrip("\n")
    with tempfile.TemporaryDirectory() as tmp_dir:
        read_back = _parse(Path(tmp_dir), _PARSE_ANN, fasta)
    assert [(sequence.id, sequence.seq) for sequence in read_back.sequences] == [
        (sequence.id, sequence.seq) for sequence in mss.sequences
    ]


def test_mss_parse_line_without_trailing_value_column_reads_value_less_qualifier(tmp_path: Path) -> None:
    ann = _PARSE_ANN + "\t\t\tenvironmental_sample\n\tTOPOLOGY\t\tcircular\n"
    mss = _parse(tmp_path, ann, ">seq1\nACGT\n//\n")
    source, topology = mss.entries[1].features
    assert source.qualifiers["environmental_sample"] == [True]
    assert topology.qualifiers == {"circular": [True]}


def test_mss_parse_line_with_only_feature_and_location_columns_makes_feature_without_qualifiers(tmp_path: Path) -> None:
    ann = _PARSE_ANN + "\tgap\t5..10\n"
    mss = _parse(tmp_path, ann, ">seq1\nACGT\n//\n")
    gap = mss.entries[1].features[1]
    assert (gap.type, gap.location, gap.qualifiers) == ("gap", "5..10", {})


def test_mss_parse_value_with_tab_keeps_the_tab(tmp_path: Path) -> None:
    ann = _PARSE_ANN + "\t\t\tnote\ta\tb\n"
    mss = _parse(tmp_path, ann, ">seq1\nACGT\n//\n")
    assert mss.entries[1].features[0].qualifiers["note"] == ["a\tb"]


def test_mss_to_tsv(mss_complete: MSS, mss_vrl: MSS) -> None:
    mss_complete.to_tsv()
    mss_vrl.to_tsv()


def test_mss_to_fasta(mss_complete: MSS, mss_vrl: MSS) -> None:
    mss_complete.to_fasta()
    mss_vrl.to_fasta()


# === MSS.from_json tests (v1 の JSON は v3 に変換してから読む) ===


def test_from_json_smoke_examples(eg_complete: dict[str, Path], eg_vrl: dict[str, Path]) -> None:
    for version in ("v1_json", "v2_json", "v3_json"):
        MSS.from_json(eg_complete[version])
        MSS.from_json(eg_vrl[version])


def test_from_json_circular_topology_adds_topology_feature_with_bool_circular_qualifier(tmp_path: Path) -> None:
    entries = [
        {
            "id": "e1",
            "name": "e1",
            "type": "chromosome",
            "topology": "circular",
            "sequence": "ATGCATGCATGC",
            "features": [
                {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
            ],
        }
    ]
    json_file = _write_record(tmp_path, _minimal_record(entries))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    topology_features = [f for f in entry.features if f.type == "TOPOLOGY"]
    assert len(topology_features) == 1
    assert topology_features[0].qualifiers == {"circular": [True]}


def test_from_json_linear_topology_does_not_add_topology_feature(tmp_path: Path) -> None:
    entries = [
        {
            "id": "e1",
            "name": "e1",
            "type": "chromosome",
            "topology": "linear",
            "sequence": "ATGCATGCATGC",
            "features": [
                {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
            ],
        }
    ]
    json_file = _write_record(tmp_path, _minimal_record(entries))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    assert [f.type for f in entry.features if f.type == "TOPOLOGY"] == []


def test_from_json_sequence_id_uses_entry_name_not_id(tmp_path: Path) -> None:
    entries = [
        {
            "id": "internal-entry-id",
            "name": "displayed-name",
            "type": "chromosome",
            "topology": "linear",
            "sequence": "ATGCATGCATGC",
            "features": [
                {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
            ],
        }
    ]
    json_file = _write_record(tmp_path, _minimal_record(entries))
    mss = MSS.from_json(json_file)

    seq_ids = [s.id for s in mss.sequences]
    assert seq_ids == ["displayed-name"]
    assert "internal-entry-id" not in seq_ids


def test_from_json_builds_locus_tag_from_prefix_and_locus_tag_id(tmp_path: Path) -> None:
    entries = [
        {
            "id": "e1",
            "name": "e1",
            "type": "chromosome",
            "topology": "linear",
            "sequence": "ATGCATGCATGC",
            "features": [
                {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
                {
                    "id": "f2",
                    "type": "CDS",
                    "location": "1..12",
                    "qualifiers": {"product": ["hypothetical protein"]},
                    "locus_tag_id": "00010",
                },
            ],
        }
    ]
    json_file = _write_record(tmp_path, _minimal_record(entries, locus_tag_prefix="PLH"))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    cds = next(f for f in entry.features if f.type == "CDS")
    assert cds.qualifiers["locus_tag"] == ["PLH_00010"]


def test_from_json_default_locus_tag_prefix_is_locus_when_missing_from_common_meta(tmp_path: Path) -> None:
    entries = [
        {
            "id": "e1",
            "name": "e1",
            "type": "chromosome",
            "topology": "linear",
            "sequence": "ATGCATGCATGC",
            "features": [
                {"id": "f1", "type": "source", "location": "1..12", "qualifiers": {"organism": ["Escherichia coli"]}},
                {
                    "id": "f2",
                    "type": "CDS",
                    "location": "1..12",
                    "qualifiers": {"product": ["hypothetical protein"]},
                    "locus_tag_id": "00020",
                },
            ],
        }
    ]
    json_file = _write_record(tmp_path, _minimal_record(entries, locus_tag_prefix=None))
    mss = MSS.from_json(json_file)

    entry = next(e for e in mss.entries if e.id == "e1")
    cds = next(f for f in entry.features if f.type == "CDS")
    assert cds.qualifiers["locus_tag"] == ["LOCUS_00020"]


# === MSS.from_record (record v3) tests ===


def test_from_record_value_less_qualifier_becomes_true_and_is_written_without_value() -> None:
    record = _v3_record(
        features=[
            {
                "alias": "f1",
                "type": "CDS",
                "location": "1..4",
                "sequence_id": "e1",
                "qualifiers": {"pseudo": [{}], "product": [{"value": "p"}]},
            }
        ]
    )
    mss = MSS.from_record(record)
    cds = _features_of(mss, "e1")[0]
    assert cds.qualifiers == {"pseudo": [True], "product": ["p"]}
    assert [row[3:] for row in cds.to_tsv()] == [["pseudo", ""], ["product", "p"]]


def test_from_record_entry_source_overrides_common_source_by_qualifier_name_and_keeps_others() -> None:
    record = _v3_record(
        sequences={
            "common_source": {
                "organism": {"name": "E. coli"},
                "mol_type": "genomic DNA",
                "qualifiers": {"strain": [{"value": "A"}], "isolation_source": [{"value": "soil"}]},
            },
            "entries": [
                {
                    "alias": "e1",
                    "name": "e1",
                    "sequence": "acgt",
                    "source_features": [
                        {
                            "location": "1..4",
                            "source": {"qualifiers": {"strain": [{"value": "B"}], "plasmid": [{"value": "pX"}]}},
                            "definition": ["@@[organism]@@ plasmid @@[plasmid]@@"],
                        },
                    ],
                }
            ],
        }
    )
    source = _features_of(MSS.from_record(record), "e1")[0]
    assert source.type == "source"
    assert source.location == "1..4"
    assert source.qualifiers == {
        "organism": ["E. coli"],
        "mol_type": ["genomic DNA"],
        "strain": ["B"],
        "isolation_source": ["soil"],
        "plasmid": ["pX"],
        "ff_definition": ["@@[organism]@@ plasmid @@[plasmid]@@"],
    }


def test_from_record_entry_without_source_override_gets_common_source_only() -> None:
    record = _v3_record(
        sequences={
            "common_source": {"organism": {"name": "E. coli"}, "qualifiers": {"strain": [{"value": "A"}]}},
            "entries": [{"alias": "e1", "name": "e1", "sequence": "acgt", "source_features": [{"location": "1..4"}]}],
        }
    )
    source = _features_of(MSS.from_record(record), "e1")[0]
    assert source.qualifiers == {"organism": ["E. coli"], "strain": ["A"]}


def test_from_record_source_features_topology_comments_then_features_in_order() -> None:
    record = _v3_record(
        sequences={
            "entries": [
                {
                    "alias": "e1",
                    "name": "e1",
                    "topology": "circular",
                    "sequence": "acgt",
                    "comments": ["c1\nc2", "c3"],
                    "source_features": [{"location": "1..4"}],
                }
            ]
        },
        features=[{"alias": "f1", "type": "gene", "location": "1..4", "sequence_id": "e1"}],
    )
    features = _features_of(MSS.from_record(record), "e1")
    assert [f.type for f in features] == ["source", "TOPOLOGY", "COMMENT", "COMMENT", "gene"]
    assert features[2].qualifiers == {"line": ["c1", "c2"]}
    assert features[3].qualifiers == {"line": ["c3"]}


def test_from_record_multiple_structured_comments_become_multiple_st_comment_features() -> None:
    record = _v3_record(
        sequences={
            "entries": [],
            "structured_comments": [
                {"tagset_id": "Genome-Assembly-Data", "fields": {"Assembly Method": "HGAP", "Genome Coverage": "60x"}},
                {"tagset_id": "Other", "fields": {"x": "1"}},
            ],
        }
    )
    st_comments = [f for f in _features_of(MSS.from_record(record), "COMMON") if f.type == "ST_COMMENT"]
    assert [f.qualifiers for f in st_comments] == [
        {"tagset_id": ["Genome-Assembly-Data"], "Assembly Method": ["HGAP"], "Genome Coverage": ["60x"]},
        {"tagset_id": ["Other"], "x": ["1"]},
    ]


def test_from_record_submission_category_wgs_writes_datatype_feature() -> None:
    record = _v3_record(provenance={"submission_category": "WGS"})
    datatype = _common_feature(MSS.from_record(record), "DATATYPE")
    assert datatype is not None
    assert datatype.qualifiers == {"type": ["WGS"]}


def test_from_record_submission_category_gnm_has_no_datatype_feature() -> None:
    record = _v3_record(provenance={"submission_category": "GNM"})
    assert _common_feature(MSS.from_record(record), "DATATYPE") is None


def test_from_record_submission_category_mag_wgs_writes_datatype_feature() -> None:
    record = _v3_record(provenance={"submission_category": "MAG-WGS"})
    datatype = _common_feature(MSS.from_record(record), "DATATYPE")
    assert datatype is not None
    assert datatype.qualifiers == {"type": ["WGS"]}


def test_from_record_submission_category_mag_gnm_has_no_datatype_feature() -> None:
    record = _v3_record(provenance={"submission_category": "MAG-GNM"})
    assert _common_feature(MSS.from_record(record), "DATATYPE") is None


def test_from_record_env_division_writes_division_feature_before_datatype() -> None:
    record = _v3_record(
        provenance={"submission_category": "MAG-WGS"},
        sequences={
            "entries": [
                {"alias": "e1", "name": "e1", "sequence": "acgt", "division": "ENV"},
                {"alias": "e2", "name": "e2", "sequence": "acgt", "division": "ENV"},
            ]
        },
    )
    common = _features_of(MSS.from_record(record), "COMMON")
    assert [f.qualifiers for f in common if f.type == "DIVISION"] == [{"division": ["ENV"]}]
    types = [f.type for f in common]
    assert types.index("DIVISION") < types.index("DATATYPE")


@pytest.mark.parametrize("division", ["BCT", "VRL", "UNK", None])
def test_from_record_taxonomic_division_has_no_division_feature(division: str | None) -> None:
    record = _v3_record(
        sequences={"entries": [{"alias": "e1", "name": "e1", "sequence": "acgt", "division": division}]}
    )
    assert _common_feature(MSS.from_record(record), "DIVISION") is None


def test_from_record_keywords_become_keyword_feature() -> None:
    record = _v3_record(sequences={"entries": [], "keywords": ["WGS", "STANDARD_DRAFT"]})
    keyword = _common_feature(MSS.from_record(record), "KEYWORD")
    assert keyword is not None
    assert keyword.qualifiers == {"keyword": ["WGS", "STANDARD_DRAFT"]}


def test_from_record_xref_relations_become_dblink_grouped_by_key_in_first_appearance_order() -> None:
    record = _v3_record(
        relations=[
            {"type": "xref", "target": {"db": "biosample", "id": "SAMD1"}},
            {"type": "xref", "target": {"db": "insdc.sra", "id": "DRR1"}},
            {"type": "xref", "target": {"db": "bioproject", "id": "PRJDB1"}},
            {"type": "xref", "target": {"db": "insdc.sra", "id": "DRR2"}},
            {"type": "part_of", "target": {"db": "bioproject", "id": "PRJDB9"}},
        ]
    )
    dblink = _common_feature(MSS.from_record(record), "DBLINK")
    assert dblink is not None
    assert list(dblink.qualifiers.items()) == [
        ("biosample", ["SAMD1"]),
        ("sequence read archive", ["DRR1", "DRR2"]),
        ("project", ["PRJDB1"]),
    ]


def test_from_record_xref_to_unknown_db_is_dropped_with_warning() -> None:
    record = _v3_record(relations=[{"type": "xref", "target": {"db": "assembly", "id": "GCA_1"}}])
    with pytest.warns(UserWarning, match="assembly"):
        mss = MSS.from_record(record)
    assert _common_feature(mss, "DBLINK") is None


def test_from_record_no_relations_has_no_dblink_feature() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _common_feature(MSS.from_record(_v3_record()), "DBLINK") is None


def test_from_record_publication_becomes_reference_with_year_from_date_and_status_in_mss_form() -> None:
    record = _v3_record(
        projects=[
            {
                "publications": [
                    {
                        "title": "T",
                        "status": "in-press",
                        "date": "2023-05-01",
                        "journal": "J",
                        "authors": [{"abbreviation": "Smith,J."}, {"name": "No Abbreviation"}],
                        "consortiums": ["C1"],
                    },
                ]
            }
        ]
    )
    reference = _common_feature(MSS.from_record(record), "REFERENCE")
    assert reference is not None
    assert list(reference.qualifiers.items()) == [
        ("title", ["T"]),
        ("ab_name", ["Smith,J."]),
        ("consrtm", ["C1"]),
        ("status", ["In press"]),
        ("year", ["2023"]),
    ]


def test_from_record_publications_of_all_projects_become_references() -> None:
    record = _v3_record(
        projects=[{"accession": "PRJDB1", "publications": [{"title": "A"}]}, {"publications": [{"title": "B"}]}]
    )
    references = [f for f in _features_of(MSS.from_record(record), "COMMON") if f.type == "REFERENCE"]
    assert [f.qualifiers["title"] for f in references] == [["A"], ["B"]]


def test_from_record_hold_date_is_written_as_yyyymmdd() -> None:
    record = _v3_record(submission={"hold_date": "2025-03-31"})
    date = _common_feature(MSS.from_record(record), "DATE")
    assert date is not None
    assert date.qualifiers == {"hold_date": ["20250331"]}


def test_from_record_submission_comments_become_comment_features_with_lines() -> None:
    record = _v3_record(submission={"comments": ["Comment line1\nline2", "Annotated by DFAST"]})
    comments = [f for f in _features_of(MSS.from_record(record), "COMMON") if f.type == "COMMENT"]
    assert [f.qualifiers for f in comments] == [{"line": ["Comment line1", "line2"]}, {"line": ["Annotated by DFAST"]}]


def test_from_record_submitter_contact_is_first_person_and_ab_names_come_from_all() -> None:
    record = _v3_record(
        submission={
            "submitters": [
                {
                    "name": "Hanako Mishima",
                    "abbreviation": "Mishima,H.",
                    "email": "m@example.com",
                    "role": "contact",
                    "organizations": [
                        {
                            "name": "NIG",
                            "url": "https://example.com",
                            "type": "institution",
                            "department": "DDBJ",
                            "address": {
                                "country": "Japan",
                                "state": "Shizuoka",
                                "city": "Mishima",
                                "street": "Yata 1111",
                                "postal_code": "411-8540",
                            },
                        },
                        {"name": "Big Consortium", "type": "consortium"},
                    ],
                },
                {"abbreviation": "Tanizawa,Y.", "role": "submitter"},
                {"name": "No Abbreviation", "role": "submitter"},
            ]
        }
    )
    submitter = _common_feature(MSS.from_record(record), "SUBMITTER")
    assert submitter is not None
    assert list(submitter.qualifiers.items()) == [
        ("ab_name", ["Mishima,H.", "Tanizawa,Y."]),
        ("contact", ["Hanako Mishima"]),
        ("email", ["m@example.com"]),
        ("url", ["https://example.com"]),
        ("institute", ["NIG"]),
        ("department", ["DDBJ"]),
        ("consrtm", ["Big Consortium"]),
        ("country", ["Japan"]),
        ("state", ["Shizuoka"]),
        ("city", ["Mishima"]),
        ("street", ["Yata 1111"]),
        ("zip", ["411-8540"]),
    ]


def test_from_record_submitter_without_any_value_has_no_submitter_feature() -> None:
    record = _v3_record(submission={"submitters": [{"role": "contact"}]})
    assert _common_feature(MSS.from_record(record), "SUBMITTER") is None


def test_from_record_common_feature_order_follows_mss_convention() -> None:
    record = _v3_record(
        provenance={"submission_category": "WGS"},
        submission={"hold_date": "2025-01-01", "comments": ["c"], "submitters": [{"abbreviation": "A,B."}]},
        projects=[{"publications": [{"title": "T"}]}],
        sequences={
            "entries": [],
            "keywords": ["WGS"],
            "structured_comments": [{"tagset_id": "X", "fields": {"a": "1"}}],
        },
        relations=[{"type": "xref", "target": {"db": "bioproject", "id": "PRJDB1"}}],
    )
    assert [f.type for f in _features_of(MSS.from_record(record), "COMMON")] == [
        "DBLINK",
        "SUBMITTER",
        "REFERENCE",
        "COMMENT",
        "ST_COMMENT",
        "DATE",
        "DATATYPE",
        "KEYWORD",
    ]


def test_from_record_feature_of_unknown_entry_is_dropped_with_warning() -> None:
    record = _v3_record(features=[{"alias": "f1", "type": "CDS", "location": "1..4", "sequence_id": "nowhere"}])
    with pytest.warns(UserWarning, match="nowhere"):
        mss = MSS.from_record(record)
    assert _features_of(mss, "e1") == []


def test_from_record_entry_without_alias_uses_name_for_entry_and_sequence() -> None:
    record = _v3_record(sequences={"entries": [{"name": "named", "sequence": "acgt"}]})
    mss = MSS.from_record(record)
    assert [(entry.id, entry.name) for entry in mss.entries[1:]] == [("named", "named")]
    assert [(s.id, s.seq) for s in mss.sequences] == [("named", "acgt")]


def test_from_record_entry_name_is_used_for_sequence_header_and_alias_for_feature_lookup() -> None:
    record = _v3_record(
        sequences={"entries": [{"alias": "internal", "name": "displayed", "sequence": "acgt"}]},
        features=[{"alias": "f1", "type": "gene", "location": "1..4", "sequence_id": "internal"}],
    )
    mss = MSS.from_record(record)
    assert mss.sequences[0].id == "displayed"
    assert [f.type for f in _features_of(mss, "internal")] == ["gene"]
    assert mss.to_tsv().startswith("displayed\tgene\t1..4\t\t\n")


def test_from_record_feature_without_type_raises_value_error() -> None:
    record = _v3_record(features=[{"alias": "notype", "location": "1..4", "sequence_id": "e1"}])
    with pytest.raises(ValueError, match="notype"):
        MSS.from_record(record)


def test_from_record_locus_tag_is_built_from_project_prefix() -> None:
    record = _v3_record(
        projects=[{"locus_tag_prefix": [{"prefix": "ABC"}]}],
        features=[
            {
                "alias": "f1",
                "type": "CDS",
                "location": "1..4",
                "sequence_id": "e1",
                "locus_tag_id": "00010",
                "qualifiers": {"product": [{"value": "p"}]},
            }
        ],
    )
    cds = _features_of(MSS.from_record(record), "e1")[0]
    assert cds.qualifiers == {"product": ["p"], "locus_tag": ["ABC_00010"]}


def test_from_record_without_sequences_has_only_common_entry() -> None:
    record = DdbjRecordV3.model_validate({"schema_version": "v3"})
    mss = MSS.from_record(record)
    assert [entry.id for entry in mss.entries] == ["COMMON"]
    assert mss.sequences == []


# === hold_date / reference status の変換 ===


@pytest.mark.parametrize(
    ("value", "expected"), [("20250331", "2025-03-31"), ("2025-03-31", "2025-03-31"), ("2025", "2025"), ("", "")]
)
def test_hold_date_to_record(value: str, expected: str) -> None:
    assert hold_date_to_record(value) == expected


@pytest.mark.parametrize(
    ("value", "expected"), [("2025-03-31", "20250331"), ("20250331", "20250331"), ("2025-03", "2025-03"), ("", "")]
)
def test_hold_date_to_mss(value: str, expected: str) -> None:
    assert hold_date_to_mss(value) == expected


@given(digits=st.text(alphabet="0123456789", min_size=8, max_size=8))
def test_hold_date_roundtrip_mss_to_record_to_mss_is_identity_for_eight_digits(digits: str) -> None:
    assert hold_date_to_mss(hold_date_to_record(digits)) == digits


@pytest.mark.parametrize(
    ("mss_status", "record_status"),
    [("Unpublished", "unpublished"), ("In press", "in-press"), ("Published", "published")],
)
def test_reference_status_roundtrip(mss_status: str, record_status: str) -> None:
    assert reference_status_to_record(mss_status) == record_status
    assert reference_status_to_mss(record_status) == mss_status
    assert reference_status_to_mss(reference_status_to_record(mss_status)) == mss_status


def test_reference_status_to_mss_is_idempotent_for_mss_form() -> None:
    assert reference_status_to_mss("In press") == "In press"
