from pathlib import Path
from typing import Dict

import pytest

from dr_tools.MSS import MSS, Entry, Feature, Sequence

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


# === Feature tests ===


def test_feature_to_dict_valid_locus_tag() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"locus_tag": ["abc_123"], "note": ["test"]})
    d = f.to_dict()
    assert d["locus_tag_id"] == "123"
    assert "locus_tag" not in d["qualifiers"]


def test_feature_to_dict_invalid_locus_tag() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"locus_tag": ["abc123"]})
    with pytest.raises(ValueError):
        f.to_dict()


def test_feature_to_tsv_bool_qualifier() -> None:
    f = Feature(type="source", id=1, qualifiers={"pseudo": [True]})
    tsv = f.to_tsv()
    assert tsv[0] == ["", "source", "", "pseudo", ""]


def test_feature_to_tsv_string_qualifier() -> None:
    f = Feature(type="gene", id=1, location="1..100", qualifiers={"note": ["test"]})
    tsv = f.to_tsv()
    assert tsv[0] == ["", "gene", "1..100", "note", "test"]


def test_feature_to_tsv_no_qualifier() -> None:
    f = Feature(type="UTR", id=1, location="1..50", qualifiers={})
    tsv = f.to_tsv()
    assert tsv[0] == ["", "UTR", "1..50", "", ""]


# === Entry tests ===


def test_entry_to_tsv_with_features() -> None:
    f1 = Feature(type="gene", id=1, location="1..100", qualifiers={"note": ["test"]})
    entry = Entry(id="entry1", name="entry1", features=[f1])
    tsv = entry.to_tsv()
    assert tsv[0] == ["entry1", "gene", "1..100", "note", "test"]


def test_entry_to_tsv_no_features() -> None:
    entry = Entry(id="entry2", name="entry2", features=[])
    tsv = entry.to_tsv()
    assert not tsv


# === MSS tests ===


def test_mss_parse_with_example_files(eg_complete: Dict[str, Path], eg_vrl: Dict[str, Path]) -> None:
    MSS.parse(eg_complete["ann"], eg_complete["seq"])
    MSS.parse(eg_vrl["ann"], eg_vrl["seq"])


def test_mss_to_tsv(mss_complete: MSS, mss_vrl: MSS) -> None:
    mss_complete.to_tsv()
    mss_vrl.to_tsv()


def test_mss_to_fasta(mss_complete: MSS, mss_vrl: MSS) -> None:
    mss_complete.to_fasta()
    mss_vrl.to_fasta()


def test_mss_from_json(eg_complete: Dict[str, Path], eg_vrl: Dict[str, Path]) -> None:
    MSS.from_json(eg_complete["v1_json"])
    MSS.from_json(eg_complete["v2_json"])
    MSS.from_json(eg_vrl["v1_json"])
    MSS.from_json(eg_vrl["v2_json"])
