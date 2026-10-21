from pathlib import Path
from urllib.parse import unquote

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import Location, SeqFeature
from Bio.SeqRecord import SeqRecord
from hypothesis import given
from hypothesis import strategies as st

from dr_tools.json2gff import (
    escape_attribute_value,
    feature_source,
    json2gff,
    json_to_gff,
    part_phases,
    seqrecords_to_gff,
)


def _feature(location: str, feature_type: str = "CDS", length: int = 100, **qualifiers: list[str]) -> SeqFeature:
    return SeqFeature(  # type: ignore[no-untyped-call]
        location=Location.fromstring(location, length=length, circular=True),  # type: ignore[no-untyped-call]
        type=feature_type,
        qualifiers=qualifiers,
    )


def _record(*features: SeqFeature, name: str = "seq1", length: int = 100) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=name, name=name, description="")
    record.features = list(features)
    return record


def _feature_lines(gff: str) -> list[list[str]]:
    body = gff.partition("##FASTA")[0]
    return [line.split("\t") for line in body.splitlines() if line and not line.startswith("#")]


# === The legacy layout ===


def test_json_to_gff_writes_the_legacy_dfast_columns(eg_complete: dict[str, Path]) -> None:
    lines = json_to_gff(eg_complete["v3_json"]).splitlines()

    assert lines[:4] == [
        "##gff-version 3",
        "##sequence-region chromosome 1 13",
        "##sequence-region pPLH-1 1 13",
        "##sequence-region pPLH-2 1 13",
    ]
    assert lines[4] == "\t".join(
        [
            "chromosome",
            "ab initio prediction:MetaGeneAnnotator",
            "CDS",
            "360",
            "1676",
            ".",
            "+",
            "0",
            (
                "ID=PLH_00010;product=chromosomal replication initiator protein DnaA;gene=dnaA;"
                "inference=COORDINATES:ab initio prediction:MetaGeneAnnotator,similar to AA sequence:RefSeq:WP_003641628.1;"
                "locus_tag=PLH_00010"
            ),
        ]
    )


def test_json_to_gff_trna_is_on_the_minus_strand_with_no_phase(eg_complete: dict[str, Path]) -> None:
    trna = next(c for c in _feature_lines(json_to_gff(eg_complete["v3_json"])) if c[2] == "tRNA")

    assert trna[1] == "Aragorn:1.2.38"
    assert trna[6:8] == ["-", "."]
    assert trna[8].startswith("ID=PLH_t00010;")


def test_json_to_gff_skips_gene_source_and_assembly_gap() -> None:
    record = _record(
        _feature("1..30", "source"),
        _feature("1..30", "gene", locus_tag=["L_1"]),
        _feature("40..50", "assembly_gap"),
        _feature("1..30", "CDS", locus_tag=["L_1"]),
    )

    assert [c[2] for c in _feature_lines(seqrecords_to_gff([record]))] == ["CDS"]


def test_json_to_gff_fasta_section_has_each_entry_unwrapped(full_record: dict[str, Path]) -> None:
    fasta = json_to_gff(full_record["v3"]).split("##FASTA\n")[1].splitlines()

    assert [line for line in fasta if line.startswith(">")] == [">chromosome", ">pPLH-1", ">pPLH-2"]
    assert len(fasta) == 6
    assert set(fasta[1]) <= set("ACGTN")


def test_json_to_gff_long_product_has_no_inserted_space(full_record: dict[str, Path]) -> None:
    """GenBank text wraps a long product; the legacy converter had to remove the space it left."""
    gff = json_to_gff(full_record["v3"])

    for columns in _feature_lines(gff):
        for attribute in columns[8].split(";"):
            if attribute.startswith("product="):
                assert "  " not in unquote(attribute)


def test_json_to_gff_same_output_from_v1_v2_v3(full_record: dict[str, Path]) -> None:
    assert json_to_gff(full_record["v1"]) == json_to_gff(full_record["v2"]) == json_to_gff(full_record["v3"])


def test_json2gff_writes_the_file(eg_complete: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path.joinpath("genome.gff3")
    json2gff(eg_complete["v3_json"], out)

    assert out.read_text(encoding="utf-8") == json_to_gff(eg_complete["v3_json"])


# === Fixes over the legacy converter ===


def test_location_of_several_parts_is_written_part_by_part_with_one_id() -> None:
    """The legacy converter wrote one line from the smallest to the largest position, which for a CDS over the
    origin of a circular genome covered almost the whole genome."""
    record = _record(_feature("join(91..100,1..20)", locus_tag=["L_1"]))

    lines = _feature_lines(seqrecords_to_gff([record]))

    assert [(c[3], c[4], c[6]) for c in lines] == [("91", "100", "+"), ("1", "20", "+")]
    assert all(c[8].startswith("ID=L_1;") for c in lines)


def test_phase_of_each_part_completes_the_codon_left_open_by_the_part_before() -> None:
    # 10 bases from phase 0 leave 1 base of a codon open, so the next part starts 2 bases into it.
    assert part_phases(_feature("join(91..100,1..20)")) == ["0", "2"]
    # codon_start 2: 9 bases after the first one make 3 codons, so the next part starts a codon.
    assert part_phases(_feature("join(1..10,21..40)", codon_start=["2"])) == ["1", "0"]


def test_phase_on_the_minus_strand_follows_the_parts_in_transcription_order() -> None:
    feature = _feature("complement(join(1..10,21..31))")

    lines = _feature_lines(seqrecords_to_gff([_record(feature)]))

    assert [(c[3], c[4], c[6], c[7]) for c in lines] == [("21", "31", "-", "0"), ("1", "10", "-", "1")]


def test_feature_other_than_cds_has_no_phase() -> None:
    assert part_phases(_feature("join(1..10,21..31)", "rRNA")) == [".", "."]


def test_phase_of_a_cds_with_codon_start_3() -> None:
    assert part_phases(_feature("1..30", codon_start=["3"])) == ["2"]


@pytest.mark.parametrize(
    ("inferences", "expected"),
    [
        (["COORDINATES:profile:Barrnap:0.8"], "Barrnap:0.8"),
        (["ab initio prediction"], "dfast"),
        (["ab initio prediction", "COORDINATES:alignment:CRT:1.2"], "CRT:1.2"),
        ([], "dfast"),
    ],
)
def test_source_column(inferences: list[str], expected: str) -> None:
    """An inference without two ":" made the legacy converter fail."""
    assert feature_source(_feature("1..30", inference=inferences)) == expected


def test_valueless_qualifier_is_written_as_true() -> None:
    """json_to_seqrecords turns /pseudo into an empty string."""
    record = _record(_feature("1..30", locus_tag=["L_1"], pseudo=[""]))

    assert "pseudo=true" in _feature_lines(seqrecords_to_gff([record]))[0][8].split(";")


def test_codon_start_and_transl_table_are_not_attributes() -> None:
    record = _record(_feature("1..30", locus_tag=["L_1"], codon_start=["1"], transl_table=["11"]))

    assert _feature_lines(seqrecords_to_gff([record]))[0][8] == "ID=L_1;locus_tag=L_1"


def test_feature_without_locus_tag_has_no_id() -> None:
    record = _record(_feature("1..30", "repeat_region", note=["CRISPR"]))

    assert _feature_lines(seqrecords_to_gff([record]))[0][8] == "note=CRISPR"


@given(st.text())
def test_escaped_attribute_value_decodes_to_the_value(value: str) -> None:
    escaped = escape_attribute_value(value)

    assert not set(escaped) & set(";=&,\t\n\r")
    assert unquote(escaped) == value


def test_attribute_value_with_reserved_characters_stays_in_its_attribute() -> None:
    record = _record(_feature("1..30", locus_tag=["L_1"], product=["a;b=c,d&e%f\ng"]))

    attributes = _feature_lines(seqrecords_to_gff([record]))[0][8].split(";")

    assert "product=a%3Bb%3Dc%2Cd%26e%25f%0Ag" in attributes


def test_seqid_with_characters_gff_does_not_allow_is_escaped_everywhere() -> None:
    gff = seqrecords_to_gff([_record(_feature("1..30"), name="seq 1")])

    assert "##sequence-region seq%201 1 100" in gff
    assert _feature_lines(gff)[0][0] == "seq%201"
    assert ">seq%201" in gff
