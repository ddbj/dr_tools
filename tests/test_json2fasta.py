from pathlib import Path

import pytest

from dr_tools.json2fasta import json2fasta

OUTPUT_FILES = ["genome.fna", "cds.fna", "misc_rnas.fna", "protein.faa"]


@pytest.mark.parametrize("version", ["v1_json", "v2_json", "v3_json"])
def test_json2fasta_complete_writes_all_output_files(
    version: str, eg_complete: dict[str, Path], tmp_path: Path
) -> None:
    json2fasta(json_file=eg_complete[version], out_dir=tmp_path)
    for name in OUTPUT_FILES:
        assert tmp_path.joinpath(name).exists()


def test_json2fasta_v3_full_record_cds_and_protein_have_locus_tag_descriptions(
    full_record: dict[str, Path], tmp_path: Path
) -> None:
    json2fasta(json_file=full_record["v3"], out_dir=tmp_path)
    cds = tmp_path.joinpath("cds.fna").read_text(encoding="utf-8").splitlines()
    protein = tmp_path.joinpath("protein.faa").read_text(encoding="utf-8").splitlines()
    assert (
        cds[0]
        == ">PLH_00010 chromosomal replication initiator protein DnaA [gene=dnaA] [type=CDS] [location=chromosome:360..1676]"
    )
    assert len(cds[1]) == 1676 - 360 + 1
    assert set(cds[1]) <= set("ACGT")
    assert protein[0] == cds[0]
    assert protein[1].startswith("M")
    genome = tmp_path.joinpath("genome.fna").read_text(encoding="utf-8")
    assert [line for line in genome.splitlines() if line.startswith(">")] == [">chromosome", ">pPLH-1", ">pPLH-2"]


def test_json2fasta_same_output_from_v1_v2_v3_of_same_ann(full_record: dict[str, Path], tmp_path: Path) -> None:
    outputs = []
    for version in ("v1", "v2", "v3"):
        out_dir = tmp_path.joinpath(version)
        json2fasta(json_file=full_record[version], out_dir=out_dir)
        outputs.append({name: out_dir.joinpath(name).read_text(encoding="utf-8") for name in OUTPUT_FILES})
    assert outputs[0] == outputs[1] == outputs[2]


def test_json2fasta_rna_features_go_to_misc_rnas(full_record: dict[str, Path], tmp_path: Path) -> None:
    json2fasta(json_file=full_record["v3"], out_dir=tmp_path)
    rna = tmp_path.joinpath("misc_rnas.fna").read_text(encoding="utf-8")
    assert "[type=tRNA]" in rna
    assert "[type=CDS]" not in rna
