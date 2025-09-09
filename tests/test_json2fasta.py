import tempfile
from pathlib import Path
from typing import Dict

from dr_tools.json2fasta import json2fasta


def test_json2fasta_complete_v1(eg_complete: Dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        json2fasta(
            json_file=eg_complete["v1_json"],
            out_dir=out_dir,
        )
        genome_fasta_file = out_dir / Path("genome.fna")
        cds_fasta_file = out_dir / Path("cds.fna")
        rna_fasta_file = out_dir / Path("misc_rnas.fna")
        protein_fasta_file = out_dir / Path("protein.faa")

        assert genome_fasta_file.exists()
        assert cds_fasta_file.exists()
        assert rna_fasta_file.exists()
        assert protein_fasta_file.exists()


def test_json2fasta_complete_v2(eg_complete: Dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        json2fasta(
            json_file=eg_complete["v2_json"],
            out_dir=out_dir,
        )
        genome_fasta_file = out_dir / Path("genome.fna")
        cds_fasta_file = out_dir / Path("cds.fna")
        rna_fasta_file = out_dir / Path("misc_rnas.fna")
        protein_fasta_file = out_dir / Path("protein.faa")

        assert genome_fasta_file.exists()
        assert cds_fasta_file.exists()
        assert rna_fasta_file.exists()
        assert protein_fasta_file.exists()
