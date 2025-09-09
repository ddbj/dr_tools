import tempfile
from pathlib import Path
from typing import Dict

from dr_tools.json2ann import json2ann


def test_json2ann_complete_v1(eg_complete: Dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        json2ann(
            mss_json_file=eg_complete["v1_json"],
            out_dir=tmpdir,
            out_prefix=None,
        )
        out_ann_file = Path(tmpdir).joinpath("mss.ann")
        out_seq_file = Path(tmpdir).joinpath("mss.fa")
        assert out_ann_file.exists()
        assert out_seq_file.exists()


def test_json2ann_complete_v2(eg_complete: Dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        json2ann(
            mss_json_file=eg_complete["v2_json"],
            out_dir=tmpdir,
            out_prefix=None,
        )
        out_ann_file = Path(tmpdir).joinpath("mss.ann")
        out_seq_file = Path(tmpdir).joinpath("mss.fa")
        assert out_ann_file.exists()
        assert out_seq_file.exists()


def test_json2ann_complete_with_prefix(eg_complete: Dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        json2ann(
            mss_json_file=eg_complete["v1_json"],
            out_dir=tmpdir,
            out_prefix="test_prefix",
        )
        out_ann_file = Path(tmpdir).joinpath("test_prefix.ann")
        out_seq_file = Path(tmpdir).joinpath("test_prefix.fa")
        assert out_ann_file.exists()
        assert out_seq_file.exists()


def test_json2ann_complete_no_outdir(eg_complete: Dict[str, Path]) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        json2ann(
            mss_json_file=eg_complete["v1_json"],
            out_dir=tmpdir,
            out_prefix="test_prefix",
        )
        out_ann_file = Path(tmpdir).joinpath("test_prefix.ann")
        out_seq_file = Path(tmpdir).joinpath("test_prefix.fa")
        assert out_ann_file.exists()
        assert out_seq_file.exists()
