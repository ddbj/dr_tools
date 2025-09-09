import tempfile
from pathlib import Path
from typing import Dict

from dr_tools.json2biopython import json_to_seqrecords


def test_json_to_seq_records_complete_v1(eg_complete: Dict[str, Path]) -> None:
    json_to_seqrecords(eg_complete["v1_json"])


def test_json_to_seq_records_complete_v2(eg_complete: Dict[str, Path]) -> None:
    json_to_seqrecords(eg_complete["v2_json"])
