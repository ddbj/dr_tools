from pathlib import Path
from typing import Dict

from dr_tools.json2stats import json2stats_for_dfast


def test_json2fasta_complete_v1(eg_complete: Dict[str, Path]) -> None:
    json2stats_for_dfast(json_file=eg_complete["v1_json"])


def test_json2fasta_complete_v2(eg_complete: Dict[str, Path]) -> None:
    json2stats_for_dfast(json_file=eg_complete["v2_json"])
