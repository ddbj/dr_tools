import json
from pathlib import Path
from typing import Any, Dict

import pytest

from dr_tools.json2stats import json2stats, json2stats_for_dfast


@pytest.mark.parametrize("version", ["v1_json", "v2_json", "v3_json"])
def test_json2stats_for_dfast_complete_each_record_version_returns_labels(version: str, eg_complete: Dict[str, Path]) -> None:
    stats = json2stats_for_dfast(json_file=eg_complete[version])
    assert stats is not None
    assert set(stats) == {
        "Total Length (bp)", "No. of Sequences", "GC Content (%)", "N50 (bp)", "No. of CDSs",
        "No. of tRNA/tmRNA", "No. of rRNA", "No. of CRISPRS", "Coding Ratio (%)",
    }


def test_json2stats_full_record_v3_counts_entries_and_features(full_record: Dict[str, Path]) -> None:
    stats = json2stats(full_record["v3"])
    assert stats["number_of_sequences"] == 3
    assert stats["genome_size"] == 2277985 + 81630 + 40971
    assert stats["N50"] == 2277985
    assert stats["CDS"] > 2000
    assert "source" not in stats
    assert 0 < stats["coding_ratio"] < 1


def test_json2stats_same_result_from_v1_v2_v3_of_same_ann(full_record: Dict[str, Path]) -> None:
    results = [json2stats(full_record[version]) for version in ("v1", "v2", "v3")]
    assert results[0] == results[1] == results[2]


def test_json2stats_for_dfast_format_true_writes_formatted_strings_to_output_file(full_record: Dict[str, Path], tmp_path: Path) -> None:
    output_file = tmp_path.joinpath("stats.json")
    returned = json2stats_for_dfast(json_file=full_record["v3"], format_=True, output_file=output_file)
    assert returned is None
    written: Dict[str, Any] = json.loads(output_file.read_text(encoding="utf-8"))
    assert written["No. of Sequences"] == "3"
    assert written["Total Length (bp)"] == "2,400,586"
    assert written["GC Content (%)"].endswith("%")
