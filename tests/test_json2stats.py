import json
from pathlib import Path
from typing import Any

import pytest
from ddbj_record.schema.v3 import Entry
from hypothesis import given
from hypothesis import strategies as st

from dr_tools.json2stats import get_genome_stats, get_N50, json2stats, json2stats_for_dfast


def _write_v3_record(tmp_path: Path, sequences: list[str]) -> Path:
    record = {
        "schema_version": "v3",
        "sequences": {"entries": [{"alias": f"seq{i}", "sequence": seq} for i, seq in enumerate(sequences)]},
    }
    json_file = tmp_path.joinpath("record.json")
    json_file.write_text(json.dumps(record), encoding="utf-8")
    return json_file


@pytest.mark.parametrize("version", ["v1_json", "v2_json", "v3_json"])
def test_json2stats_for_dfast_complete_each_record_version_returns_labels(
    version: str, eg_complete: dict[str, Path]
) -> None:
    stats = json2stats_for_dfast(json_file=eg_complete[version])
    assert stats is not None
    assert set(stats) == {
        "Total Length (bp)",
        "No. of Sequences",
        "GC Content (%)",
        "N50 (bp)",
        "No. of CDSs",
        "No. of tRNA/tmRNA",
        "No. of rRNA",
        "No. of CRISPRS",
        "Coding Ratio (%)",
    }


def test_json2stats_full_record_v3_counts_entries_and_features(full_record: dict[str, Path]) -> None:
    stats = json2stats(full_record["v3"])
    assert stats["number_of_sequences"] == 3
    assert stats["genome_size"] == 2277985 + 81630 + 40971
    assert stats["N50"] == 2277985
    assert stats["CDS"] > 2000
    assert "source" not in stats
    assert 0 < stats["coding_ratio"] < 1


def test_json2stats_same_result_from_v1_v2_v3_of_same_ann(full_record: dict[str, Path]) -> None:
    results = [json2stats(full_record[version]) for version in ("v1", "v2", "v3")]
    assert results[0] == results[1] == results[2]


def test_json2stats_for_dfast_format_true_writes_formatted_strings_to_output_file(
    full_record: dict[str, Path], tmp_path: Path
) -> None:
    output_file = tmp_path.joinpath("stats.json")
    returned = json2stats_for_dfast(json_file=full_record["v3"], format_=True, output_file=output_file)
    assert returned is None
    written: dict[str, Any] = json.loads(output_file.read_text(encoding="utf-8"))
    assert written["No. of Sequences"] == "3"
    assert written["Total Length (bp)"] == "2,400,586"
    assert written["GC Content (%)"].endswith("%")


# --- GC 含量・N の割合・N50・coding ratio ---


def test_get_genome_stats_gc_content_counts_g_and_c_over_length_without_n() -> None:
    stats = get_genome_stats([Entry(alias="s1", sequence="GGNN")])
    assert stats["gc_content"] == 1.0
    assert stats["gap_ratio"] == 0.5


def test_get_genome_stats_lowercase_sequence_is_counted_like_uppercase() -> None:
    stats = get_genome_stats([Entry(alias="s1", sequence="gcatnn")])
    assert stats["gc_content"] == 0.5


@given(st.lists(st.text(alphabet="ACGTN", max_size=30), min_size=1, max_size=4))
def test_get_genome_stats_gc_content_is_gc_over_non_n_bases(sequences: list[str]) -> None:
    stats = get_genome_stats([Entry(alias=f"s{i}", sequence=seq) for i, seq in enumerate(sequences)])
    joined = "".join(sequences)
    non_n = len(joined) - joined.count("N")
    expected = (joined.count("G") + joined.count("C")) / non_n if non_n else 0.0
    assert stats["gc_content"] == pytest.approx(expected)
    assert 0.0 <= stats["gc_content"] <= 1.0
    assert 0.0 <= stats["gap_ratio"] <= 1.0


def test_get_genome_stats_all_n_returns_zero_gc_content_and_full_gap_ratio() -> None:
    stats = get_genome_stats([Entry(alias="s1", sequence="NNNN")])
    assert stats["gc_content"] == 0.0
    assert stats["gap_ratio"] == 1.0


def test_get_genome_stats_no_entries_returns_zero_ratios() -> None:
    stats = get_genome_stats([])
    assert stats == {"genome_size": 0, "number_of_sequences": 0, "gc_content": 0.0, "gap_ratio": 0.0}


def test_get_n50_no_entries_returns_zero() -> None:
    assert get_N50([]) == 0


def test_json2stats_record_without_entries_returns_zero_values(tmp_path: Path) -> None:
    stats = json2stats(_write_v3_record(tmp_path, []))
    assert stats["genome_size"] == 0
    assert stats["N50"] == 0
    assert stats["gc_content"] == 0.0
    assert stats["coding_ratio"] == 0.0


def test_json2stats_for_dfast_format_true_record_without_entries_formats_zeros(tmp_path: Path) -> None:
    stats = json2stats_for_dfast(_write_v3_record(tmp_path, []), format_=True)
    assert stats is not None
    assert stats["N50 (bp)"] == "0"
    assert stats["GC Content (%)"] == "0.0%"
    assert stats["Coding Ratio (%)"] == "0.0%"


def test_json2stats_for_dfast_sequence_with_gap_reports_gc_of_non_gap_bases(tmp_path: Path) -> None:
    stats = json2stats_for_dfast(_write_v3_record(tmp_path, ["GCATNNNN"]))
    assert stats is not None
    assert stats["GC Content (%)"] == 50.0
