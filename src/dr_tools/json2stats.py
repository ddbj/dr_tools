"""
JSON ファイルを読み込んで、ゲノムサイズ、遺伝子数などの統計情報を取得する
"""

import argparse
import json
import sys
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio.SeqFeature import Location
from ddbj_record.schema.v3 import Entry, Feature

from dr_tools.json_utils import features_by_entry, get_entry_id, load_json_to_ddbj_record_instance


def _sequence_length(entry: Entry) -> int:
    return len(entry.sequence or "")


def count_features(features: list[Feature]) -> dict[str, int]:
    """
    feature リストから、各種遺伝子の数をカウントする
    """
    feature_types = [feature.type or "" for feature in features]
    feature_counts = Counter(feature_types)
    return dict(feature_counts)


def get_genome_stats(entries: list[Entry]) -> dict[str, Any]:
    """
    エントリーリストからゲノムサイズ、N数のカウント等を行う
    GC 含量は N を除いた塩基に対する G と C の割合。N を除いた塩基が無ければ 0 にする
    """
    number_of_entries = len(entries)
    all_sequences = "".join([entry.sequence or "" for entry in entries])
    all_sequences = all_sequences.upper()
    genome_size = len(all_sequences)
    number_of_n = all_sequences.count("N")
    number_of_gc = all_sequences.count("G") + all_sequences.count("C")
    number_of_non_n = genome_size - number_of_n
    gc_content = number_of_gc / number_of_non_n if number_of_non_n else 0.0
    n_ratio = number_of_n / genome_size if genome_size else 0.0
    return {
        "genome_size": genome_size,
        "number_of_sequences": number_of_entries,
        "gc_content": gc_content,
        "gap_ratio": n_ratio,
    }


def get_N50(entries: list[Entry]) -> int:  # noqa: N802  公開している関数の名前を変えない
    """
    エントリーリストからN50を計算する。entry が無ければ 0 を返す
    """
    seq_lengths = [_sequence_length(entry) for entry in entries]
    seq_lengths.sort(reverse=True)
    total_length = sum(seq_lengths)
    half_length = total_length / 2
    cum_length = 0
    for length in seq_lengths:
        cum_length += length
        if cum_length >= half_length:
            return length

    return 0


def get_coding_ratio(entries: list[Entry], grouped_features: Mapping[str, list[Feature]]) -> float:
    """
    エントリーリストと entry ごとの feature からコーディング領域の割合を計算する
    CDSがオーバーラップしている可能性もあるが、考慮していないので簡易的な実装
    配列の長さの合計が 0 なら 0 を返す
    """

    def _get_feature_length(feature: Feature, seq_length: int, circular: bool) -> int:
        return len(Location.fromstring(feature.location or "", seq_length, circular))  # type: ignore[no-untyped-call]

    coding_length = 0
    total_length = 0
    for entry in entries:
        circular = entry.topology == "circular"
        seq_length = _sequence_length(entry)
        cds_features = [
            _get_feature_length(feature, seq_length, circular)
            for feature in grouped_features.get(get_entry_id(entry), [])
            if feature.type == "CDS"
        ]
        total_length += seq_length
        coding_length += sum(cds_features)
    return coding_length / total_length if total_length else 0.0


def json2stats(json_file: Path) -> dict[str, Any]:
    """
    JSON ファイルから統計情報を取得する
    JSON は record v1 / v2 / v3 のどれでもよい (v1 / v2 は v3 に変換してから読む)
    """
    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v3")

    entries = (record.sequences.entries if record.sequences is not None else None) or []
    grouped_features = features_by_entry(record)
    genome_stats = get_genome_stats(entries)
    feature_stats = count_features(record.features or [])
    n50 = get_N50(entries)
    coding_ratio = get_coding_ratio(entries, grouped_features)
    stats = {}
    stats.update(genome_stats)
    stats["N50"] = n50
    stats.update(feature_stats)
    stats["coding_ratio"] = coding_ratio

    return stats


def json2stats_for_dfast(
    json_file: Path, format_: bool = False, output_file: Path | None = None
) -> dict[str, Any] | None:
    """
    DFAST 用の統計情報を取得する
    format=True の場合、書式を整える
    Falseの場合、数値データとして返す

    # DFAST の統計情報の表示名　下記のラベルに変更
    {
        'genome_size': 'Total Length (bp)',
        'number_of_sequences': 'No. of Sequences',
        'gc_content': 'GC Content (%)',
        'gap_ratio': '',
        'N50': 'N50 (bp)',
        'CDS': 'No. of CDSs',
        'tRNA': 'No. of tRNA/tmRNA',
        'rRNA': 'No. of rRNA',
        'tmRNA': 'No. of tRNA/tmRNA',
        'repeat_region': 'No. of CRISPRS',
        'coding_ratio': 'Coding Ratio (%)'
    }
    """
    stats = json2stats(json_file)
    if not format_:
        dfast_stats = {
            "Total Length (bp)": stats["genome_size"],
            "No. of Sequences": stats["number_of_sequences"],
            "GC Content (%)": stats["gc_content"] * 100,
            "N50 (bp)": stats["N50"],
            "No. of CDSs": stats.get("CDS", 0),
            "No. of tRNA/tmRNA": stats.get("tRNA", 0) + stats.get("tmRNA", 0),
            "No. of rRNA": stats.get("rRNA", 0),
            "No. of CRISPRS": stats.get("repeat_region", 0),
            "Coding Ratio (%)": stats["coding_ratio"] * 100,
        }
    else:
        dfast_stats = {
            "Total Length (bp)": f"{stats['genome_size']:,d}",
            "No. of Sequences": f"{stats['number_of_sequences']:,d}",
            "GC Content (%)": f"{stats['gc_content'] * 100:0.1f}%",
            "N50 (bp)": f"{stats['N50']:,d}",
            "No. of CDSs": f"{stats.get('CDS', 0):,d}",
            "No. of tRNA/tmRNA": f"{stats.get('tRNA', 0) + stats.get('tmRNA', 0):,d}",
            "No. of rRNA": f"{stats.get('rRNA', 0):,d}",
            "No. of CRISPRS": f"{stats.get('repeat_region', 0):,d}",
            "Coding Ratio (%)": f"{stats['coding_ratio'] * 100:.1f}%",
        }
    if output_file:
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(dfast_stats, f, indent=4)
    else:
        return dfast_stats

    return None


def json2stats_for_dfast_main() -> None:
    """
    json2stats_for_dfast のコマンドラインインターフェース
    """
    # 引数のパース
    parser = argparse.ArgumentParser(description="Create statistics from JSON file for DFAST")
    parser.add_argument("json_file", type=str, help="JSON file")
    parser.add_argument("-f", "--format", action="store_true", help="When specified, format the output into strings")
    parser.add_argument("-o", "--output", type=str, help="Output file")

    # 引数のパース
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    stats = json2stats_for_dfast(args.json_file, format_=args.format, output_file=args.output)
    if stats:
        print(json.dumps(stats, indent=4))


if __name__ == "__main__":
    json2stats_for_dfast_main()
