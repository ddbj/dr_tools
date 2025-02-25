import json
from pathlib import Path
from typing import Dict, Tuple, List
from collections import Counter

from Bio.SeqFeature import Location

"""
JSON ファイルを読み込んで、ゲノムサイズ、遺伝子数などの統計情報を取得する
"""

def count_features(entries: List) -> Dict:
    """
    エントリーリストから、各種遺伝子の数をカウントする
    """

    featute_types = [feature["type"] for entry in entries for feature in entry["features"] if feature["type"] != "source"]
    feature_counts = Counter(featute_types)
    return dict(feature_counts)

def get_genome_stats(entries: List) -> Dict:
    """
    エントリーリストからゲノムサイズ、N数のカウント等を行う
    """
    number_of_entries = len(entries)
    all_sequences = "".join([entry["sequence"] for entry in entries])
    all_sequences = all_sequences.upper()
    genome_size = len(all_sequences)
    number_of_N = all_sequences.count("N")
    number_of_gc = all_sequences.count("G") + all_sequences.count("C")
    gc_content = (number_of_gc - number_of_N) / (genome_size - number_of_N)
    n_ratio = number_of_N / genome_size
    D = {
        "genome_size": genome_size,
        "number_of_sequences": number_of_entries,
        "gc_content": gc_content,
        "gap_ratio": n_ratio
    }
    return D

def get_N50(entries: List) -> int:
    """
    エントリーリストからN50を計算する
    """
    def _get_length(entry: Dict) -> int:
        return entry.get("length") or len(entry["sequence"])
    
    seq_lengths = [_get_length(entry) for entry in entries]
    seq_lengths.sort(reverse=True)
    total_length = sum(seq_lengths)
    half_length = total_length / 2
    cum_length = 0
    for length in seq_lengths:
        cum_length += length
        if cum_length >= half_length:
            return length

def get_coding_ratio(entries: List) -> float:
    """
    エントリーリストからコーディング領域の割合を計算する
    CDSがオーバーラップしている可能性もあるが、考慮していないので簡易的な実装
    """
    def _get_feature_length(feature, seq_length, circular):
        return len(Location.fromstring(feature["location"], seq_length, circular))

    coding_length = 0
    total_length = 0
    for entry in entries:
        circular = entry.get("topology", "linear") == "circular"
        seq_length = entry.get("length") or len(entry["sequence"])
        cds_features = [_get_feature_length(feature, seq_length, circular) for feature in entry["features"] if feature["type"] == "CDS"]
        total_length += seq_length
        coding_length += sum(cds_features)
    coding_ratio = coding_length / total_length
    return coding_ratio

def json2stats(json_file: Path) -> Dict:
    """
    JSON ファイルから統計情報を取得する
    """
    json_dat = json.load(open(json_file))
    entries = json_dat.get("ENTRIES", [])
    genome_stats = get_genome_stats(entries)
    feature_stats = count_features(entries)
    n50 = get_N50(entries)
    coding_ratio = get_coding_ratio(entries)
    stats = {}
    stats.update(genome_stats)
    stats["N50"] = n50
    stats.update(feature_stats)
    stats["coding_ratio"] = coding_ratio

    return stats

def json2stats_for_dfast(json_file: Path, format=False, output_file: Path | None =None) -> Dict | None:
    """
    DFAST 用の統計情報を取得する
    format=True の場合、書式を整える
    Falseの場合、数値データとして返す
    """
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
        'coding_ratio':'Coding Ratio (%)'
    }
    stats = json2stats(json_file)
    if not format:
        D = {
            'Total Length (bp)': stats['genome_size'],
            'No. of Sequences': stats['number_of_sequences'],
            'GC Content (%)': stats['gc_content'] * 100,
            'N50 (bp)': stats['N50'],
            'No. of CDSs': stats.get('CDS', 0),
            'No. of tRNA/tmRNA': stats.get('tRNA', 0) + stats.get('tmRNA', 0),
            'No. of rRNA': stats.get('rRNA', 0),
            'No. of CRISPRS': stats.get('repeat_region', 0),
            'Coding Ratio (%)': stats['coding_ratio'] * 100
        }
    else:
        D = {
            'Total Length (bp)': f"{stats['genome_size']:,d}",
            'No. of Sequences': f"{stats['number_of_sequences']:,d}",
            'GC Content (%)': f"{stats['gc_content']*100:0.1f}%",
            'N50 (bp)': f"{stats['N50']:,d}",
            'No. of CDSs': f"{stats.get('CDS', 0):,d}",
            'No. of tRNA/tmRNA': f"{stats.get('tRNA', 0) + stats.get('tmRNA', 0):,d}",
            'No. of rRNA': f"{stats.get('rRNA', 0):,d}",
            'No. of CRISPRS': f"{stats.get('repeat_region', 0):,d}",
            'Coding Ratio (%)': f"{stats['coding_ratio']*100:.1f}%"
        }
    if output_file:
        with open(output_file, "w") as f:
            json.dump(D, f, indent=4)
    else:
        return D

def json2stats_for_dfast_main():
    """
    json2stats_for_dfast のコマンドラインインターフェース
    """
    import argparse
    import sys

    # 引数のパース
    parser = argparse.ArgumentParser(description='Create statistics from JSON file for DFAST')
    parser.add_argument('json_file', type=str, help='JSON file')
    parser.add_argument('-f', '--format', action='store_true', help='When specified, format the output into strings')
    parser.add_argument('-o', '--output', type=str, help='Output file')

    # 引数のパース
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    stats = json2stats_for_dfast(args.json_file, format=args.format, output_file=args.output)
    if stats:
        print(json.dumps(stats, indent=4))

if __name__ == "__main__":
    json2stats_for_dfast_main()
