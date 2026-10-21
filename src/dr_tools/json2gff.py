"""
record JSON から GFF3 を作る

DFAST (旧 DFAST の genbank2gff.py と dfast_vrl の dfv/genbank2gff.py) が出してきた形の GFF3 を、
json_to_seqrecords が作る SeqRecord から書く。GenBank のテキストを経由しないので、旧 DFAST が
fix_blank_in_products で直していた、折り返しで product に入る空白は生じない

- entry ごとに、##sequence-region の行と、その entry の feature の行を続けて書く
- 1 列目は entry の名前 (SeqRecord.name)。##sequence-region と ##FASTA の見出しも同じ名前を使う。
  GFF3 の seqid に使えない文字は %XX にする
- 2 列目は /inference の最後の 2 つの要素 ("ab initio prediction:MetaGeneAnnotator" など)。無ければ "dfast"
- gene、source、assembly_gap は書かない
- 9 列目は ID (locus_tag があるとき) と、codon_start と transl_table を除く全ての qualifier
- 複数の区間からなる位置 (join、原点をまたぐ位置) は、区間ごとに同じ ID の行を書く。CDS の phase は区間ごとに求める
"""

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from dr_tools.json2biopython import json_to_seqrecords

GFF_VERSION_LINE = "##gff-version 3"
# 書かない feature (gene は CDS などと同じ locus_tag を持ち、source と assembly_gap は配列の性質を表すため)
SKIPPED_FEATURE_TYPES = frozenset({"gene", "source", "assembly_gap"})
# GFF の属性にしない qualifier (phase と配列の翻訳に使うもの)
SKIPPED_QUALIFIERS = frozenset({"codon_start", "transl_table"})
DEFAULT_SOURCE = "dfast"
# 2 列目に使う /inference の要素の数 (末尾から)
SOURCE_ELEMENTS = 2
# 値を持たない qualifier (/pseudo など) の値。json_to_seqrecords は空文字列にする
VALUELESS_QUALIFIER_VALUE = "true"

# 列にも属性の値にも書けない文字 (タブと改行を含む制御文字)。%XX にする
_CONTROL_CHARACTERS = frozenset([chr(i) for i in range(0x20)] + ["\x7f"])
# GFF3 の seqid にそのまま書ける文字。それ以外は %XX にする
_SEQID_SAFE = frozenset("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|")


def _escape(value: str, reserved: str) -> str:
    """% と制御文字と reserved の文字を %XX にする"""
    return "".join(f"%{ord(c):02X}" if c == "%" or c in _CONTROL_CHARACTERS or c in reserved else c for c in value)


def escape_attribute_value(value: str) -> str:
    """属性の値に使えない文字 (; = & , と、タブ、改行などの制御文字、%) を %XX にする"""
    return _escape(value, ";=&,")


def escape_column(value: str) -> str:
    """9 列目以外の列に使えない文字 (タブ、改行などの制御文字、%) を %XX にする"""
    return _escape(value, "")


def escape_seqid(value: str) -> str:
    return "".join(c if c in _SEQID_SAFE else f"%{ord(c):02X}" for c in value)


def feature_source(feature: SeqFeature) -> str:
    """2 列目。/inference の最後の 2 つの要素を ":" でつなぐ。要素が 2 つに満たなければ次の /inference を見る"""
    for inference in feature.qualifiers.get("inference", []):
        parts = inference.split(":")
        if len(parts) >= SOURCE_ELEMENTS:
            return ":".join(parts[-SOURCE_ELEMENTS:])
    return DEFAULT_SOURCE


def first_phase(feature: SeqFeature) -> int:
    """CDS の最初の区間の phase (codon_start - 1)。読めなければ 0"""
    try:
        codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    except ValueError:
        return 0
    return codon_start - 1 if codon_start in (1, 2, 3) else 0


def part_phases(feature: SeqFeature) -> list[str]:
    """区間ごとの phase。CDS 以外は "."

    区間は転写の向きに並ぶ (biopython の CompoundLocation.parts)。次の区間の phase は、それまでの区間の
    最後の codon を完成させるのに要る塩基の数である
    """
    parts = feature.location.parts if feature.location is not None else []
    if feature.type != "CDS":
        return ["."] * len(parts)
    phases = []
    phase = first_phase(feature)
    for part in parts:
        phases.append(str(phase))
        leftover = (len(part) - phase) % 3
        phase = (3 - leftover) % 3
    return phases


def feature_attributes(feature: SeqFeature) -> str:
    """9 列目。locus_tag があれば ID にし、続けて qualifier を名前=値,値 の形で並べる"""
    attributes = []
    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
    if locus_tag:
        attributes.append(f"ID={escape_attribute_value(locus_tag)}")
    for key, values in feature.qualifiers.items():
        if key in SKIPPED_QUALIFIERS:
            continue
        written = [escape_attribute_value(value or VALUELESS_QUALIFIER_VALUE) for value in values]
        attributes.append(f"{escape_attribute_value(key)}={','.join(written)}")
    return ";".join(attributes)


def feature_lines(seqid: str, feature: SeqFeature) -> list[str]:
    """1 つの feature の行。区間ごとに 1 行で、全ての行が同じ ID と属性を持つ。位置の無い feature は書かない"""
    if feature.location is None:
        return []
    source = escape_column(feature_source(feature))
    attributes = feature_attributes(feature)
    lines = []
    for part, phase in zip(feature.location.parts, part_phases(feature), strict=True):
        if part.strand == -1:
            strand = "-"
        elif part.strand == 1:
            strand = "+"
        else:
            strand = "."
        columns = [seqid, source, escape_column(feature.type), str(int(part.start) + 1), str(int(part.end))]
        lines.append("\t".join([*columns, ".", strand, phase, attributes]))
    return lines


def seqrecords_to_gff(seq_records: Iterable[SeqRecord]) -> str:
    """
    SeqRecord の一覧を GFF3 の文字列にする
    entry ごとに、##sequence-region の行を書き、続けてその entry の feature を書く。配列は末尾の ##FASTA に、
    大文字で、折り返さずに書く
    """
    records = list(seq_records)
    lines = [GFF_VERSION_LINE]
    for record in records:
        seqid = escape_seqid(record.name)
        lines.append(f"##sequence-region {seqid} 1 {len(record.seq or '')}")
        for feature in record.features:
            if feature.type not in SKIPPED_FEATURE_TYPES:
                lines.extend(feature_lines(seqid, feature))
    lines.append("##FASTA")
    for record in records:
        lines.append(f">{escape_seqid(record.name)}")
        # json2fasta の genome.fna と同じく大文字にする
        lines.append(str(record.seq or "").upper())
    return "\n".join(lines) + "\n"


def json_to_gff(json_file: Path) -> str:
    """
    record JSON を GFF3 の文字列にする
    JSON は record v1 / v2 / v3 のどれでもよい (v1 / v2 は v3 に変換してから読む)
    """
    return seqrecords_to_gff(json_to_seqrecords(json_file))


def json2gff(json_file: Path, out_file: Path) -> None:
    """record JSON を GFF3 のファイルにする"""
    Path(out_file).write_text(json_to_gff(json_file), encoding="utf-8")


def json2gff_main() -> None:
    """
    JSONファイルをGFF3ファイルに変換する
    出力 (-o, --output_file) が指定されていない場合は、標準出力に出力する
    """
    parser = argparse.ArgumentParser(description="Convert JSON file to GFF3 file")
    parser.add_argument("input_file", type=str, help="Input JSON file")
    parser.add_argument("-o", "--output_file", type=str, help="Output GFF3 file", default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    if args.output_file:
        json2gff(Path(args.input_file), Path(args.output_file))
    else:
        sys.stdout.write(json_to_gff(Path(args.input_file)))
