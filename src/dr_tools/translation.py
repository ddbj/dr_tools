"""
CDS の翻訳

dr_tools の中で CDS をアミノ酸配列にするのはこのモジュールだけ。
Bio.Seq.translate を翻訳の本体に使い、transl_except は翻訳後の後処理として当てる。

実装は https://github.com/nigyta/translate_with_exception の
translate_cds_with_transl_except を取り込んだもの。NCBI efetch がリージョン指定で返す
GenBank ファイル向けの座標補正は落としてある (dr_tools は record JSON から
SeqFeature を組み立てるので、その経路が存在しない)。

対応する transl_except の書式:
  (pos:START..END,aa:TERM)            - 末端の部分ストップコドン (+ 鎖)
  (pos:POS,aa:TERM)                   - 1 塩基だけの末端ストップコドン (+ 鎖)
  (pos:complement(START..END),aa:Sec) - - 鎖の内部特殊コドン
  (pos:START..END,aa:OTHER)           - 非標準アミノ酸 -> X
  (pos:START..END,aa:Met)             - AUG 以外の開始コドン -> M
"""
import logging
import re
import warnings
from typing import cast

from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

logger = logging.getLogger(__name__)

# 3 文字表記 (と特殊キーワード) から 1 文字表記へのアミノ酸の対応
_AA_3TO1: dict[str, str] = {**IUPACData.protein_letters_3to1_extended}
_AA_3TO1.update({
    "Term": "*",   # ストップコドン (部分コドンの場合を含む)
    "TERM": "*",
    "Other": "X",  # 非標準・未指定のアミノ酸
    "OTHER": "X",
})

_TRANSL_EXCEPT_PATTERN = re.compile(r"\(pos:(complement\()?(\d+)(?:\.\.(\d+))?(?(1)\)),aa:(\w+)\)")


def _get_codon_table(table_id: int | str) -> CodonTable.CodonTable:
    """NCBI の数値 ID または名前から CodonTable を返す"""
    try:
        table = CodonTable.ambiguous_generic_by_id[int(table_id)]
    except (ValueError, TypeError, KeyError):
        table = CodonTable.ambiguous_generic_by_name[str(table_id)]

    return cast(CodonTable.CodonTable, table)


def _parent_pos_to_extracted_pos(parent_pos_0based: int, feature: SeqFeature) -> int:
    """
    親配列上の 0-based の位置を、feature.extract() が返す配列上の 0-based の位置に変換する

    SimpleLocation と CompoundLocation (join / order) のどちらにも対応する。
    CompoundLocation では extract() が連結する順 (location.parts の順) で辿る。
    """
    loc = feature.location

    if isinstance(loc, SimpleLocation):
        if loc.strand == -1:
            return int(loc.end) - 1 - parent_pos_0based
        return parent_pos_0based - int(loc.start)

    if isinstance(loc, CompoundLocation):
        cumulative = 0
        for part in loc.parts:
            p_start, p_end = int(part.start), int(part.end)
            if p_start <= parent_pos_0based < p_end:
                if part.strand == -1:
                    return cumulative + (p_end - 1 - parent_pos_0based)
                return cumulative + (parent_pos_0based - p_start)
            cumulative += p_end - p_start
        raise ValueError(f"Parent position {parent_pos_0based} is not within any part of feature location {loc}")

    raise TypeError(f"Unsupported location type: {type(loc)}")


def _parse_transl_except(qualifier_value: str, feature: SeqFeature, start_offset: int) -> tuple[int, str]:
    """
    transl_except の値を 1 つ解釈し、(コドンの位置, 1 文字表記のアミノ酸) を返す
    コドンの位置は start_offset を当てたあとのアミノ酸配列上の 0-based の位置

    Args:
        qualifier_value: qualifier の生の値。例: '(pos:4263..4264,aa:TERM)'
        feature:         その qualifier を持つ SeqFeature (位置の変換に使う)
        start_offset:    codon_start - 1 (0-based のフレームのずれ)
    """
    match = _TRANSL_EXCEPT_PATTERN.fullmatch(qualifier_value.strip())
    if not match:
        raise ValueError(f"Cannot parse transl_except qualifier: {qualifier_value!r}")

    is_complement = bool(match.group(1))
    pos1 = int(match.group(2))  # 1-based, 親配列上の座標
    pos2 = int(match.group(3)) if match.group(3) else pos1
    aa_name = match.group(4)

    aa = _AA_3TO1.get(aa_name) or _AA_3TO1.get(aa_name.capitalize())
    if aa is None:
        raise ValueError(f"Unknown amino acid in transl_except: {aa_name!r}")

    # 切り出した配列の中で「コドンの先頭の塩基」に当たるのはどちらの端かを決める
    #   + 鎖 / complement 無し -> 左端 (pos1)
    #   - 鎖 / complement 有り -> 右端 (pos2)。切り出した配列は逆相補になっているため
    parent_pos_0based = (pos2 if is_complement else pos1) - 1

    pos_in_feature = _parent_pos_to_extracted_pos(parent_pos_0based, feature)

    return (pos_in_feature - start_offset) // 3, aa


def translate_cds_with_transl_except(feature: SeqFeature, parent_seq: Seq, stop_symbol: str = "*") -> Seq:
    """
    CDS の SeqFeature を翻訳する。transl_except があればその位置のアミノ酸を差し替える

    SeqFeature.translate(cds=True) と比べて次を扱える。
      - 末端の部分ストップコドン (ゲノム上 1-2 塩基, aa:TERM)
      - セレノシステイン等の内部特殊コドン (aa:Sec -> 'U')
      - 非標準アミノ酸 (aa:OTHER -> 'X')
      - AUG 以外の開始コドン (aa:Met)
      - CompoundLocation (join / order) の feature

    解釈できない transl_except は警告を出して読み飛ばす。翻訳自体は続ける
    (出力形式を 1 つ増やすたびに record 全体が 500 になるのを避けるため)。

    戻り値は cds=True と同じく、開始コドンを M にし、末尾のストップ記号を落としたもの。
    """
    start_offset = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
    feat_seq: Seq = feature.extract(parent_seq)[start_offset:]  # type: ignore[no-untyped-call]

    # 末端の部分コドン (aa:TERM の transl_except が指すもの) で Seq.translate が
    # 警告を出さないよう、3 の倍数に合わせておく
    remainder = len(feat_seq) % 3
    if remainder:
        feat_seq = feat_seq + Seq("N" * (3 - remainder))

    codon_table_id = feature.qualifiers.get("transl_table", ["Standard"])[0]
    codon_table = _get_codon_table(codon_table_id)

    exceptions: dict[int, str] = {}
    for transl_except in feature.qualifiers.get("transl_except", []):
        try:
            index, aa = _parse_transl_except(transl_except, feature, start_offset)
        except (ValueError, TypeError) as err:
            logger.warning("Ignoring transl_except on feature %s: %s", feature.id, err)
            continue
        exceptions[index] = aa

    # CDS としての検証は自前で行うので cds=False で翻訳する
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        translated = feat_seq.translate(table=codon_table_id, stop_symbol=stop_symbol, to_stop=False, cds=False)  # type: ignore[no-untyped-call]
        protein = list(str(translated))

    for index, aa in exceptions.items():
        if 0 <= index < len(protein):
            protein[index] = aa
        else:
            logger.warning(
                "transl_except codon index %d is out of range (protein length %d) on feature %s; qualifier ignored",
                index, len(protein), feature.id,
            )

    # cds=True と同じく開始コドンを M にする (ATG は元から M なので、それ以外の開始コドンのため)
    first_codon = str(feat_seq[:3]).upper()
    if protein and first_codon in codon_table.start_codons and 0 not in exceptions:
        protein[0] = "M"

    # cds=True と同じく末尾のストップ記号を落とす
    if protein and protein[-1] == stop_symbol:
        protein.pop()

    return Seq("".join(protein))
