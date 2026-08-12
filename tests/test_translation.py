import logging

import pytest
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from dr_tools.json2biopython import create_seqfeature
from dr_tools.translation import translate_cds_with_transl_except

# 共通のベース配列 (30 bp, 10 codon)。ATG AAA CGC ATT AGC ACC ACC ATT ACC ACC
# 標準遺伝暗号表での翻訳は手計算で "MKRISTTITT" (M K R I S T T I T T)。
BASE_SEQ = "ATGAAACGCATTAGCACCACCATTACCACC"
BASE_PROTEIN = "MKRISTTITT"


def _make_cds(location: str, seq_length: int, qualifiers: dict[str, list[str]] | None = None) -> SeqFeature:
    return create_seqfeature("CDS", location, qualifiers or {}, seq_length, "linear")


def test_translate_cds_with_transl_except_no_exceptions_translates_directly() -> None:
    feature = _make_cds("1..30", 30)
    result = translate_cds_with_transl_except(feature, Seq(BASE_SEQ))
    assert str(result) == BASE_PROTEIN


def test_translate_cds_with_transl_except_alt_start_codon_table11_forces_m() -> None:
    # GTG は table 11 (Bacterial) の開始コドンだが、Val (V) をコードするコドンでもある。
    # 先頭に来たときだけ M に強制される。
    seq = "GTG" + BASE_SEQ[3:]
    feature = _make_cds("1..30", 30, {"transl_table": ["11"]})
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == "M" + BASE_PROTEIN[1:]


def test_translate_cds_with_transl_except_non_start_codon_not_forced_to_m() -> None:
    # CCC (Pro) は table 1 の start_codons に含まれないので、先頭でも M にならない。
    seq = "CCC" + BASE_SEQ[3:]
    feature = _make_cds("1..30", 30)
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == "P" + BASE_PROTEIN[1:]


def test_translate_cds_with_transl_except_codon_start_shifts_frame() -> None:
    # 先頭に 1 塩基足し、codon_start=2 でそれを読み飛ばすとベース配列と同じフレームになる。
    seq = "C" + BASE_SEQ
    feature = _make_cds("1..31", 31, {"codon_start": ["2"]})
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == BASE_PROTEIN


def test_translate_cds_with_transl_except_first_position_exception_overrides_m_forcing() -> None:
    # 先頭コドン (GTG, table 11) は M 強制の対象だが、そこに transl_except があれば
    # そちらが優先され、M 強制は行われない。
    seq = "GTG" + BASE_SEQ[3:]
    feature = _make_cds(
        "1..30", 30,
        {"transl_table": ["11"], "transl_except": ["(pos:1..3,aa:OTHER)"]},
    )
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == "X" + BASE_PROTEIN[1:]


def test_translate_cds_with_transl_except_single_position_partial_stop_term() -> None:
    # feature の長さが 28 bp (9 codon + 1 bp) で、その 1 bp を末端の部分ストップとして指定する。
    # 28 bp 目だけを指す (範囲でない) transl_except を解釈できることを確認する。
    seq = BASE_SEQ[:28]
    feature = _make_cds("1..28", 28, {"transl_except": ["(pos:28,aa:TERM)"]})
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == BASE_PROTEIN[:9]


def test_translate_cds_with_transl_except_range_position_partial_stop_term_capitalized() -> None:
    # 29 bp (9 codon + 2 bp) の部分ストップを範囲指定で表す。aa 名は "Term" 表記を使う。
    seq = BASE_SEQ[:29]
    feature = _make_cds("1..29", 29, {"transl_except": ["(pos:28..29,aa:Term)"]})
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == BASE_PROTEIN[:9]


def test_translate_cds_with_transl_except_compound_location_join_applies_other_exception() -> None:
    # join(1..9,20..40): 1-9 bp と 20-40 bp を連結すると BASE_SEQ と同じ 30 bp になる
    # (10-19 bp は feature に含まれないギャップ)。2 番目の part の先頭コドン
    # (親配列上 20-22 bp = codon index 3 = "ATT") を aa:OTHER で X に差し替える。
    parent = Seq(BASE_SEQ[:9] + "N" * 10 + BASE_SEQ[9:])
    feature = _make_cds("join(1..9,20..40)", len(parent), {"transl_except": ["(pos:20..22,aa:OTHER)"]})
    result = translate_cds_with_transl_except(feature, parent)
    assert str(result) == BASE_PROTEIN[:3] + "X" + BASE_PROTEIN[4:]


def test_translate_cds_with_transl_except_complement_location_sec_maps_to_correct_codon() -> None:
    # feature 全体が - 鎖 (complement(1..30))。親配列は BASE_SEQ の逆相補なので、
    # extract() した結果は BASE_SEQ と一致する。
    # codon index 1 ("AAA", 親配列上では相補鎖の 25-27 bp) を aa:Sec で U に差し替える。
    parent = Seq(BASE_SEQ).reverse_complement()  # type: ignore[no-untyped-call]
    feature = _make_cds(
        "complement(1..30)", 30,
        {"transl_except": ["(pos:complement(25..27),aa:Sec)"]},
    )
    result = translate_cds_with_transl_except(feature, parent)
    assert str(result) == BASE_PROTEIN[0] + "U" + BASE_PROTEIN[2:]


def test_translate_cds_with_transl_except_multiple_exceptions_all_applied_without_raising() -> None:
    # transl_except を 3 個 (Sec, Pyl, TERM) 同時に持つ feature でも、例外を投げずに
    # 全部が反映される。末尾を TERM にした場合、末尾のストップ記号は 1 つ落ちる。
    feature = _make_cds(
        "1..30", 30,
        {"transl_except": [
            "(pos:4..6,aa:Sec)",
            "(pos:25..27,aa:Pyl)",
            "(pos:28..30,aa:TERM)",
        ]},
    )
    result = translate_cds_with_transl_except(feature, Seq(BASE_SEQ))
    assert str(result) == "MURISTTIO"


def test_translate_cds_with_transl_except_unparseable_value_warns_and_is_ignored(
    caplog: pytest.LogCaptureFixture,
) -> None:
    feature = _make_cds("1..30", 30, {"transl_except": ["not a valid transl_except"]})
    with caplog.at_level(logging.WARNING, logger="dr_tools.translation"):
        result = translate_cds_with_transl_except(feature, Seq(BASE_SEQ))
    assert str(result) == BASE_PROTEIN
    assert "Ignoring transl_except" in caplog.text


def test_translate_cds_with_transl_except_unknown_amino_acid_warns_and_is_ignored(
    caplog: pytest.LogCaptureFixture,
) -> None:
    feature = _make_cds("1..30", 30, {"transl_except": ["(pos:1..3,aa:Foo)"]})
    with caplog.at_level(logging.WARNING, logger="dr_tools.translation"):
        result = translate_cds_with_transl_except(feature, Seq(BASE_SEQ))
    assert str(result) == BASE_PROTEIN
    assert "Unknown amino acid" in caplog.text


def test_translate_cds_with_transl_except_out_of_range_position_warns_and_is_ignored(
    caplog: pytest.LogCaptureFixture,
) -> None:
    feature = _make_cds("1..30", 30, {"transl_except": ["(pos:1000..1002,aa:TERM)"]})
    with caplog.at_level(logging.WARNING, logger="dr_tools.translation"):
        result = translate_cds_with_transl_except(feature, Seq(BASE_SEQ))
    assert str(result) == BASE_PROTEIN
    assert "out of range" in caplog.text


def test_translate_cds_with_transl_except_non_multiple_of_three_length_pads_with_n() -> None:
    # 31 bp (10 codon + 1 bp)。末尾を N で 2 つ埋めて 11 番目の codon "ANN" を作る。
    # ANN は標準遺伝暗号表では確定しないアミノ酸になるため X になる。
    seq = BASE_SEQ + "A"
    feature = _make_cds("1..31", 31)
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == BASE_PROTEIN + "X"


def test_translate_cds_with_transl_except_trailing_stop_symbol_dropped_once() -> None:
    # "ATG AAA CGC TAA TAA" は連続する 2 つのストップコドンを含む。
    # 落ちるのは末尾の 1 つだけで、内側のストップ記号は残る。
    seq = "ATGAAACGCTAATAA"
    feature = _make_cds("1..15", 15)
    result = translate_cds_with_transl_except(feature, Seq(seq))
    assert str(result) == "MKR*"


def test_translate_cds_with_transl_except_non_nucleotide_characters_raise_translation_error() -> None:
    # codon index 1 を不正な文字にすると、Bio.Seq.translate が
    # Bio.Data.CodonTable.TranslationError を送出する。呼び出し側で握りつぶさない。
    seq = "ATG" + "ZZZ" + BASE_SEQ[6:]
    feature = _make_cds("1..30", 30)
    with pytest.raises(TranslationError):
        translate_cds_with_transl_except(feature, Seq(seq))
