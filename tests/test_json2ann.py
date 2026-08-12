import re
from pathlib import Path
from typing import Dict, List

from hypothesis import example, given
from hypothesis import strategies as st

from dr_tools.json2ann import generate_out_prefix, json2ann
from dr_tools.MSS import MSS, Entry, Feature


def _build_mss(
    *,
    biosample: str | None = None,
    strain: str | None = None,
    isolate: str | None = None,
) -> MSS:
    """COMMON entry の DBLINK/biosample と、別 entry の source/strain, source/isolate を
    任意に持つ MSS を組み立てる。None の場合はその qualifier を持たせない。"""
    common_entry = Entry(id="COMMON", name="COMMON")
    if biosample is not None:
        dblink_qualifiers: Dict[str, List[str | bool]] = {"biosample": [biosample]}
        common_entry.features.append(Feature(type="DBLINK", id="feature_1", qualifiers=dblink_qualifiers))

    source_qualifiers: Dict[str, List[str | bool]] = {}
    if strain is not None:
        source_qualifiers["strain"] = [strain]
    if isolate is not None:
        source_qualifiers["isolate"] = [isolate]
    sample_entry = Entry(
        id="entry1",
        name="entry1",
        features=[Feature(type="source", id="feature_2", location="1..10", qualifiers=source_qualifiers)],
    )

    return MSS(ann_file=None, seq_file=None, entries=[common_entry, sample_entry])


# --- generate_out_prefix: 実データ ---


def test_generate_out_prefix_complete_genome_mss_returns_biosample_and_strain(mss_complete: MSS) -> None:
    assert generate_out_prefix(mss_complete) == "SAMD999999_LOOC260"


def test_generate_out_prefix_vrl_mss_without_biosample_returns_isolate_value(mss_vrl: MSS) -> None:
    assert generate_out_prefix(mss_vrl) == "hCov-19_Japan_SZ-NIG-12345_2021"


# --- generate_out_prefix: どちらも無い / 片方だけ / strain と isolate ---


def test_generate_out_prefix_no_biosample_no_strain_no_isolate_returns_mss() -> None:
    mss = _build_mss()
    assert generate_out_prefix(mss) == "mss"


def test_generate_out_prefix_biosample_only_returns_biosample_value() -> None:
    mss = _build_mss(biosample="SAMD000123")
    assert generate_out_prefix(mss) == "SAMD000123"


def test_generate_out_prefix_strain_only_returns_strain_value() -> None:
    mss = _build_mss(strain="StrainA1")
    assert generate_out_prefix(mss) == "StrainA1"


def test_generate_out_prefix_isolate_only_returns_isolate_value() -> None:
    mss = _build_mss(isolate="IsoB2")
    assert generate_out_prefix(mss) == "IsoB2"


def test_generate_out_prefix_strain_and_isolate_both_present_returns_strain_value() -> None:
    mss = _build_mss(strain="StrainX", isolate="IsolateY")
    assert generate_out_prefix(mss) == "StrainX"


# --- generate_out_prefix: サニタイズ ---


def test_generate_out_prefix_special_characters_are_sanitized_to_underscore() -> None:
    mss = _build_mss(biosample="SAM D/123", strain="A B*C")
    assert generate_out_prefix(mss) == "SAM_D_123_A_B_C"


def test_generate_out_prefix_leading_and_trailing_special_characters_are_stripped() -> None:
    mss = _build_mss(biosample=" SAMD1 ")
    assert generate_out_prefix(mss) == "SAMD1"


def test_generate_out_prefix_each_part_is_stripped_before_join_not_only_the_combined_result() -> None:
    # biosample の末尾と strain の先頭にそれぞれ special character がある場合、
    # 各値を個別にサニタイズしてから join することで、区切りの "_" が二重にならない
    mss = _build_mss(biosample="SAM ", strain=" A1")
    assert generate_out_prefix(mss) == "SAM_A1"


def test_generate_out_prefix_sanitization_removes_all_characters_returns_mss() -> None:
    mss = _build_mss(biosample="!!!", strain="###")
    assert generate_out_prefix(mss) == "mss"


def test_generate_out_prefix_biosample_sanitizes_to_empty_returns_strain_only() -> None:
    mss = _build_mss(biosample="@@@", strain="StrainOK")
    assert generate_out_prefix(mss) == "StrainOK"


# --- generate_out_prefix: 100 文字境界 ---


def test_generate_out_prefix_combined_length_over_100_is_truncated_to_100() -> None:
    mss = _build_mss(strain="A" * 105)
    result = generate_out_prefix(mss)
    assert result == "A" * 100
    assert len(result) == 100


def test_generate_out_prefix_combined_length_exactly_100_is_unchanged() -> None:
    mss = _build_mss(strain="B" * 100)
    result = generate_out_prefix(mss)
    assert result == "B" * 100
    assert len(result) == 100


def test_generate_out_prefix_combined_length_101_is_truncated_to_100() -> None:
    mss = _build_mss(strain="C" * 101)
    result = generate_out_prefix(mss)
    assert result == "C" * 100
    assert len(result) == 100


def test_generate_out_prefix_truncation_drops_trailing_underscore() -> None:
    # biosample を 99 文字にすると、biosample と strain の間の "_" がちょうど 100 文字目に来る
    mss = _build_mss(biosample="D" * 99, strain="E" * 50)
    result = generate_out_prefix(mss)
    assert result == "D" * 99
    assert not result.endswith("_")


# --- generate_out_prefix: PBT ---


@given(
    biosample=st.text(max_size=60),
    strain=st.text(max_size=60),
)
@example(biosample="", strain="")
@example(biosample="!!!###???", strain="   ")
@example(biosample="日本語のビオサンプル", strain="株名/日本語")
@example(biosample="-----", strain="-----")
@example(biosample="A" * 60, strain="B" * 60)
def test_generate_out_prefix_pbt_result_is_always_a_safe_filename(biosample: str, strain: str) -> None:
    mss = _build_mss(biosample=biosample, strain=strain)
    result = generate_out_prefix(mss)
    assert result != ""
    assert len(result) <= 100
    assert re.fullmatch(r"[A-Za-z0-9_-]+", result) is not None


# --- json2ann: out_prefix=None の自動命名 ---


def test_json2ann_complete_v1_out_prefix_none_names_files_by_biosample_and_strain(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json2ann(mss_json_file=eg_complete["v1_json"], out_dir=tmp_path, out_prefix=None)
    assert tmp_path.joinpath("SAMD999999_LOOC260.ann").exists()
    assert tmp_path.joinpath("SAMD999999_LOOC260.fa").exists()


def test_json2ann_complete_v2_out_prefix_none_names_files_by_biosample_and_strain(
    eg_complete: Dict[str, Path], tmp_path: Path
) -> None:
    json2ann(mss_json_file=eg_complete["v2_json"], out_dir=tmp_path, out_prefix=None)
    assert tmp_path.joinpath("SAMD999999_LOOC260.ann").exists()
    assert tmp_path.joinpath("SAMD999999_LOOC260.fa").exists()


def test_json2ann_vrl_v1_out_prefix_none_names_files_by_isolate(eg_vrl: Dict[str, Path], tmp_path: Path) -> None:
    json2ann(mss_json_file=eg_vrl["v1_json"], out_dir=tmp_path, out_prefix=None)
    assert tmp_path.joinpath("hCov-19_Japan_SZ-NIG-12345_2021.ann").exists()
    assert tmp_path.joinpath("hCov-19_Japan_SZ-NIG-12345_2021.fa").exists()


def test_json2ann_vrl_v2_out_prefix_none_names_files_by_isolate(eg_vrl: Dict[str, Path], tmp_path: Path) -> None:
    json2ann(mss_json_file=eg_vrl["v2_json"], out_dir=tmp_path, out_prefix=None)
    assert tmp_path.joinpath("hCov-19_Japan_SZ-NIG-12345_2021.ann").exists()
    assert tmp_path.joinpath("hCov-19_Japan_SZ-NIG-12345_2021.fa").exists()


def test_json2ann_explicit_out_prefix_overrides_auto_naming(eg_complete: Dict[str, Path], tmp_path: Path) -> None:
    json2ann(mss_json_file=eg_complete["v1_json"], out_dir=tmp_path, out_prefix="test_prefix")
    assert tmp_path.joinpath("test_prefix.ann").exists()
    assert tmp_path.joinpath("test_prefix.fa").exists()
    assert not tmp_path.joinpath("SAMD999999_LOOC260.ann").exists()


def test_json2ann_missing_out_dir_is_created(eg_complete: Dict[str, Path], tmp_path: Path) -> None:
    out_dir = tmp_path.joinpath("nested", "sub")
    assert not out_dir.exists()
    json2ann(mss_json_file=eg_complete["v1_json"], out_dir=out_dir, out_prefix="test_prefix")
    assert out_dir.joinpath("test_prefix.ann").exists()
    assert out_dir.joinpath("test_prefix.fa").exists()
