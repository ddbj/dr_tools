from pathlib import Path
from typing import Dict, Tuple

import pytest

from dr_tools.ann2json import ann2json_for_dfast
from dr_tools.MSS import MSS

HERE = Path(__file__).parent.resolve()
REPO_ROOT = HERE.parent
EXAMPLES_DIR = REPO_ROOT.joinpath("examples")
EG_COMPLETE = {
    "ann": EXAMPLES_DIR.joinpath("complete_genome.ann"),
    "seq": EXAMPLES_DIR.joinpath("complete_genome.fa"),
    "v1_json": EXAMPLES_DIR.joinpath("complete_genome.json"),
    "v2_json": EXAMPLES_DIR.joinpath("complete_genome.v2.json"),
}
EG_VRL = {
    "ann": EXAMPLES_DIR.joinpath("vrl_result.ann"),
    "seq": EXAMPLES_DIR.joinpath("vrl_result.fa"),
    "v1_json": EXAMPLES_DIR.joinpath("vrl_result.json"),
    "v2_json": EXAMPLES_DIR.joinpath("vrl_result.v2.json"),
}

# WGS (draft) の MSS。examples には GNM しかないので、登録区分で挙動が変わる箇所のために置く。
# contig01 は circular、contig02 は linear。どちらの配列も 30 bp で、CDS は 1 つ。
WGS_ANN = "\n".join([
    "COMMON\tDATATYPE\t\ttype\tWGS",
    "\tKEYWORD\t\tkeyword\tWGS",
    "\tSUBMITTER\t\tab_name\tMishima,H.",
    "\t\t\tcontact\tHanako Mishima",
    "\t\t\temail\tmishima@ddbj.nig.ac.jp",
    "\t\t\tinstitute\tNational Institute of Genetics",
    "\t\t\tcountry\tJapan",
    "\t\t\tcity\tMishima",
    "\t\t\tstreet\tYata 1111",
    "\t\t\tzip\t411-8540",
    "\tREFERENCE\t\ttitle\tDraft genome of a test organism",
    "\t\t\tab_name\tMishima,H.",
    "\t\t\tstatus\tUnpublished",
    "\t\t\tyear\t2026",
    "\tCOMMENT\t\tline\tAnnotated by DFAST 1.4.2 https://dfast.ddbj.nig.ac.jp/",
    "contig01\tsource\t1..30\torganism\tEscherichia coli",
    "\t\t\tmol_type\tgenomic DNA",
    "\t\t\tff_definition\t@@[organism]@@ DNA, @@[submitter_seqid]@@",
    "\tTOPOLOGY\t\tcircular\t",
    "\tCDS\t1..30\tproduct\thypothetical protein",
    "\t\t\ttransl_table\t11",
    "\t\t\tcodon_start\t1",
    "\t\t\tlocus_tag\tTEST_00010",
    "contig02\tsource\t1..30\torganism\tEscherichia coli",
    "\t\t\tmol_type\tgenomic DNA",
    "\t\t\tff_definition\t@@[organism]@@ DNA, @@[submitter_seqid]@@",
    "\tCDS\t1..30\tproduct\thypothetical protein",
    "\t\t\ttransl_table\t11",
    "\t\t\tcodon_start\t1",
    "\t\t\tlocus_tag\tTEST_00020",
]) + "\n"
WGS_SEQ = "\n".join([
    ">contig01",
    "atgaaacgcattagcaccaccattaccacc",
    "//",
    ">contig02",
    "atgaaacgcattagcaccaccattaccacc",
    "//",
]) + "\n"


@pytest.fixture
def eg_complete() -> Dict[str, Path]:
    return EG_COMPLETE


@pytest.fixture
def eg_vrl() -> Dict[str, Path]:
    return EG_VRL


@pytest.fixture
def mss_complete() -> MSS:
    return MSS.parse(EG_COMPLETE["ann"], EG_COMPLETE["seq"])


@pytest.fixture
def mss_vrl() -> MSS:
    return MSS.parse(EG_VRL["ann"], EG_VRL["seq"])


@pytest.fixture(scope="session")
def full_record(tmp_path_factory: pytest.TempPathFactory) -> Dict[str, Path]:
    """
    完全な配列を持つ record JSON を ann / fa から作って返す

    examples/complete_genome.json は配列を 13 文字に切り詰めた確認用のデータなので、
    翻訳や配列の長さを見るテストの入力には使えない。
    """
    out_dir = tmp_path_factory.mktemp("full_record")
    result = {}
    for version in ("v1", "v2"):
        json_file = out_dir.joinpath(f"complete_genome.{version}.json")
        ann2json_for_dfast(EG_COMPLETE["ann"], EG_COMPLETE["seq"], json_file, record_version=version)
        result[version] = json_file

    return result


@pytest.fixture
def wgs_mss_files(tmp_path: Path) -> Tuple[Path, Path]:
    """WGS (draft) の ann / fa を書き出してそのパスを返す"""
    ann_file = tmp_path.joinpath("wgs.ann")
    seq_file = tmp_path.joinpath("wgs.fa")
    ann_file.write_text(WGS_ANN, encoding="utf-8")
    seq_file.write_text(WGS_SEQ, encoding="utf-8")

    return ann_file, seq_file
