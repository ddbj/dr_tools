import argparse
import re
import sys
from pathlib import Path

from dr_tools.MSS import MSS

DEFAULT_OUT_PREFIX = "mss"
MAX_OUT_PREFIX_LENGTH = 100


def _sanitize_for_filename(value: str) -> str:
    """
    qualifier の値をファイル名の一部として使える形にする
    英数字とハイフン以外は _ に潰し、連続と前後の _ を落とす
    """
    return re.sub(r"[^A-Za-z0-9-]+", "_", value).strip("_")


def _find_qualifier(mss: MSS, in_common: bool, feature_type: str, key: str) -> str:
    """
    MSS から qualifier の値を 1 つ探す。見つからない場合は空文字列を返す
    in_common が True なら COMMON entry を、False ならそれ以外の entry を見る
    """
    for entry in mss.entries:
        if (entry.id == "COMMON") is not in_common:
            continue
        for feature in entry.features:
            if feature.type != feature_type:
                continue
            for value in feature.qualifiers.get(key, []):
                if isinstance(value, str) and value:
                    return value
    return ""


def generate_out_prefix(mss: MSS) -> str:
    """
    出力ファイル名の prefix を {biosample}_{strain or isolate} の形で組み立てる
    片方しか取れない場合はその値だけを使い、どちらも取れない場合は "mss" を返す
    """
    biosample = _find_qualifier(mss, in_common=True, feature_type="DBLINK", key="biosample")
    strain = _find_qualifier(mss, in_common=False, feature_type="source", key="strain")
    if not strain:
        strain = _find_qualifier(mss, in_common=False, feature_type="source", key="isolate")

    parts = [part for part in (_sanitize_for_filename(biosample), _sanitize_for_filename(strain)) if part]
    prefix = "_".join(parts)[:MAX_OUT_PREFIX_LENGTH].strip("_")

    return prefix or DEFAULT_OUT_PREFIX


def json2ann(mss_json_file: Path, out_dir: Path | str | None = None, out_prefix: str | None = None) -> None:
    """
    jsonファイルからMSSオブジェクトを生成し、そのMSSオブジェクトを使ってMSS登録ファイルを出力する
    out_prefixが指定されていない場合、自動で出力ファイル名を生成する
    {biosample}_{strain or isolate}.ann, {biosample}_{strain or isolate}.fa
    biosampleやstrain or isolateが指定されていない場合は、出力ファイル名はmss.annとmss.faとなる。
    """
    mss = MSS.from_json(mss_json_file)
    if out_dir:
        out_dir = Path(out_dir)
    else:
        out_dir = Path(".")

    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    if out_prefix is None:
        out_prefix = generate_out_prefix(mss)
    out_ann_file = out_dir / Path(f"{out_prefix}.ann")
    out_seq_file = out_dir / Path(f"{out_prefix}.fa")
    mss.write(out_ann_file, out_seq_file)


def main() -> None:
    """
    mss_json_file, out_ann_file, out_seq_file をコマンドライン引数として受け取り、
    jsonファイルからMSSオブジェクトを生成し、そのMSSオブジェクトを使ってMSS登録ファイルを出力する
    """
    # 引数のパース
    parser = argparse.ArgumentParser(description='Convert MSS json file to MSS annotation and sequence files')
    parser.add_argument('json_file', type=str, help='MSS json file')
    parser.add_argument('-O', '--out_dir', type=str, help='Output directory. By default, output files are saved in the current directory.', default=None)
    parser.add_argument('-o', '--out_prefix', type=str, help='Output file prefix. By default, prefix is automatically generated.', default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    json2ann(args.json_file, args.out_dir, args.out_prefix)


if __name__ == "__main__":
    main()
