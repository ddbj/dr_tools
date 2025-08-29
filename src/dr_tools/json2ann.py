import argparse
import sys
from pathlib import Path

from dr_tools.MSS import MSS


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
        out_ann_file = out_dir / Path("mss.ann")
        out_seq_file = out_dir / Path("mss.fa")
    else:
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
