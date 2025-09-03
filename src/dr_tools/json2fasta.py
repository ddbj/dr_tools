import argparse
import sys
from pathlib import Path
from typing import Any, Dict, Tuple

from dr_tools.json_utils import (get_locus_tag_prefix, json_to_seqfeature,
                                 load_json_to_ddbj_record_instance,
                                 set_locus_tag)


def json_to_gene_seq(feature_json: Dict[str, Any], entry_json: Dict[str, Any]) -> Tuple[str, str, str]:
    """
    Feature データの JSON から遺伝子の FASTA の description、塩基配列、タンパク質配列を取得する
    description の書式は "{locus_tag} {product} [type={product}] [location={location}]"
    実行前には set_locus_tag で locus_tag クオリファイアを設定しておくこと
    """
    seq_feature, nucleotide = json_to_seqfeature(feature_json, entry_json)
    feature_type = seq_feature.type
    locus_tag = seq_feature.qualifiers["locus_tag"][0]
    product = seq_feature.qualifiers.get("product", ["unknown product"])[0]
    gene = seq_feature.qualifiers.get("gene", [""])[0]
    location = feature_json.get("location", "unknown location")
    if feature_type == "CDS":
        translate = seq_feature.qualifiers.get("translation", [""])[0]
    else:
        translate = ""
    nucleotide = str(nucleotide).upper()
    if gene:
        description = f"{locus_tag} {product} [gene={gene}] [type={feature_type}] [location={entry_json['name']}:{location}]"
    else:
        description = f"{locus_tag} {product} [type={feature_type}] [location={entry_json['name']}:{location}]"
    return description, nucleotide, translate


def json2fasta(json_file: Path, out_dir: Path) -> None:
    """
    JSON ファイルからゲノムおよび遺伝子のFASTAファイルを生成する
    """

    out_dir = Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    genome_fasta, cds_fasta, rna_fasta, protein_fasta = [], [], [], []

    # === ddbj record v2 対応 ===
    record_instance = load_json_to_ddbj_record_instance(json_file, to_record_version="v1")
    json_dat = record_instance.model_dump(exclude_none=True, by_alias=True)  # dict形式に変換

    for entry in json_dat.get("ENTRIES", []):
        entry_id = entry["name"]
        seq = entry["sequence"].upper()
        genome_fasta.append(f">{entry_id}\n{seq}\n")
        locus_tag_prefix = get_locus_tag_prefix(json_dat)
        for feature in entry.get("features", []):
            feature_type = feature["type"]
            if feature_type not in ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"]:
                continue
            # feature_id = feature["id"]
            set_locus_tag(feature, locus_tag_prefix)
            description, nucleotide, translate = json_to_gene_seq(feature, entry)
            if feature["type"] == "CDS":
                cds_fasta.append(f">{description}\n{nucleotide}\n")
                protein_fasta.append(f">{description}\n{translate}\n")
            elif feature["type"] in ["rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"]:
                rna_fasta.append(f">{description}\n{nucleotide}\n")

    genome_fasta_file = out_dir / Path("genome.fna")
    cds_fasta_file = out_dir / Path("cds.fna")
    rna_fasta_file = out_dir / Path("misc_rnas.fna")
    protein_fasta_file = out_dir / Path("protein.faa")
    genome_fasta_file.write_text("".join(genome_fasta))
    cds_fasta_file.write_text("".join(cds_fasta))
    rna_fasta_file.write_text("".join(rna_fasta))
    protein_fasta_file.write_text("".join(protein_fasta))


def json2fasta_main() -> None:
    # 引数のパース
    parser = argparse.ArgumentParser(description='Convert MSS json file to genome and gene FASTA files')
    parser.add_argument('json_file', type=str, help='MSS json file')
    parser.add_argument('-O', '--out_dir', type=str, help='Output directory. By default, output files are saved in the current directory.', default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    out_dir = Path(args.out_dir) if args.out_dir else Path(".")
    json2fasta(Path(args.json_file), Path(out_dir))


if __name__ == "__main__":
    json2fasta_main()
