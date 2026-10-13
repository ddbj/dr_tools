import argparse
import sys
from pathlib import Path
from typing import Tuple

from ddbj_record.schema.v3 import Entry, Feature

from dr_tools.json_utils import (features_by_entry, get_entry_id,
                                 get_locus_tag_prefix, json_to_seqfeature,
                                 load_json_to_ddbj_record_instance)

GENE_FEATURE_TYPES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"]


def json_to_gene_seq(feature: Feature, entry: Entry, locus_tag_prefix: str) -> Tuple[str, str, str]:
    """
    record v3 の feature から遺伝子の FASTA の description、塩基配列、タンパク質配列を取得する
    description の書式は "{locus_tag} {product} [type={product}] [location={location}]"
    locus_tag は locus_tag_prefix と feature の locus_tag_id から組み立てる。無ければ feature の alias を使う
    """
    seq_feature, nucleotide_seq = json_to_seqfeature(feature, entry, locus_tag_prefix)
    feature_type = seq_feature.type
    locus_tag = seq_feature.qualifiers.get("locus_tag", [feature.alias or ""])[0]
    product = seq_feature.qualifiers.get("product", ["unknown product"])[0]
    gene = seq_feature.qualifiers.get("gene", [""])[0]
    location = feature.location or "unknown location"
    entry_name = entry.name or get_entry_id(entry)
    if feature_type == "CDS":
        translate = seq_feature.qualifiers.get("translation", [""])[0]
    else:
        translate = ""
    nucleotide = str(nucleotide_seq).upper()
    if gene:
        description = f"{locus_tag} {product} [gene={gene}] [type={feature_type}] [location={entry_name}:{location}]"
    else:
        description = f"{locus_tag} {product} [type={feature_type}] [location={entry_name}:{location}]"
    return description, nucleotide, str(translate)


def json2fasta(json_file: Path, out_dir: Path) -> None:
    """
    JSON ファイルからゲノムおよび遺伝子のFASTAファイルを生成する
    JSON は record v1 / v2 / v3 のどれでもよい (v1 / v2 は v3 に変換してから読む)
    """

    out_dir = Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    genome_fasta, cds_fasta, rna_fasta, protein_fasta = [], [], [], []

    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v3")
    locus_tag_prefix = get_locus_tag_prefix(record)
    grouped_features = features_by_entry(record)

    for entry in (record.sequences.entries if record.sequences is not None else None) or []:
        entry_id = get_entry_id(entry)
        entry_name = entry.name or entry_id
        seq = (entry.sequence or "").upper()
        genome_fasta.append(f">{entry_name}\n{seq}\n")
        for feature in grouped_features.get(entry_id, []):
            if feature.type not in GENE_FEATURE_TYPES:
                continue
            description, nucleotide, translate = json_to_gene_seq(feature, entry, locus_tag_prefix)
            if feature.type == "CDS":
                cds_fasta.append(f">{description}\n{nucleotide}\n")
                protein_fasta.append(f">{description}\n{translate}\n")
            else:
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
