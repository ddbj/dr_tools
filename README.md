# mss_tools
Python module for handling MSS format and its json data model

## Installation

```
# installation
pip install "git+https://github.com/ddbj/mss_tools.git"
```

Requirements:  
- Python >= 3.9
- Biopython >= 1.82


## Usage
```
from mss_tools import mss_ann2json, mss_json2ann, json_to_seqrecords, mss_json2fasta, mss_json2stats

# DFAST が生成した MSS 登録ファイル (ann, seq) を DFAST results JSON に変換
# DFAST 以外の MSS 登録ファイルにも今後対応予定
mss_ann2json("examples/complete_genome.ann", "examples/complete_genome.fa", "dfast_results.json")

# JSON ファイルを MSS 登録ファイル (ann, seq) に変換
# out_dir はデフォルトではカレントディレクトリ
# out_prefix はデフォルトで None で、BioSample や strain の値を反映して自動で生成される。
# 出力ファイルは {out_dir}/{out_prefix}.ann と .fa
mss_json2ann("examples/complete_genome.json", "OUTPUT", "DDBJ-MSS")


# JSON ァイルから BioPython の SeqRecord オブジェクトに変換 (List[SeqRecord])
records = json_to_seqrecords("examples/complete_genome.json")
# BioPython の機能を使って GenBank 形式に変換
from Bio import SeqIO
with open("out.gbk", "w") as f:
    SeqIO.write(records, f, "genbank")


# JSON ファイルから各種 FASTA ファイルを生成 (ゲノム、遺伝子塩基配列、タンパク質配列)
# 出力ファイル名: genome.fna, cds.fna, misc_rnas.fna, protein.faa
mss_json2fasta("dfast_results.json", "out_dir")

# JSON ファイルからゲノムサイズ、遺伝子数等の統計情報を取得し JSON で保存
# format=Falaseの場合、数値として保存
# output_fileを指定しない場合、辞書を返す
# keyはDFAST webのUIで表示するときと同じ名称
mss_json2stats("complete_genome.json", format=True, output_file="genome_stats.json")


# JSON ファイルと、遺伝子の feature.id を取得して遺伝子詳細情報を辞書として得る
# (DFAST web サービスで遺伝子詳細ページで表示する内容を取得)
from mss_tools.json_utils import get_feature_json
data = get_feature_json("dfast_results.json", "feature_11")
# 出力例:
print(json.dumps(data, indent=2))
{
  "id": "feature_11",
  "type": "CDS",
  "location": "1852..2991",
  "qualifiers": {
    "product": ["DNA polymerase III subunit beta"],
    "transl_table": ["11"],
    "codon_start": ["1"],
    "gene": ["dnaN"],
    "inference": ["COORDINATES:ab initio prediction:MetaGeneAnnotator",
      "similar to AA sequence:RefSeq:WP_003641629.1"
    ],
    "locus_tag": ["PLH_00020"]
  },
  "locus_tag_id": "00020",
  "nucleotide": "ATGA...",
  "translation": "MKFT..."
}
```

その他、`mss_tools.MSS.MSS` に MSS 登録ファイル情報を格納するクラス、`mss_tools.json_utils` に JSON データを扱うための関数を定義している。



## Scripts
```
mss_ann2json [-o out.json] input_mss_file.ann input_mss_file.fa
mss_json2ann [-O out_dir] [-o out_prefix] input_file.json
mss_json2gbk [-o out.gbk] input_file.json
mss_json2fasta [-O out_dir] input_file.json
mss_json2stats_for_dfast [-o out.json] [-f] input_file.json
```


## Development
```
pip install -e .[dev]
or 
pip install --break-system-packages -e .[dev]
```
