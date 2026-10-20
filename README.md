# DFAST Record Tools (dr_tools)

Python module for handling MSS format and its json data model

## Installation

```bash
# installation
pip install "git+https://github.com/ddbj/dr_tools.git"
```

Requirements:

- Python >= 3.10
- Biopython >= 1.84

## Record の版

JSON は [ddbj-record-specifications](https://github.com/ddbj/ddbj-record-specifications) の DDBJ Record で、v1 / v2 / v3 の 3 つの版がある。

- 読み込み (`drt_json2ann` / `json_to_seqrecords` / `drt_json2fasta` / `drt_json2stats` / `get_feature_json`) は、どの版の JSON も受け付ける。v1 / v2 の JSON は ddbj-record の converter で v3 に変換してから読む (v1 → v2 → v3)
- 書き出し (`drt_ann2json`) は `record_version` で版を選ぶ。既定は v2 で、v3 を指定すると MSS から v3 を直接組み立てる (`provenance.source_format` は `"MSS"`)。v1 / v2 は v1 を組み立ててから converter で変換する
- v3 の JSON を v1 / v2 に戻す converter は無いので、v3 の JSON を v1 / v2 として読む要求はエラーになる

MSS の ann と v3 の対応:

| ann | v3 |
| --- | --- |
| DBLINK `project` / `biosample` / `sequence read archive` | `relations[]` の xref (`db` は `bioproject` / `biosample` / `insdc.sra`) |
| SUBMITTER | `submission.submitters[]` (先頭が連絡先。所属は `organizations[]`) |
| REFERENCE `title` / `ab_name` / `consrtm` / `status` / `year` | `projects[0].publications[]` の `title` / `authors[].abbreviation` / `consortiums` / `status` / `date` |
| COMMENT | `submission.comments[]` (行は改行でつなぐ) |
| ST_COMMENT | `sequences.structured_comments[]` |
| DATE `hold_date` (YYYYMMDD) | `submission.hold_date` (YYYY-MM-DD) |
| KEYWORD | `sequences.keywords[]` |
| DIVISION `division` ENV | `entries[].division` (全ての entry) |
| DATATYPE `type` WGS | `provenance.submission_category` (読むときは `WGS` / `GNM`。書くときは `WGS` と `MAG-WGS` のように末尾が `-WGS` のもの) |
| 全ての entry で同じ source の qualifier | `sequences.common_source` |
| entry 固有の source の qualifier (`/plasmid` など) | `entries[].source_features[].source` (`common_source` を名前ごとに上書きする) |
| source の `ff_definition` | `entries[].source_features[].definition` |
| TOPOLOGY `circular` | `entries[].topology` |
| feature の `locus_tag` | `projects[0].locus_tag_prefix[0].prefix` と `features[].locus_tag_id` |
| 値の無い qualifier (`/pseudo` など) | `value` を省いた qualifier (`{}`) |

`entries[].type` は source の qualifier から決める (`/plasmid` → plasmid、`/chromosome` → chromosome、`/organelle` → organelle、`/segment` → segment、無ければ WGS は unplaced、division が VRL なら segment、それ以外は chromosome)。`entries[].completeness` は `ff_definition` の文から決める (complete / nearly complete / partial。読めなければ省く)。`entries[].division` は COMMON の DIVISION で、DIVISION が無ければ `drt_ann2json` の `division` 引数で、それも省くと DFAST の出力では BCT、それ以外は UNK になる。v3 から ann を書くときは、entry の division が ENV なら DIVISION を書き、BCT や VRL などの分類の division (DDBJ が決める) は書かない。

## Usage

```python
from dr_tools import drt_ann2json, drt_json2ann, json_to_seqrecords, drt_json2fasta, drt_json2stats

# DFAST が生成した MSS 登録ファイル (ann, seq) を DDBJ Record JSON に変換
# record_version は "v1" / "v2" / "v3" で、既定は "v2"
drt_ann2json("examples/complete_genome.ann", "examples/complete_genome.fa", "ddbj_record.json", record_version="v3")

# JSON ファイルを MSS 登録ファイル (ann, fasta) に変換
# out_dir はデフォルトではカレントディレクトリ
# out_prefix はデフォルトで None で、BioSample や strain の値を反映して自動で生成される。
# 出力ファイルは {out_dir}/{out_prefix}.ann と .fasta
drt_json2ann("examples/complete_genome.v3.json", "OUTPUT", "DDBJ-MSS")


# JSON ファイルから BioPython の SeqRecord オブジェクトに変換 (List[SeqRecord])
records = json_to_seqrecords("examples/complete_genome.v3.json")
# BioPython の機能を使って GenBank 形式に変換
from Bio import SeqIO
with open("out.gbk", "w") as f:
    SeqIO.write(records, f, "genbank")


# JSON ファイルから各種 FASTA ファイルを生成 (ゲノム、遺伝子塩基配列、タンパク質配列)
# 出力ファイル名: genome.fna, cds.fna, misc_rnas.fna, protein.faa
drt_json2fasta("ddbj_record.json", "out_dir")

# JSON ファイルからゲノムサイズ、遺伝子数等の統計情報を取得し JSON で保存
# format_=False の場合、数値として保存
# output_fileを指定しない場合、辞書を返す
# keyはDFAST webのUIで表示するときと同じ名称
# GC Content は N を除いた塩基に対する G と C の割合。配列が無いときと、全ての塩基が N のときは 0
drt_json2stats("ddbj_record.json", format_=True, output_file="genome_stats.json")


# JSON ファイルと、遺伝子の feature の alias を指定して遺伝子詳細情報を辞書として得る
# (DFAST web サービスで遺伝子詳細ページで表示する内容を取得)
# v3 の feature の項目に、locus_tag の qualifier、nucleotide、translation を足したもの
from dr_tools.json_utils import get_feature_json
data = get_feature_json("ddbj_record.json", "feature_11")
# 出力例:
print(json.dumps(data, indent=2))
{
  "alias": "feature_11",
  "type": "CDS",
  "location": "1852..2991",
  "sequence_id": "chromosome",
  "qualifiers": {
    "product": [{"value": "DNA polymerase III subunit beta"}],
    "transl_table": [{"value": "11"}],
    "codon_start": [{"value": "1"}],
    "gene": [{"value": "dnaN"}],
    "inference": [{"value": "COORDINATES:ab initio prediction:MetaGeneAnnotator"},
      {"value": "similar to AA sequence:RefSeq:WP_003641629.1"}
    ],
    "locus_tag": [{"value": "PLH_00020"}]
  },
  "locus_tag_id": "00020",
  "nucleotide": "ATGA...",
  "translation": "MKFT..."
}
```

その他、`dr_tools.MSS.MSS` に MSS 登録ファイル情報を格納するクラス (`MSS.parse` が ann / seq を、`MSS.from_json` / `MSS.from_record` が record を読む)、`dr_tools.json_utils` に record v3 を扱うための関数を定義している。`MSS.parse` は、列が 5 つに満たない ann の行を、足りない列を空として読む。seq は `>` で始まる行を配列名にし、`//` だけの行を読み飛ばす。

`examples/` には `complete_genome` (DFAST の完全ゲノム) と `vrl_result` (DFAST_VRL のウイルスゲノム) の ann / fa と、それぞれを v1 (`*.json`)、v2 (`*.v2.json`)、v3 (`*.v3.json`) にした JSON がある。JSON の配列は 10 文字に切り詰めてある。

## Scripts

```bash
drt_ann2json [-o out.json] [-d division] [-r v1|v2|v3] input_mss_file.ann input_mss_file.fa
drt_json2ann [-O out_dir] [-o out_prefix] input_file.json
drt_json2gbk [-o out.gbk] input_file.json
drt_json2fasta [-O out_dir] input_file.json
drt_json2stats_for_dfast [-o out.json] [-f] input_file.json
```

## Development

開発環境は docker compose の `python-dev` である。image は `uv.lock` どおりの依存を `/opt/venv` に入れ、container は host の UID / GID (既定は 1000) で動く。

```bash
docker compose up -d --build
docker compose exec python-dev pytest
docker compose exec python-dev mypy
docker compose exec python-dev ruff check
docker compose exec python-dev ruff format
```

`uv.lock` を変えたら、`docker compose up -d --build` で image を作り直す。

container を使わない場合は、uv で環境を作る。

```bash
uv sync --extra dev
uv run pytest
```

`ddbj-record` を進めるときは `uv lock -P ddbj-record` を実行する。

GitHub Actions (`.github/workflows/test.yml`) は、Python 3.10 から 3.14 のそれぞれで pytest と mypy を実行する。pyproject.toml の下限 (Python 3.10 と biopython 1.84) でも pytest を実行し、ほかに `ruff check` と `ruff format --check` を実行する。下限での pytest は、手元では次のコマンドで実行する。

```bash
uv run --isolated --python 3.10 --extra dev --with biopython==1.84 pytest
```

`tests/fixtures/trad_gnm.json` は ddbj-record-specifications の `tests/fixtures/v3/records/trad_gnm.json` の写しで、`examples/complete_genome.ann` から `drt_ann2json` (v3) が作る record と一致することを確かめるのに使う。
