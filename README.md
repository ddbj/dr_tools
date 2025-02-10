# mss_tools
Python module for handling MSS format and its json data model

## Installation

```
# installation
pip install "git+https://github.com/ddbj/mss_tools.git"
```

## Usage
```
from mss_tools import mss_ann2json, mss_json2ann, mss_json2gbk, json_to_seqrecords

# DFASTが生成したMSS登録ファイル (ann, seq) をJSONに変換
# DFAST以外のMSS登録ファイルにも今後対応予定
mss_ann2json(ann_file, seq_file, out_json_file)

# JSONファイルをMSS登録ファイル (ann, seq) に変換
# out_dirはデフォルトではカレントディレクトリ
# out_prefixはデフォルトでNoneで、BioSampleやstrainの値を反映して自動で生成される。
# 出力ファイルは {out_dir}/{out_prefix}.ann と .fa
mss_json2ann(json_file, out_dir, out_prefix)


# JSONファイルからBioPythonのSeqRecordに変換
records = json_to_seqrecords(json_file)
# BioPythonの機能を使ってGenBank形式に変換
with open(output_file, "w") as f:
    SeqIO.write(records, f, "genbank")

# JSONファイルと、遺伝子の feature.id を取得して遺伝子詳細情報を辞書として得る
# (DFAST web サービスで遺伝子詳細ページで表示する内容を取得)
from mss_tools.json_utils import get_feature_json
data = get_feature_json(json_file, feature_id)
json_data = json.dumps(data)

# JSONファイルから各種FASTAファイルを生成 (ゲノム、遺伝子塩基配列、タンパク質配列)
# 未実装
# from mss_tools import mss_json2fasta
# mss_json2fasta(json_file, out_fasta_file, format)

```

その他、`mss_tools.MSS.MSS` にMSS登録ファイル情報を格納するクラス、`mss_tools.json_utils` にJSONデータを扱うための関数を定義している。



## Scripts
```
mss_ann2json [-o out.json] input_mss_file.ann input_mss_file.fa
mss_json2ann [-O out_dir] [-o out_prefix] input_file.json
mss_json2gbk [-o out.gbk] input_file.json
```


## Development
```
pip install -e .[dev]
or 
pip install --break-system-packages -e .[dev]
```
