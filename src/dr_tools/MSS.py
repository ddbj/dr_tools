#!/usr/bin/env python3

from dataclasses import dataclass, field
from typing import List, Dict, Optional
import json
from pathlib import Path

BOOL_QUALIFIERS = ["pseudo", "environmental_sample", "ribosomal_slippage", "circular_RNA", "proviral", "focus", "germline", "macronuclear", "circular"]  # circular is only for topology
# the list above may need to be updated. See https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html

@dataclass
class Sequence:
    """
    DDBJアノテーションファイルのFASTAファイルの情報を格納する
    """
    id: str
    seq: str

    def to_fasta(self, width=60, separator=False) -> str:
        """
        FASTA形式の文字列に変換する
        デフォルトでは60文字で改行
        separator=Trueにした場合、DDBJ登録用に // を末尾に加える
        """
        fasta = []
        fasta.append(f">{self.id}")
        for i in range(0, len(self.seq), width):
            fasta.append(self.seq[i:i+width])
        fasta = "\n".join(fasta)
        if separator:
            fasta += "\n//"
        return fasta + "\n"

@dataclass
class Feature:
    """DDBJアノテーションのFeatureを表すクラス

    Attributes:
        type (str): Featureのタイプ（例：'source', 'CDS', 'gene'など）
        id (int | str): Featureの一意な識別子
        location (str): 配列上の位置情報（例：'1..300'）。COMMONタイプの場合は空文字列
        qualifiers (Dict[str, List[str | bool]]): Featureの修飾子（qualifier）を格納する辞書
            - キー：修飾子の名前（例：'product', 'note'など）
            - 値：修飾子の値のリスト。真偽値の修飾子の場合はbool型
    """
    type: str
    id: int | str
    location: str = ""  # COMMONタイプの場合は空文字列
    qualifiers: Dict[str, List[str | bool]] = field(default_factory=dict)

    def to_dict(self):
        """Convert the feature to a dictionary"""
        d = {
            "id": self.id,
            "type": self.type,
            "location": self.location,
            "qualifiers": self.qualifiers
        }
        locus_tag = self.qualifiers.get("locus_tag", [""])[0]
        if locus_tag:
            if "_" not in locus_tag:
                raise ValueError(f"Invalid locus_tag: {locus_tag}")
            del d["qualifiers"]["locus_tag"]
            d["locus_tag_id"] = locus_tag.split("_", 1)[1]
        return d

    def to_tsv(self) -> List[List[str]]:
        # ５列TSV形式に変換
        tsv = []
        for key, values in self.qualifiers.items():
            for value in values:
                if isinstance(value, bool):
                    tsv.append(["", "", "", key, ""])
                else:
                    tsv.append(["", "", "", key, value])
        # location がある場合はlocationを追加 (MSSファイルではlocationを持たないfeatureもある e.g. TOPOLOGYやDDBJ登録に関するもの)
        if self.location:
            if len(tsv) == 0:
                tsv.append(["", "", "", "", ""])  # qualifierを持たないfeatureがある (UTR等)。その場合は空行を追加
            tsv[0][2] = self.location
        tsv[0][1] = self.type
        return tsv

    def show(self):
        """Print the feature in a human-readable format"""
        indent = 2
        print(" " * indent + f"Feature {self.id}:")
        indent += 2
        print(" " * indent + f"Type: {self.type}")
        if self.location:
            print(" " * indent + f"Location: {self.location}")
        for key, values in self.qualifiers.items():
            for value in values:
                print(" " * indent + f"{key}: {value}")

@dataclass
class Entry:
    """
    MSS登録ファイルのEntry (1行目) の情報を格納する
    """
    id: int | str
    name: str = ""
    features: List[Feature] = field(default_factory=list)

    def to_tsv(self) -> List[List[str]]:
        # ５列TSV形式に変換
        tsv = []
        for feature in self.features:
            tsv.extend(feature.to_tsv())
        if tsv:
            tsv[0][0] = self.id
        return tsv

    def show(self):
        """Print the entry in a human-readable format"""
        print(f"Entry ID: {self.id} Name: {self.name}")
        print("Features:")
        for feature in self.features:
            feature.show()

@dataclass
class MSS:
    ann_file: Path
    seq_file: Path
    entries: List[Entry] = field(default_factory=list)
    sequences: List[Sequence] = field(default_factory=list)

    @staticmethod
    def parse(ann_file: Path, seq_file: Path) -> "MSS":
        """
        MSSのannファイルとseqファイルをパースしてMSSインスタンスを作成する       
        """
        current_entry = None
        current_feature: Optional[Feature] = None
        feature_counter = 1  # 連番カウンター
        mss = MSS(ann_file=ann_file, seq_file=seq_file)
        with open(mss.ann_file, "r") as f:
            for line in f:
                line = line.strip("\n")
                if not line:
                    continue
                # ５列に分割. entry, feature, location, qualifier_key, qualifier_value
                col_entry, col_feature, col_location, col_qkey, col_qvalue = line.split("\t", 4)
                
                if col_entry:
                    # 列1に値を持つ場合は新しいEntryを作成
                    current_entry = Entry(id=col_entry, name=col_entry)
                    mss.entries.append(current_entry)

                if col_feature:
                    # 列2に値を持つ場合は新しいFeatureを作成
                    current_feature = Feature(type=col_feature, id=f"feature_{feature_counter}")
                    feature_counter += 1
                    current_entry.features.append(current_feature)

                if col_location:
                    # 列3に値を持つ場合はcurrent_featureのlocationを設定
                    current_feature.location = col_location

                if col_qkey:
                    # current_featureにqualifierを追加
                    if col_qkey in BOOL_QUALIFIERS and col_qvalue == "":
                        # 値が空の場合はTrueを追加
                        current_feature.qualifiers.setdefault(col_qkey, []).append(True)
                    else:
                        current_feature.qualifiers.setdefault(col_qkey, []).append(col_qvalue)
        mss.parse_seq()
        return mss

    def parse_seq(self) -> None:
        """Parse the DDBJ sequence file"""
        current_seq = None
        
        with open(self.seq_file, "r") as f:
            entries = f.read()
            entries = entries.replace("//\n", "") # remove trailing slash (// in DDBJ format)
            entries = entries.split(">")
            entries = entries[1:]  # remove the first empty entry
            for entry in entries:
                lines = entry.split("\n")
                seq_id = lines[0]
                seq = "".join(lines[1:])
                self.sequences.append(Sequence(id=seq_id, seq=seq))



    def to_tsv(self) -> str:
        """Convert the parsed data to TSV format"""
        tsv = []
        for entry in self.entries:
            tsv.extend(entry.to_tsv())
        return "\n".join(["\t".join(row) for row in tsv]) + "\n"

    def to_fasta(self, width=60, separator=False) -> str:
        """Convert the parsed sequence data to FASTA format"""
        fasta = []
        for sequence in self.sequences:
            fasta.append(sequence.to_fasta(width=width, separator=separator))
        return "".join(fasta)
    
    def write(self, out_ann_file: Path, out_seq_file: Path) -> None:
        with open(out_ann_file, "w") as f:
            f.write(self.to_tsv())
        with open(out_seq_file, "w") as f:
            f.write(self.to_fasta(separator=True))

    def show(self) -> None:
        """[For debugging] Print the parsed data in a human-readable format"""
        print("\nEntries:")
        for entry in self.entries:
            entry.show()
        
        print("\nSequences:")
        for seq_id, sequence in self.sequences.items():
            print(f"\n{seq_id}:")
            print(f"  Length: {len(sequence)}")
            if len(sequence) > 10:
                print(f"  Sequence: {sequence[:10]}...")
            else:
                print(f"  Sequence: {sequence}")

    @staticmethod
    def from_json(json_file: Path) -> "MSS":
        """
        JSONファイルを読み込みMSSインスタンスを作成する
        """

        with open(json_file) as f:
            data = json.load(f)

        mss = MSS(ann_file="", seq_file="")  # 空のMSSインスタンスを作成
        
        # COMMON エントリーの作成
        common_entry = Entry(id="COMMON", name="COMMON")
        feature_counter = 1
        
        # DBLINKの作成
        if "DBLINK" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="DBLINK"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["DBLINK"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # SUBMITTERの作成
        if "SUBMITTER" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="SUBMITTER"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["SUBMITTER"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # REFERENCEの作成
        if "REFERENCE" in data["COMMON"]:
            for ref in data["COMMON"]["REFERENCE"]:
                feature = Feature(
                    id=f"feature_{feature_counter}",
                    type="REFERENCE"
                )
                feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                    for k, v in ref.items()}
                common_entry.features.append(feature)
                feature_counter += 1
        
        # COMMENTの作成
        if "COMMENT" in data["COMMON"]:
            for comment in data["COMMON"]["COMMENT"]:
                feature = Feature(
                    id=f"feature_{feature_counter}",
                    type="COMMENT"
                )
                feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                    for k, v in comment.items()}
                common_entry.features.append(feature)
                feature_counter += 1
        
        # ST_COMMENTの作成
        if "ST_COMMENT" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="ST_COMMENT"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["ST_COMMENT"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # DATEの作成
        if "DATE" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="DATE"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["DATE"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # DATATYPEの作成
        if "DATATYPE" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="DATATYPE"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["DATATYPE"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # KEYWORDの作成
        if "KEYWORD" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="KEYWORD"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["KEYWORD"].items()}
            common_entry.features.append(feature)
        
        # common_entryを先頭に追加
        mss.entries.append(common_entry)
        
        # COMMON_SOURCEの内容を取得
        common_source = data.get("COMMON_SOURCE", {})
        
        # エントリーの作成
        for entry_data in data["ENTRIES"]:
            entry_id = entry_data["id"]
            entry = Entry(id=entry_id, name=entry_data["name"])
            features = []
            
            # source featureを作成し、COMMON_SOURCEの内容を追加
            for feature_data in entry_data["features"]:
                if feature_data["type"] == "source":
                    feature_id = feature_data["id"]
                    feature = Feature(
                        id=feature_id,
                        type="source",
                        location=feature_data["location"]
                    )
                    # COMMON_SOURCEの内容を追加
                    qualifiers = dict(common_source)
                    # エントリー固有のqualifierを追加
                    qualifiers.update({k: v for k, v in feature_data["qualifiers"].items()})
                    feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                       for k, v in qualifiers.items()}
                    features.append(feature)
                    break
            
            # topologyがcircularの場合、TOPOLOGY featureを追加
            if entry_data.get("_topology") == "circular":
                topology_feature = Feature(
                    id=f"topology_{entry_id}",
                    type="TOPOLOGY",
                    location=""
                )
                topology_feature.qualifiers = {"circular": [True]}
                features.append(topology_feature)
            
            # その他のfeatureを追加
            for feature_data in entry_data["features"]:
                feature_id = feature_data["id"]
                if feature_data["type"] != "source":  # sourceは既に追加済み
                    feature = Feature(
                        id=feature_id,
                        type=feature_data["type"],
                        location=feature_data["location"]
                    )
                    feature.qualifiers = feature_data["qualifiers"]
                    features.append(feature)
            
            # 作成したfeaturesをentryに設定
            entry.features = features
            
            # シーケンスの作成
            sequence = Sequence(id=entry_id, seq=entry_data["sequence"])
            mss.sequences.append(sequence)
            mss.entries.append(entry)
        
        return mss

def test():

    input_ann, input_seq = "examples/complete_genome.ann", "examples/complete_genome.fa"
    output_ann, output_seq = "complete_genome_test.ann", "complete_genome_test.fa"
    # MSSファイルをパース
    mss = MSS.parse(input_ann, input_seq)
    mss.write(output_ann, output_seq)
    # print(mss.to_fasta(separator=True), end="")
    # 結果を表示
    
    # 内容が一致することを確認
    with open(input_ann, "r") as f_in, open(output_ann, "r") as f_out:
        assert f_out.read() == f_in.read()

    # テストファイルを削除
    import os
    os.remove("complete_genome_test.ann")
    os.remove("complete_genome_test.fa")

if __name__ == "__main__":
    test() 