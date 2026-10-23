import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ddbj_record.schema.v3 import DdbjRecord as DdbjRecordV3
from ddbj_record.schema.v3 import Person, Publication, Relation, StructuredComment

from dr_tools.json_utils import (
    QualifierValues,
    feature_qualifier_values,
    features_by_entry,
    get_entry_id,
    get_locus_tag_prefix,
    load_json_to_ddbj_record_instance,
    merge_source_qualifiers,
)

# COMMON の DIVISION に書く division。BCT や VRL などの分類の division は DDBJ が決めるので書かない。
# MAG (metagenome-assembled genome) の entry は ENV になる。
MSS_DIVISIONS = ("ENV",)

# 値を持たない qualifier の名前。ann では値なしの qualifier と空文字列の qualifier を区別できないので、
# ann を読むときはこの名前のときだけ True (値なし) にする。record v3 から読むときは value の無い
# qualifier を全て True にするので、この一覧は使わない。
BOOL_QUALIFIERS = [
    "pseudo",
    "environmental_sample",
    "ribosomal_slippage",
    "circular_RNA",
    "proviral",
    "focus",
    "germline",
    "macronuclear",
    "circular",
]  # circular is only for topology
# the list above may need to be updated. See https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html

# DBLINK の qualifier と、record の relations (xref) の db の対応
DBLINK_KEY_TO_DB = {
    "project": "bioproject",
    "biosample": "biosample",
    "sequence read archive": "insdc.sra",
}
DB_TO_DBLINK_KEY = {db: key for key, db in DBLINK_KEY_TO_DB.items()}

# SUBMITTER の qualifier の並び。record から ann を書くときにこの順で並べる
SUBMITTER_KEYS = [
    "ab_name",
    "contact",
    "email",
    "url",
    "institute",
    "department",
    "consrtm",
    "country",
    "state",
    "city",
    "street",
    "zip",
]

# show() で表示する配列の先頭の長さ
SHOW_SEQUENCE_LENGTH = 10

_YYYYMMDD = re.compile(r"\d{8}")
_YYYY_MM_DD = re.compile(r"\d{4}-\d{2}-\d{2}")


def hold_date_to_record(hold_date: str) -> str:
    """DATE の hold_date (YYYYMMDD) を record の hold_date (YYYY-MM-DD) にする。その形でなければそのまま返す"""
    if _YYYYMMDD.fullmatch(hold_date):
        return f"{hold_date[:4]}-{hold_date[4:6]}-{hold_date[6:8]}"
    return hold_date


def hold_date_to_mss(hold_date: str) -> str:
    """record の hold_date (YYYY-MM-DD) を DATE の hold_date (YYYYMMDD) にする。その形でなければそのまま返す"""
    if _YYYY_MM_DD.fullmatch(hold_date):
        return hold_date.replace("-", "")
    return hold_date


def reference_status_to_record(status: str) -> str:
    """REFERENCE の status ("In press" など) を record の書き方 ("in-press") にする"""
    return "-".join(status.lower().split(" "))


def reference_status_to_mss(status: str) -> str:
    """record の status ("in-press" など) を REFERENCE の書き方 ("In press") にする"""
    words = " ".join(status.split("-"))
    return words[:1].upper() + words[1:]


@dataclass
class Sequence:
    """
    DDBJアノテーションファイルのFASTAファイルの情報を格納する
    """

    id: str
    seq: str

    def to_fasta(self, width: int = 60, separator: bool = False) -> str:
        """
        FASTA形式の文字列に変換する
        デフォルトでは60文字で改行
        separator=Trueにした場合、DDBJ登録用に // を末尾に加える
        """
        fasta = [f">{self.id}"]
        fasta.extend(self.seq[i : i + width] for i in range(0, len(self.seq), width))
        fasta_str = "\n".join(fasta)
        if separator:
            fasta_str += "\n//"
        return fasta_str + "\n"


@dataclass
class Feature:
    """DDBJアノテーションのFeatureを表すクラス

    Attributes:
        type (str): Featureのタイプ（例：'source', 'CDS', 'gene'など）
        id (int | str): Featureの一意な識別子
        location (str): 配列上の位置情報（例：'1..300'）。COMMONタイプの場合は空文字列
        qualifiers (Dict[str, List[str | bool]]): Featureの修飾子（qualifier）を格納する辞書
            - キー：修飾子の名前（例：'product', 'note'など）
            - 値：修飾子の値のリスト。値を持たない修飾子の場合は True
    """

    type: str
    id: int | str
    location: str = ""  # COMMONタイプの場合は空文字列
    qualifiers: dict[str, list[str | bool]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the feature to a dictionary
        locus_tag は qualifiers から外して locus_tag_id にする。self.qualifiers は変えない
        """
        qualifiers = dict(self.qualifiers)
        d: dict[str, Any] = {"id": self.id, "type": self.type, "location": self.location, "qualifiers": qualifiers}
        locus_tag = qualifiers.get("locus_tag", [""])[0]
        if isinstance(locus_tag, str) and locus_tag:
            if "_" not in locus_tag:
                raise ValueError(f"Invalid locus_tag: {locus_tag}")
            del qualifiers["locus_tag"]
            d["locus_tag_id"] = locus_tag.split("_", 1)[1]
        return d

    def to_tsv(self) -> list[list[str]]:
        """
        ５列TSV形式に変換する

        値を持たない qualifier の扱いは 2 通りある。
        - 値を持たない qualifier (/pseudo や TOPOLOGY の /circular。True で表す) は、値の列を空にした行として出す
        - 値が空文字列の qualifier は行ごと落とす。MSS では値の無い qualifier は
          上の True 側だけを指すので、空文字列の行を出すと DDBJ 側で値なしの
          qualifier として解釈されてしまう
        """
        tsv = []
        for key, values in self.qualifiers.items():
            for value in values:
                if isinstance(value, bool):
                    tsv.append(["", "", "", key, ""])
                elif value != "":
                    tsv.append(["", "", "", key, value])
        if len(tsv) == 0:
            # 出力する qualifier を持たない feature がある (UTR 等)。location 用に空行を置く
            tsv.append(["", "", "", "", ""])
        # location を持たない feature もある (TOPOLOGY や DDBJ 登録に関するもの)
        if self.location:
            tsv[0][2] = self.location
        tsv[0][1] = self.type
        return tsv

    def show(self) -> None:
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

    id は record JSON 上の識別子、name は利用者が付けた配列名。
    MSS 登録ファイルと FASTA に書き出す識別子は name に統一する
    (mss.ann の1列目、mss.fasta のヘッダ、genome.fna のヘッダ)。
    """

    id: str
    name: str
    features: list[Feature] = field(default_factory=list)

    def to_tsv(self) -> list[list[str]]:
        # ５列TSV形式に変換
        tsv = []
        for feature in self.features:
            tsv.extend(feature.to_tsv())
        if tsv:
            tsv[0][0] = self.name
        return tsv

    def show(self) -> None:
        """Print the entry in a human-readable format"""
        print(f"Entry ID: {self.id} Name: {self.name}")
        print("Features:")
        for feature in self.features:
            feature.show()


@dataclass
class MSS:
    ann_file: Path | None
    seq_file: Path | None
    entries: list[Entry] = field(default_factory=list)
    sequences: list[Sequence] = field(default_factory=list)

    @staticmethod
    def parse(ann_file: Path, seq_file: Path) -> "MSS":
        """
        MSSのannファイルとseqファイルをパースしてMSSインスタンスを作成する
        """
        current_entry = None
        current_feature: Feature | None = None
        feature_counter = 1  # 連番カウンター
        mss = MSS(ann_file=ann_file, seq_file=seq_file)
        if mss.ann_file is None or mss.seq_file is None:
            raise ValueError("MSS instance must have ann_file and seq_file")

        with open(mss.ann_file, encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.strip("\n")
                if not line:
                    continue
                # ５列に分割. entry, feature, location, qualifier_key, qualifier_value
                # 末尾の空の列のタブが省かれた行もあるので、足りない列は空文字列で埋める
                columns = line.split("\t", 4)
                columns += [""] * (5 - len(columns))
                col_entry, col_feature, col_location, col_qkey, col_qvalue = columns

                if col_entry:
                    # 列1に値を持つ場合は新しいEntryを作成
                    current_entry = Entry(id=col_entry, name=col_entry)
                    mss.entries.append(current_entry)

                if col_feature:
                    # 列2に値を持つ場合は新しいFeatureを作成
                    current_feature = Feature(type=col_feature, id=f"feature_{feature_counter}")
                    feature_counter += 1
                    if current_entry is not None:
                        current_entry.features.append(current_feature)

                if col_location and current_feature is not None:
                    # 列3に値を持つ場合はcurrent_featureのlocationを設定
                    current_feature.location = col_location

                if col_qkey and current_feature is not None:
                    # current_featureにqualifierを追加
                    if col_qkey in BOOL_QUALIFIERS and col_qvalue == "":
                        # 値が空の場合はTrueを追加
                        current_feature.qualifiers.setdefault(col_qkey, []).append(True)
                    else:
                        current_feature.qualifiers.setdefault(col_qkey, []).append(col_qvalue)
        mss.parse_seq()
        return mss

    def parse_seq(self) -> None:
        """
        Parse the DDBJ sequence file
        > で始まる行が配列名で、次の > までの行をつないだものを配列にする。// だけの行 (DDBJ の区切り) は読み飛ばす
        """
        if self.seq_file is None:
            raise ValueError("MSS instance must have seq_file")
        seq_id: str | None = None
        seq_lines: list[str] = []
        with open(self.seq_file, encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.rstrip("\n")
                if line.startswith(">"):
                    if seq_id is not None:
                        self.sequences.append(Sequence(id=seq_id, seq="".join(seq_lines)))
                    seq_id, seq_lines = line[1:], []
                elif line.strip() == "//":
                    continue
                elif seq_id is not None:
                    seq_lines.append(line)
        if seq_id is not None:
            self.sequences.append(Sequence(id=seq_id, seq="".join(seq_lines)))

    def to_tsv(self) -> str:
        """Convert the parsed data to TSV format"""
        tsv = []
        for entry in self.entries:
            tsv.extend(entry.to_tsv())
        return "\n".join(["\t".join(row) for row in tsv]) + "\n"

    def to_fasta(self, width: int = 60, separator: bool = False) -> str:
        """Convert the parsed sequence data to FASTA format"""
        return "".join(sequence.to_fasta(width=width, separator=separator) for sequence in self.sequences)

    def write(self, out_ann_file: Path, out_seq_file: Path) -> None:
        with open(out_ann_file, "w", encoding="utf-8") as f:
            f.write(self.to_tsv())
        with open(out_seq_file, "w", encoding="utf-8") as f:
            f.write(self.to_fasta(separator=True))

    def show(self) -> None:
        """[For debugging] Print the parsed data in a human-readable format"""
        print("\nEntries:")
        for entry in self.entries:
            entry.show()

        print("\nSequences:")
        for sequence in self.sequences:
            seq_id = sequence.id
            print(f"\n{seq_id}:")
            print(f"  Length: {len(sequence.seq)}")
            if len(sequence.seq) > SHOW_SEQUENCE_LENGTH:
                print(f"  Sequence: {sequence.seq[:SHOW_SEQUENCE_LENGTH]}...")
            else:
                print(f"  Sequence: {sequence.seq}")

    @staticmethod
    def from_json(json_file: Path) -> "MSS":
        """
        JSONファイルを読み込みMSSインスタンスを作成する
        v1 / v2 の JSON は record v3 に変換してから読む
        """
        record = load_json_to_ddbj_record_instance(json_file, to_record_version="v3")
        return MSS.from_record(record)

    @staticmethod
    def from_record(record: DdbjRecordV3) -> "MSS":
        """
        record v3 から MSS インスタンスを作成する

        entry の source は sequences.common_source を organism / mol_type / qualifier の名前ごとに上書きして展開する。
        locus_tag は projects[].locus_tag_prefix[].prefix と feature の locus_tag_id から組み立てる。
        value の無い qualifier は値の無い行 (True) にする。
        """
        mss = MSS(ann_file=None, seq_file=None)
        mss.entries.append(_common_entry_from_record(record))

        sequences = record.sequences
        common_source = sequences.common_source if sequences is not None else None
        grouped_features = features_by_entry(record)
        locus_tag_prefix = get_locus_tag_prefix(record)

        for entry_data in (sequences.entries if sequences is not None else None) or []:
            entry_id = get_entry_id(entry_data)
            entry = Entry(id=entry_id, name=entry_data.name or entry_id)
            features: list[Feature] = []

            # source feature を作成し、common_source の内容を entry の source で上書きして展開する
            for i, source_feature in enumerate(entry_data.source_features or []):
                qualifiers = merge_source_qualifiers(common_source, source_feature.source)
                if source_feature.definition:
                    qualifiers["ff_definition"] = list(source_feature.definition)
                features.append(
                    Feature(
                        id=source_feature.alias or f"source_{entry_id}_{i + 1}",
                        type="source",
                        location=source_feature.location or "",
                        qualifiers=qualifiers,
                    )
                )

            # topologyがcircularの場合、TOPOLOGY featureを追加
            if entry_data.topology == "circular":
                features.append(Feature(id=f"topology_{entry_id}", type="TOPOLOGY", qualifiers={"circular": [True]}))

            # entry の COMMENT
            for i, comment in enumerate(entry_data.comments or []):
                comment_lines: list[str | bool] = list(comment.split("\n"))
                features.append(
                    Feature(id=f"{entry_id}_comment_{i + 1}", type="COMMENT", qualifiers={"line": comment_lines})
                )

            # その他のfeatureを追加
            for j, feature_data in enumerate(grouped_features.pop(entry_id, [])):
                if not feature_data.type:
                    raise ValueError(f"feature {feature_data.alias!r} of entry {entry_id} has no type")
                features.append(
                    Feature(
                        id=feature_data.alias or f"feature_{entry_id}_{j + 1}",
                        type=feature_data.type,
                        location=feature_data.location or "",
                        qualifiers=feature_qualifier_values(feature_data, locus_tag_prefix),
                    )
                )

            entry.features = features
            # シーケンスの作成 (FASTA のヘッダは ann の1列目と揃えて entry.name を使う)
            mss.sequences.append(Sequence(id=entry.name, seq=entry_data.sequence or ""))
            mss.entries.append(entry)

        for sequence_id, features_left in grouped_features.items():
            warnings.warn(
                f"{len(features_left)} feature(s) refer to unknown entry {sequence_id!r}; dropped",
                UserWarning,
                stacklevel=2,
            )

        return mss


def _common_entry_from_record(record: DdbjRecordV3) -> Entry:
    """record v3 の登録全体の情報から COMMON entry を作る"""
    common_entry = Entry(id="COMMON", name="COMMON")

    def _add(feature_type: str, qualifiers: QualifierValues) -> None:
        if qualifiers:
            common_entry.features.append(
                Feature(
                    id=f"feature_{len(common_entry.features) + 1}",
                    type=feature_type,
                    qualifiers=qualifiers,
                )
            )

    submission = record.submission
    sequences = record.sequences

    _add("DBLINK", _dblink_qualifiers(record.relations))
    _add("SUBMITTER", _submitter_qualifiers(submission.submitters if submission is not None else None))
    for project in record.projects or []:
        for publication in project.publications or []:
            _add("REFERENCE", _reference_qualifiers(publication))
    for comment in (submission.comments if submission is not None else None) or []:
        _add("COMMENT", {"line": list(comment.split("\n"))})
    for structured_comment in (sequences.structured_comments if sequences is not None else None) or []:
        _add("ST_COMMENT", _structured_comment_qualifiers(structured_comment))
    if submission is not None and submission.hold_date:
        _add("DATE", {"hold_date": [hold_date_to_mss(submission.hold_date)]})
    division = _mss_division(record)
    if division is not None:
        _add("DIVISION", {"division": [division]})
    if record.provenance is not None and is_wgs_category(record.provenance.submission_category):
        _add("DATATYPE", {"type": ["WGS"]})
    keywords = sequences.keywords if sequences is not None else None
    if keywords:
        _add("KEYWORD", {"keyword": list(keywords)})

    return common_entry


def is_wgs_category(submission_category: str | None) -> bool:
    """WGS と、MAG-WGS のように末尾が -WGS の登録区分を WGS とする"""
    return submission_category is not None and submission_category.split("-")[-1] == "WGS"


def _mss_division(record: DdbjRecordV3) -> str | None:
    """
    COMMON の DIVISION に書く division。entry の division のうち MSS_DIVISIONS にあるものを返す
    DIVISION は登録全体で 1 つなので、entry の division が 2 つ以上あれば (ENV と BCT など) 警告を出して、
    MSS_DIVISIONS の先頭のものを使う。ann から読み戻すと、全ての entry がその division になる
    """
    entries = (record.sequences.entries if record.sequences is not None else None) or []
    divisions = list(dict.fromkeys(entry.division for entry in entries))
    written = [division for division in divisions if division in MSS_DIVISIONS]
    if not written:
        return None
    if len(divisions) > 1:
        warnings.warn(f"entries have divisions {divisions}; DIVISION is {written[0]!r}", UserWarning, stacklevel=3)
    return written[0]


def _dblink_qualifiers(relations: list[Relation] | None) -> QualifierValues:
    """relations の xref を DBLINK の qualifier にする。DBLINK に書けない db は警告を出して落とす"""
    qualifiers: QualifierValues = {}
    for relation in relations or []:
        if relation.type != "xref" or relation.target is None or not relation.target.id:
            continue
        key = DB_TO_DBLINK_KEY.get(relation.target.db or "")
        if key is None:
            warnings.warn(
                f"xref to {relation.target.db!r} ({relation.target.id}) has no place in DBLINK; dropped",
                UserWarning,
                stacklevel=3,
            )
            continue
        qualifiers.setdefault(key, []).append(relation.target.id)
    return qualifiers


def _submitter_qualifiers(submitters: list[Person] | None) -> QualifierValues:
    """
    submitters を SUBMITTER の qualifier にする
    先頭の人が連絡先で、contact / email / 所属は先頭の人から取る。ab_name は abbreviation を持つ全員から取る
    """
    if not submitters:
        return {}
    values: dict[str, list[str]] = {}
    ab_names = [person.abbreviation for person in submitters if person.abbreviation]
    if ab_names:
        values["ab_name"] = ab_names

    contact = submitters[0]
    if contact.name:
        values["contact"] = [contact.name]
    if contact.email:
        values["email"] = [contact.email]

    organizations = contact.organizations or []
    institutions = [org for org in organizations if org.type != "consortium"]
    consortiums = [org for org in organizations if org.type == "consortium"]
    if institutions:
        institution = institutions[0]
        if institution.url:
            values["url"] = [institution.url]
        if institution.name:
            values["institute"] = [institution.name]
        if institution.department:
            values["department"] = [institution.department]
        address = institution.address
        if address is not None:
            for key, value in (
                ("country", address.country),
                ("state", address.state),
                ("city", address.city),
                ("street", address.street),
                ("zip", address.postal_code),
            ):
                if value:
                    values[key] = [value]
    if consortiums and consortiums[0].name:
        values["consrtm"] = [consortiums[0].name]

    qualifiers: QualifierValues = {}
    for key in SUBMITTER_KEYS:
        if key in values:
            qualifiers[key] = list(values[key])
    return qualifiers


def _reference_qualifiers(publication: Publication) -> QualifierValues:
    """publication を REFERENCE の qualifier にする。year は date の先頭 4 文字"""
    qualifiers: QualifierValues = {}
    if publication.title:
        qualifiers["title"] = [publication.title]
    ab_names = [author.abbreviation for author in publication.authors or [] if author.abbreviation]
    if ab_names:
        qualifiers["ab_name"] = list(ab_names)
    if publication.consortiums:
        qualifiers["consrtm"] = list(publication.consortiums)
    if publication.status:
        qualifiers["status"] = [reference_status_to_mss(publication.status)]
    if publication.date:
        qualifiers["year"] = [publication.date[:4]]
    return qualifiers


def _structured_comment_qualifiers(structured_comment: StructuredComment) -> QualifierValues:
    """structured_comment を ST_COMMENT の qualifier にする。tagset_id を先頭に置く"""
    qualifiers: QualifierValues = {}
    if structured_comment.tagset_id:
        qualifiers["tagset_id"] = [structured_comment.tagset_id]
    for key, value in (structured_comment.fields or {}).items():
        qualifiers[key] = [value]
    return qualifiers
