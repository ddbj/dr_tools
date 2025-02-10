#!/usr/bin/env python3

import json
from pathlib import Path
import re
import logging
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference
from Bio.SeqFeature import Location
from Bio.Data.CodonTable import TranslationError


logger = logging.getLogger(__name__)

def collect_metadata(json_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    JSONデータからメタデータを収集する
    """
    metadata = {}
    
    # COMMONの情報を収集
    common = json_data.get("COMMON", {})
    metadata["common"] = common
    
    # COMMON_SOURCEの情報を収集
    common_source = json_data.get("COMMON_SOURCE", {})
    metadata["common_source"] = common_source
    
    return metadata

def create_submitter(submitter_data: Dict[str, Any]) -> Reference|None:
    """
    jsonのSUBMITTER情報をBioPythonのReferenceオブジェクトに変換する
    (GenBank形式では形式的に登録者情報はReferenceとして扱う)
    入力データの例:
        "SUBMITTER": {
            "ab_name": [
                "Mishima,H.",
                "Tanizawa,Y.",
                "Nakamura,Y."
            ],
            "contact": "Hanako Mishima",
            "email": "mishima@ddbj.nig.ac.jp",
            "url": "https://ddbj.nig.ac.jp",
            "institute": "National Institute of Genetics",
            "department": "DNA Data Bank of Japan",
            "country": "Japan",
            "state": "Shizuoka",
            "city": "Mishima",
            "street": "Yata 1111",
            "zip": "411-8540"
        },
 
    """
    submitters = submitter_data.get("ab_name", [])
    consrtm = submitter_data.get("consrtm", "")
    contact = submitter_data.get("contact", "Contact Person")
    # email = submitter_data.get("email", "email@example.com") email is not visible in GenBank format
    url = submitter_data.get("url")
    institute = submitter_data.get("institute", "Institute Name")
    department = submitter_data.get("department")
    country = submitter_data.get("country", "Country")
    state = submitter_data.get("state")
    city = submitter_data.get("city", "City")
    street = submitter_data.get("street", "Street Name")
    zip = submitter_data.get("zip", "###-####")
    if submitters or consrtm:
        ref = Reference()
        ref.title = "Direct Submission"
        if len(submitters) == 1:
            ref.authors = submitters[0]
        elif len(submitters) > 1:
            ref.authors = ", ".join(submitters[:-1]) + " and " + submitters[-1]
        if consrtm:
            ref.consrtm = consrtm
        contact = "Contact:" + contact
        if department:
            affiliation = f"{institute}, {department};"
        else:
            affiliation = institute
        address = f"{street}, {city}"
        if state:
            address += f", {state}"
        address = f"{address} {zip}, {country}"
        country = submitter_data.get("country")
        if url:
            address += f" URL:{url}"
        ref.journal = f"{contact}, {affiliation}, {address}"
        return ref
    else:
        return None

def create_reference(reference_data: Dict[str, Any]) -> Reference|None:
    """
    ReferenceのjsonデータからRefereeceオブジェクトを作成する
    入力データの例:
            {
            "title": "Complete genome sequencing of Paucilactobacillus hokkaidonensis",
            "ab_name": [
                "Mishima,H.",
                "Tanizawa,Y.",
                "Nakamura,Y."
            ],
            "status": "Unpublished",
            "year": "2023"
        }
    statusがIn pressの場合、journalは"In press, year"となる
    statusがPublishedの場合、journalは"journal volume:pages (year)"となる
    """
    title = reference_data.get("title")
    if title is None:
        return None
    status = reference_data.get("status", "Unpublished")
    ref = Reference()
    ref.title = title
    authors = reference_data.get("ab_name", ["Author,1.", "Author,2.", "Author,3."])
    if len(authors) == 1:
        ref.authors = authors[0]
    elif len(authors) > 1:
        ref.authors = ", ".join(authors[:-1]) + " and " + authors[-1]
    year = reference_data.get("year", "2000")
    if status == "In press":
        journal = reference_data.get("journal", "Journal Name")
        ref.journal = f"{journal} (In press, {year})"
    elif status == "Published":
        volume = reference_data.get("volume", "vol#")
        pages = reference_data.get("start_page", "###")
        end_page = reference_data.get("end_page")
        if end_page:
            pages += "-" + end_page
        ref.journal = f"{journal} {volume}:{pages}, ({year})"
    else:
        ref.journal = f"Unpublished. ({year})"
    refconsrtm = reference_data.get("consrtm")
    if refconsrtm:
        ref.consrtm = refconsrtm
    pubmed = reference_data.get("pubmed")
    if pubmed:
        ref.pubmed = pubmed
    return ref

def create_definition(ff_definition, **kwargs):
    """
    @@で囲まれる部分を正規表現で抽出し、同名の変数で置換します。
    
    :param text: 置換対象のテキスト
    :param kwargs: 置換する変数名とその値
    :return: 置換後のテキスト
    """
    def replacer(match):
        placeholder = match.group(1)  # @@[ ]@@の中の文字列を取得
        return kwargs.get(placeholder, match.group(0))  # 変数があれば置換、なければそのまま

    # 正規表現で@@[ ]@@を検索し、置換処理を行う
    return re.sub(r'@@\[(.*?)\]@@', replacer, ff_definition)

def create_seqfeature(feature_type: str, location: str, qualifiers: Dict[str, List[str]], seq_length: int, topology: str) -> SeqFeature:
    """
    featureの情報からBioPythonのSeqFeatureオブジェクトを作成する
    """
    circular = topology == "circular"  # topoloygyがcircularの場合、true
    feature_location = Location.fromstring(location, length=seq_length, circular=circular)
    
    # qualifiersの変換
    converted_qualifiers = {}
    for key, values in qualifiers.items():
        # Booleanの場合は空文字列に変換
        if len(values) == 0:
            converted_qualifiers[key] = [""]
        elif isinstance(values[0], bool):
            converted_qualifiers[key] = [""]
        else:
            converted_qualifiers[key] = values
    
    return SeqFeature(
        location=feature_location,
        type=feature_type,
        qualifiers=converted_qualifiers
    )

def add_dbxrefs(record: SeqRecord, common_dict: dict):
    """
    DBLINK情報を追加する
    """
    record.dbxrefs = []
    common_dblink = common_dict.get("DBLINK", {})
    bioproject = common_dblink.get("project", "")
    biosample = common_dblink.get("biosample", "")
    sra = common_dblink.get("sequence read archive", "")
    if isinstance(bioproject, list):
        bioproject = ", ".join(bioproject)
    if isinstance(biosample, list):
        biosample = ", ".join(biosample)
    if isinstance(sra, list):
        sra = ", ".join(sra)
    if bioproject:
        record.dbxrefs.append("BioProject:" + bioproject)
    if biosample:
        record.dbxrefs.append("BioSample:" + biosample)
    if sra:
        record.dbxrefs.append("Sequence Read Archive:" + sra)

def add_comment(record: SeqRecord, common_dict: dict):
    """
    ST_COMMENTおよびCOMMENT情報を追加する

    ST_COMMENTは以下の通り GenBankファイルでは


        "ST_COMMENT": {
            "tagset_id": "Genome-Assembly-Data",
            "Assembly Method": "HGAP v. x.x",
            "Genome Coverage": "60x",
            "Sequencing Technology": "Illumina MiSeq; PacBio RSII"
        },    
    GenBankファイルでは以下のような構造になっている
    'structured_comment': defaultdict(<class 'dict'>, {'Genome-Assembly-Data': {'Assembly Method': 'HGAP v. x.x', 'Genome Coverage': '60x', 'Sequencing Technology': 'Illumina MiSeq; PacBio RSII'}}
        
    COMMENTは複数記載できるのでjsonではリストになっている。list[dict[str, list[str]]]

        "COMMENT": [
            {
                "line": ["Comment line1", "line2"]
            },
            {
                "line": ["Annotated by DFAST https://dfast.ddbj.nig.ac.jp/"]
            }
        ],

    
    """
    comments = common_dict.get("COMMENT", [])
    comment_str = "\n".join(["\n".join(comment["line"]) for comment in comments])
    record.annotations["comment"] = comment_str

    st_comment = common_dict.get("ST_COMMENT", {})
    tagset_id = st_comment.get("tagset_id")
    if tagset_id:
        st_comment = {k: v for k, v in st_comment.items() if k != "tagset_id"}
        structured_comment = {tagset_id: st_comment}
        record.annotations["structured_comment"] = structured_comment

def set_date(record: SeqRecord, date: Optional[str] = None):
    """
    日付を設定する
    """
    if date is None:
        date = datetime.now().strftime("%d-%b-%Y").upper()
    record.annotations["date"] = date

def add_translate_qualifier(feature: SeqFeature, seq: Seq) -> None:
    """
    CDS featureにtranslation qualifierを追加する

    transl_exceptの処理を追加
    記載例: (pos:5272379..5272381,aa:Sec)
    現状では Sec: U, Pyl: O にのみ対応している
    aaにTERMが記載されている場合や1つのSeqFeatureに複数のtransl_exceptが記載されている場合、未実装
    """
    def _parse_transl_except(transl_except_value: str) -> Tuple[int, int, str]:
        """
        transl_exceptの開始位置、終了位置、翻訳後のアミノ酸を取得する
        開始位置、終了位置は0-based

        """
        pattern = r'\(pos:(?:complement\()?(\d+)\.\.(\d+)(?:\))?,aa:(\w+)\)'
        match = re.search(pattern, transl_except_value)
        
        if match:
            start = match.group(1)
            end = match.group(2)
            aa = match.group(3)
            return int(start) - 1, int(end), aa
        else:
            raise ValueError(f"Invalid transl_except format: {transl_except_value}")

    def try_translate(feature: SeqFeature, seq: Seq) -> Seq:
        """
        BiopythonのSeqFeature.translate()を試し、codon_startとcodon_tableを考慮してCDSを翻訳する
        デフォルトでは cds=True が指定されており、開始コドンが ATG 以外の場合も M に置換される
        partial locationの場合はエラーが発生するので、エラーが発生した場合は cds=Falseを考慮して翻訳する
        """
        try:
            return feature.translate(seq)
        except TranslationError as e:
            try:
                # start_offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
                # table = int(feature.qualifiers.get("transl_table", [1])[0])
                return feature.translate(seq, cds=False).rstrip("*")
            except TranslationError as e:
                logger.warning("TranslationError: Cannot translate the feature", feature.id)
                logger.warning(feature)
                return None
    dict_aa = {"Sec": "U", "Pyl": "O"} # {"TERM": "*"}  # TERM not implemented

    # transl_except が記載されている場合、その箇所を @@@ に置換して CDS を切り出し、
    # CDS とアミノ酸配列中での位置を確認しておく
    # translate 実行前に @@@ をNNNに置換しておくと該当箇所は X に翻訳されるので、翻訳後にアミノ酸を置換する
    if "transl_except" in feature.qualifiers:
        transl_except = feature.qualifiers["transl_except"]
        if len(transl_except) > 1:
            raise ValueError("Multiple transl_except qualifiers are not supported")
        transl_except = transl_except[0]
        start, end, aa = _parse_transl_except(transl_except)
        assert aa in dict_aa
        seq = seq[:start] + Seq("@@@") + seq[end:]  # replace the amino acid with @@@
        cds = feature.location.extract(seq)
        start_offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1

         # CDSとアミノ酸配列中でのtranl_exceptの開始位置 (0-based)
        except_index_nuc = cds.index("@@@")
        except_index_aa = (except_index_nuc - start_offset) // 3

        # print(except_index_nuc, except_index_aa)
        # print(feature.qualifiers)
        # print(cds[except_index_nuc:except_index_nuc+3])
        seq = seq.replace("@@@", "NNN")

    translation = try_translate(feature, seq)

    if "transl_except" in feature.qualifiers:
        assert translation[except_index_aa] == "X"
        translation = translation[:except_index_aa] + dict_aa[aa] + translation[except_index_aa+1:]

    if translation:
        feature.qualifiers["translation"] = [str(translation)]

def json_to_seqrecords(json_file: Path) -> List[SeqRecord]:
    """
    JSONデータをBioPythonのSeqRecordオブジェクトのリストに変換する
    """
    with open(json_file) as f:
        json_data = json.load(f)
    records = []
    
    common = json_data.get("COMMON", {})

    # COMMON_SOURCE, COMMON_METAの内容を取得
    common_source = json_data.get("COMMON_SOURCE", {})
    common_meta = json_data.get("COMMON_META", {})

    # エントリーごとにSeqRecordを作成
    for entry_data in json_data["ENTRIES"]:
        # シーケンスの作成
        seq = Seq(entry_data["sequence"])
        
        # SeqRecordの作成
        record = SeqRecord(
            seq=seq,
            id=entry_data["id"],
            name=entry_data["name"],
            description="test"
        )
        
        # molecule_typeの設定
        mol_type_value = common_source.get("mol_type", "DNA")
        if "DNA" in mol_type_value:
            record.annotations["molecule_type"] = "DNA"
        elif "RNA" in mol_type_value:
            record.annotations["molecule_type"] = "RNA"
        else:
            raise ValueError(f"Invalid molecule_type: {mol_type_value}. Must contain 'DNA' or 'RNA'.")
        
        # Division
        division = common_meta.get("division", "UNK")
        record.annotations["data_file_division"] = division

        # topologyの設定
        topology = entry_data.get("topology", "linear")
        record.annotations["topology"] = topology
 
        # organismの設定
        organism = common_source.get("organism", "Unknown organism")
        record.annotations["organism"] = organism
        record.annotations["source"] = organism

        # keywordsの設定
        keywords = common.get("KEYWORD", {}).get("keyword", [])
        if keywords:
            record.annotations["keywords"] = keywords

        # dbxrefsの設定
        add_dbxrefs(record, common)

        # commentの設定
        add_comment(record, common)
        
        # dateの設定
        set_date(record, common_meta.get("date"))

        # submitter, referencesの設定
        references = []
        submitter = create_submitter(common.get("SUBMITTER", {}))
        if submitter:
            references.append(submitter)
        for ref_data in common.get("REFERENCE", []):
            ref = create_reference(ref_data)
            if ref:
                references.append(ref)
        record.annotations["references"] = references

        seq_length = len(seq)

        features = []
        
        # source featureの作成
        for feature_data in entry_data["features"]:
            if feature_data["type"] == "source":
                # COMMON_SOURCEの内容とマージ
                qualifiers = {k: [v] for k, v in common_source.items()}
                qualifiers.update(feature_data["qualifiers"])
                feature = create_seqfeature(
                    "source",
                    feature_data["location"],
                    qualifiers,
                    seq_length, topology
                )  # seq_lengthとtopologyはlocationを解析するために必要
                features.append(feature)
                ff_definition = qualifiers.get("ff_definition", ["@@[organism]@@"])[0]
                source_qualifiers = {k: v[0] for k, v in qualifiers.items() if k != "ff_definition" and len(v) == 1}
                source_qualifiers["submitter_seqid"] = entry_data["name"]
                definition = create_definition(ff_definition, **source_qualifiers)
                record.description = definition
            else:
                # source feature以外のfeatureを作成
                feature = create_seqfeature(
                    feature_data["type"],
                    feature_data["location"],
                    feature_data["qualifiers"],
                    seq_length, topology
                )

                if feature.type == "CDS":
                    # print(feature)
                    # print(seq[:10])
                    # print(len(seq))
                    add_translate_qualifier(feature, seq)

                features.append(feature)

        
        # featuresをrecordに設定
        record.features = features
        records.append(record)
    
    return records



def json2gbk_main():
    """
    JSONファイルをGenBankファイルに変換する
    出力 (-o, --output_file) が指定されていない場合は、標準出力に出力する
    """

    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Convert JSON file to GenBank file')
    parser.add_argument('input_file', type=str, help='Input JSON file')
    parser.add_argument('-o', '--output_file', type=str, help='Output GenBank file', default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file

    # JSONファイルの読み込み
    # with open(input_file) as f:
    # # with open("draft_genome.json") as f:
    #     json_data = json.load(f)
    
    # JSONデータをSeqRecordに変換
    records = json_to_seqrecords(input_file)
    
    # GenBank形式で出力
    if output_file:
        with open(output_file, "w") as f:
            SeqIO.write(records, f, "genbank")
    else:
        for r in records:
            print(r.format("genbank"))

# if __name__ == "__main__":
#     json2gbk_main() 