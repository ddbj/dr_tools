import argparse
import logging
import re
import sys
from collections.abc import Mapping, Sequence
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqFeature import Location, Reference, SeqFeature
from Bio.SeqRecord import SeqRecord
from ddbj_record.schema.v3 import Organization, Person, Publication, Relation, StructuredComment

from dr_tools.json_utils import (
    feature_qualifier_values,
    features_by_entry,
    get_entry_id,
    get_locus_tag_prefix,
    load_json_to_ddbj_record_instance,
    merge_source_qualifiers,
)
from dr_tools.MSS import DB_TO_DBLINK_KEY, reference_status_to_mss
from dr_tools.translation import translate_cds_with_transl_except

logger = logging.getLogger(__name__)

# GenBank の DBLINK に書く順と見出し
_DBXREF_LABELS = {
    "project": "BioProject",
    "biosample": "BioSample",
    "sequence read archive": "Sequence Read Archive",
}


def _join_authors(authors: list[str]) -> str:
    """著者の略記を GenBank の AUTHORS の書き方 ("A, B and C") にする"""
    if len(authors) == 1:
        return authors[0]
    if len(authors) > 1:
        return ", ".join(authors[:-1]) + " and " + authors[-1]
    return ""


def _institution(organizations: Sequence[Organization] | None) -> Organization | None:
    """所属のうち consortium でない最初のものを返す"""
    for organization in organizations or []:
        if organization.type != "consortium":
            return organization
    return None


def _consortium_name(organizations: Sequence[Organization] | None) -> str | None:
    """所属のうち consortium の最初のものの名前を返す"""
    for organization in organizations or []:
        if organization.type == "consortium" and organization.name:
            return organization.name
    return None


def create_submitter(submitters: Sequence[Person] | None) -> Reference | None:
    """
    record の submitters を BioPython の Reference オブジェクトに変換する
    (GenBank形式では形式的に登録者情報はReferenceとして扱う)
    先頭の人が連絡先で、連絡先の名前と所属 (institution の名前、部署、住所、URL) を JOURNAL に書く。
    abbreviation を持つ人が 1 人もおらず、consortium も無ければ None を返す
    """
    if not submitters:
        return None
    ab_names = [person.abbreviation for person in submitters if person.abbreviation]
    contact = submitters[0]
    consrtm = _consortium_name(contact.organizations)
    if not ab_names and not consrtm:
        return None

    institution = _institution(contact.organizations)
    address = institution.address if institution is not None else None

    ref = Reference()  # type: ignore[no-untyped-call]
    ref.title = "Direct Submission"
    ref.authors = _join_authors(ab_names)
    if consrtm:
        ref.consrtm = consrtm
    contact_text = "Contact:" + (contact.name or "Contact Person")
    institute = (institution.name if institution is not None else None) or "Institute Name"
    department = institution.department if institution is not None else None
    affiliation = f"{institute}, {department};" if department else institute
    street = (address.street if address is not None else None) or "Street Name"
    city = (address.city if address is not None else None) or "City"
    state = address.state if address is not None else None
    postal_code = (address.postal_code if address is not None else None) or "###-####"
    country = (address.country if address is not None else None) or "Country"
    address_text = f"{street}, {city}"
    if state:
        address_text += f", {state}"
    address_text = f"{address_text} {postal_code}, {country}"
    url = institution.url if institution is not None else None
    if url:
        address_text += f" URL:{url}"
    ref.journal = f"{contact_text}, {affiliation}, {address_text}"
    return ref


def create_reference(publication: Publication) -> Reference | None:
    """
    record の publication から Reference オブジェクトを作成する
    title が無ければ None を返す
    statusがIn pressの場合、journalは"journal (In press, year)"となる
    statusがPublishedの場合、journalは"journal volume:pages, (year)"となる
    year は date の先頭 4 文字
    """
    if publication.title is None:
        return None
    status = reference_status_to_mss(publication.status or "unpublished")
    ref = Reference()  # type: ignore[no-untyped-call]
    ref.title = publication.title
    authors = [author.abbreviation for author in publication.authors or [] if author.abbreviation]
    if not authors:
        authors = ["Author,1.", "Author,2.", "Author,3."]
    ref.authors = _join_authors(authors)
    year = (publication.date or "2000")[:4]
    journal = publication.journal or "Journal Name"
    if status == "In press":
        ref.journal = f"{journal} (In press, {year})"
    elif status == "Published":
        volume = publication.volume or "vol#"
        pages = publication.pages_from or "###"
        if publication.pages_to:
            pages += "-" + publication.pages_to
        ref.journal = f"{journal} {volume}:{pages}, ({year})"
    else:
        ref.journal = f"Unpublished. ({year})"
    if publication.consortiums:
        ref.consrtm = ", ".join(publication.consortiums)
    if publication.pubmed_id:
        ref.pubmed_id = publication.pubmed_id
    return ref


def create_definition(ff_definition: str, **kwargs: str) -> str:
    """
    @@で囲まれる部分を正規表現で抽出し、同名の変数で置換します。

    :param text: 置換対象のテキスト
    :param kwargs: 置換する変数名とその値
    :return: 置換後のテキスト
    """

    def replacer(match: re.Match[str]) -> str:
        placeholder = match.group(1)  # @@[ ]@@の中の文字列を取得
        return kwargs.get(placeholder, match.group(0))  # 変数があれば置換、なければそのまま

    # 正規表現で@@[ ]@@を検索し、置換処理を行う
    return re.sub(r"@@\[(.*?)\]@@", replacer, ff_definition)


def create_seqfeature(
    feature_type: str, location: str, qualifiers: Mapping[str, Sequence[str | bool]], seq_length: int, topology: str
) -> SeqFeature:
    """
    featureの情報からBioPythonのSeqFeatureオブジェクトを作成する
    値を持たない qualifier (True) は空文字列にする
    """
    circular = topology == "circular"  # topoloygyがcircularの場合、true
    feature_location = Location.fromstring(location, length=seq_length, circular=circular)  # type: ignore[no-untyped-call]

    # qualifiersの変換
    converted_qualifiers: dict[str, list[str]] = {}
    for key, values in qualifiers.items():
        # Booleanの場合は空文字列に変換
        if len(values) == 0 or isinstance(values[0], bool):
            converted_qualifiers[key] = [""]
        else:
            converted_qualifiers[key] = [str(value) for value in values]

    return SeqFeature(  # type: ignore[no-untyped-call]
        location=feature_location, type=feature_type, qualifiers=converted_qualifiers
    )


def add_dbxrefs(record: SeqRecord, relations: Sequence[Relation] | None) -> None:
    """
    relations の xref (bioproject / biosample / insdc.sra) を DBLINK として追加する
    同じ DB の id が複数あれば ", " でつなぐ
    """
    record.dbxrefs = []
    ids_by_key: dict[str, list[str]] = {}
    for relation in relations or []:
        if relation.type != "xref" or relation.target is None or not relation.target.id:
            continue
        key = DB_TO_DBLINK_KEY.get(relation.target.db or "")
        if key is None:
            continue
        ids_by_key.setdefault(key, []).append(relation.target.id)
    for key, label in _DBXREF_LABELS.items():
        if key in ids_by_key:
            record.dbxrefs.append(f"{label}:" + ", ".join(ids_by_key[key]))


def add_comment(
    record: SeqRecord, comments: Sequence[str] | None, structured_comments: Sequence[StructuredComment] | None
) -> None:
    """
    COMMENT および ST_COMMENT の情報を追加する

    comments は submission.comments で、1 つのコメントは行を改行でつないだ文字列。全てを改行でつないで COMMENT にする
    structured_comments は sequences.structured_comments で、GenBank ファイルでは以下のような構造になる
        'structured_comment': {'Genome-Assembly-Data': {'Assembly Method': 'HGAP v. x.x', 'Genome Coverage': '60x', 'Sequencing Technology': 'Illumina MiSeq; PacBio RSII'}}
    """
    record.annotations["comment"] = "\n".join(comments or [])

    structured_comment: dict[str, dict[str, str]] = {}
    for item in structured_comments or []:
        if item.tagset_id:
            structured_comment[item.tagset_id] = dict(item.fields or {})
    if structured_comment:
        record.annotations["structured_comment"] = structured_comment  # type: ignore[assignment]


def set_date(record: SeqRecord, date: str | None = None) -> None:
    """
    日付を設定する
    """
    if date is None:
        date = datetime.now().astimezone().strftime("%d-%b-%Y").upper()
    record.annotations["date"] = date


def add_translate_qualifier(feature: SeqFeature, seq: Seq) -> None:
    """
    CDS featureにtranslation qualifierを追加する
    翻訳と transl_except の解釈は dr_tools.translation に一本化してある
    翻訳できない feature (塩基以外の文字を含む配列、無い transl_table、数字でない codon_start) では
    translation を付けず、警告だけ出す
    """
    try:
        translation = translate_cds_with_transl_except(feature, seq)
    except (TranslationError, KeyError, ValueError) as err:
        logger.warning("Cannot translate the feature %s: %s", feature.id, err)
        return
    if translation:
        feature.qualifiers["translation"] = [str(translation)]


def json_to_seqrecords(json_file: Path) -> list[SeqRecord]:
    """
    JSONデータをBioPythonのSeqRecordオブジェクトのリストに変換する
    JSON は record v1 / v2 / v3 のどれでもよい (v1 / v2 は v3 に変換してから読む)
    """
    record = load_json_to_ddbj_record_instance(json_file, to_record_version="v3")

    records = []

    submission = record.submission
    sequences = record.sequences
    common_source = sequences.common_source if sequences is not None else None
    locus_tag_prefix = get_locus_tag_prefix(record)
    grouped_features = features_by_entry(record)

    # submitter, referencesの設定 (全ての entry で同じ)
    references = []
    submitter = create_submitter(submission.submitters if submission is not None else None)
    if submitter:
        references.append(submitter)
    for project in record.projects or []:
        for publication in project.publications or []:
            ref = create_reference(publication)
            if ref:
                references.append(ref)

    # エントリーごとにSeqRecordを作成
    for entry_data in (sequences.entries if sequences is not None else None) or []:
        entry_id = get_entry_id(entry_data)
        entry_name = entry_data.name or entry_id
        seq = Seq(entry_data.sequence or "")

        seq_record = SeqRecord(seq=seq, id=entry_id, name=entry_name, description="")

        # source featureの作成 (COMMON_SOURCEの内容を entry の source で上書きして展開する)
        source_features = entry_data.source_features or []
        first_source = (
            merge_source_qualifiers(common_source, source_features[0].source)
            if source_features
            else merge_source_qualifiers(common_source, None)
        )

        # molecule_typeの設定
        mol_type_values = first_source.get("mol_type", ["DNA"])
        mol_type_value = str(mol_type_values[0]) if mol_type_values else "DNA"
        if "DNA" in mol_type_value:
            seq_record.annotations["molecule_type"] = "DNA"
        elif "RNA" in mol_type_value:
            seq_record.annotations["molecule_type"] = "RNA"
        else:
            raise ValueError(f"Invalid molecule_type: {mol_type_value}. Must contain 'DNA' or 'RNA'.")

        # Division
        seq_record.annotations["data_file_division"] = entry_data.division or "UNK"

        # topologyの設定
        topology = entry_data.topology or "linear"
        seq_record.annotations["topology"] = topology

        # organismの設定
        organism_values = first_source.get("organism", [])
        organism = str(organism_values[0]) if organism_values else "Unknown organism"
        seq_record.annotations["organism"] = organism
        seq_record.annotations["source"] = organism

        # keywordsの設定
        keywords = sequences.keywords if sequences is not None else None
        if keywords:
            seq_record.annotations["keywords"] = list(keywords)  # type: ignore[assignment]

        # dbxrefsの設定
        add_dbxrefs(seq_record, record.relations)

        # commentの設定
        add_comment(
            seq_record,
            submission.comments if submission is not None else None,
            sequences.structured_comments if sequences is not None else None,
        )

        # dateの設定
        set_date(seq_record, None)

        seq_record.annotations["references"] = list(references)  # type: ignore[assignment]

        seq_length = len(seq)
        features = []

        for i, source_feature in enumerate(source_features):
            qualifiers = merge_source_qualifiers(common_source, source_feature.source)
            feature = create_seqfeature(
                "source", source_feature.location or "", qualifiers, seq_length, topology
            )  # seq_lengthとtopologyはlocationを解析するために必要
            features.append(feature)
            if i == 0:
                # DEFINITION は最初の source feature の ff_definition から作る
                # @@[organism]@@ などは source の qualifier の値で、@@[submitter_seqid]@@ と @@[entry]@@ は entry の名前で展開する
                ff_definition = source_feature.definition[0] if source_feature.definition else "@@[organism]@@"
                source_qualifiers = {
                    key: str(values[0])
                    for key, values in qualifiers.items()
                    if len(values) == 1 and not isinstance(values[0], bool)
                }
                source_qualifiers["submitter_seqid"] = entry_name
                source_qualifiers["entry"] = entry_name
                seq_record.description = create_definition(ff_definition, **source_qualifiers)

        for feature_data in grouped_features.get(entry_id, []):
            # source feature以外のfeatureを作成
            feature = create_seqfeature(
                feature_data.type or "",
                feature_data.location or "",
                feature_qualifier_values(feature_data, locus_tag_prefix),
                seq_length,
                topology,
            )

            if feature.type == "CDS":
                add_translate_qualifier(feature, seq)

            features.append(feature)

        # featuresをrecordに設定
        seq_record.features = features
        records.append(seq_record)

    return records


def json2gbk_main() -> None:
    """
    JSONファイルをGenBankファイルに変換する
    出力 (-o, --output_file) が指定されていない場合は、標準出力に出力する
    """
    parser = argparse.ArgumentParser(description="Convert JSON file to GenBank file")
    parser.add_argument("input_file", type=str, help="Input JSON file")
    parser.add_argument("-o", "--output_file", type=str, help="Output GenBank file", default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file

    # JSONデータをSeqRecordに変換
    records = json_to_seqrecords(input_file)

    # GenBank形式で出力
    if output_file:
        with open(output_file, "w", encoding="utf-8") as f:
            SeqIO.write(records, f, "genbank")
    else:
        for r in records:
            print(r.format("genbank"))
