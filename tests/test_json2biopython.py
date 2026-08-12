import io
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, cast

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import LocationParserError, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from dr_tools.json2biopython import (add_comment, add_dbxrefs,
                                     add_translate_qualifier,
                                     create_definition, create_reference,
                                     create_seqfeature, create_submitter,
                                     json_to_seqrecords, set_date)

# --- JSON builders (直接 ddbj_record v1 スキーマを満たす最小 JSON を組み立てる) ---


def _build_source_feature(
    *,
    location: str = "1..12",
    ff_definition: str = "@@[organism]@@",
    extra_qualifiers: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, Any]:
    qualifiers: Dict[str, List[str]] = {"ff_definition": [ff_definition]}
    if extra_qualifiers:
        qualifiers.update(extra_qualifiers)
    return {"id": "feature_source", "type": "source", "location": location, "qualifiers": qualifiers}


def _build_cds_feature(
    location: str,
    *,
    product: str = "hypothetical protein",
    transl_table: str = "11",
    codon_start: str = "1",
) -> Dict[str, Any]:
    return {
        "id": "feature_cds",
        "type": "CDS",
        "location": location,
        "qualifiers": {
            "product": [product],
            "transl_table": [transl_table],
            "codon_start": [codon_start],
        },
    }


def _build_entry(
    *,
    entry_id: str = "entry1",
    name: str = "entry1",
    entry_type: str = "chromosome",
    topology: str = "linear",
    sequence: str = "atgaaacgctaa",
    features: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    if features is None:
        features = [_build_source_feature(location=f"1..{len(sequence)}")]
    return {
        "id": entry_id,
        "name": name,
        "type": entry_type,
        "topology": topology,
        "sequence": sequence,
        "features": features,
    }


def _build_record_json(
    *,
    common: Optional[Dict[str, Any]] = None,
    common_source: Optional[Dict[str, Any]] = None,
    common_meta: Optional[Dict[str, Any]] = None,
    entries: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    base_common: Dict[str, Any] = {
        "SUBMITTER": {
            "ab_name": ["Smith,J."],
            "contact": "John Smith",
            "email": "smith@example.com",
            "institute": "Example Institute",
            "country": "Japan",
            "city": "Tokyo",
            "street": "1-1-1 Test",
            "zip": "100-0001",
        },
        "trad_submission_category": "GNM",
    }
    if common:
        base_common.update(common)

    base_common_source: Dict[str, Any] = {
        "organism": "Escherichia coli",
        "mol_type": "genomic DNA",
    }
    if common_source:
        base_common_source.update(common_source)

    base_common_meta: Dict[str, Any] = {"division": "BCT"}
    if common_meta:
        base_common_meta.update(common_meta)

    if entries is None:
        entries = [_build_entry()]

    return {
        "schema_version": "v1.0",
        "COMMON": base_common,
        "COMMON_SOURCE": base_common_source,
        "COMMON_META": base_common_meta,
        "ENTRIES": entries,
    }


def _write_record_json(tmp_path: Path, data: Dict[str, Any]) -> Path:
    path = tmp_path / "record.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


# --- add_translate_qualifier ---


def test_add_translate_qualifier_valid_coding_sequence_sets_single_element_translation_list() -> None:
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 12), type="CDS", qualifiers={"transl_table": ["11"], "codon_start": ["1"]}  # type: ignore[no-untyped-call]
    )
    add_translate_qualifier(feature, Seq("atgaaacgctaa"))
    assert feature.qualifiers["translation"] == ["MKR"]


def test_add_translate_qualifier_non_nucleotide_characters_no_exception_and_no_translation_qualifier() -> None:
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 9), type="CDS", qualifiers={"transl_table": ["11"], "codon_start": ["1"]}  # type: ignore[no-untyped-call]
    )
    # 例外を投げないことそのものが検証対象なので、try/except で囲まない
    add_translate_qualifier(feature, Seq("atgXYZaaa"))
    assert "translation" not in feature.qualifiers


def test_add_translate_qualifier_translation_reduces_to_empty_string_does_not_set_qualifier() -> None:
    # 単一のストップコドンのみの CDS はストップ記号を除くと空文字列になる
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 3), type="CDS", qualifiers={"transl_table": ["11"], "codon_start": ["1"]}  # type: ignore[no-untyped-call]
    )
    add_translate_qualifier(feature, Seq("taa"))
    assert "translation" not in feature.qualifiers


# --- json_to_seqrecords: id / name ---


def test_json_to_seqrecords_id_is_entry_id_and_name_is_entry_name(tmp_path: Path) -> None:
    entries = [_build_entry(entry_id="chr1", name="mychromosome", features=[])]
    data = _build_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].id == "chr1"
    assert records[0].name == "mychromosome"


# --- json_to_seqrecords: description / ff_definition ---


def test_json_to_seqrecords_source_feature_placeholders_replaced_with_source_qualifier_values(tmp_path: Path) -> None:
    features = [_build_source_feature(ff_definition="@@[organism]@@ DNA, complete genome")]
    entries = [_build_entry(features=features)]
    data = _build_record_json(common_source={"organism": "Vibrio cholerae"}, entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].description == "Vibrio cholerae DNA, complete genome"


def test_json_to_seqrecords_no_source_feature_description_is_empty_string(tmp_path: Path) -> None:
    entries = [_build_entry(features=[])]
    data = _build_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].description == ""


def test_json_to_seqrecords_submitter_seqid_placeholder_uses_entry_name_even_if_qualifier_overrides_it(
    tmp_path: Path,
) -> None:
    features = [
        _build_source_feature(
            ff_definition="@@[submitter_seqid]@@",
            extra_qualifiers={"submitter_seqid": ["should_be_ignored"]},
        )
    ]
    entries = [_build_entry(name="realname", features=features)]
    data = _build_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].description == "realname"


def test_json_to_seqrecords_source_qualifier_with_multiple_values_placeholder_left_unreplaced(tmp_path: Path) -> None:
    features = [
        _build_source_feature(
            ff_definition="@@[note]@@ region",
            extra_qualifiers={"note": ["a", "b"]},
        )
    ]
    entries = [_build_entry(features=features)]
    data = _build_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].description == "@@[note]@@ region"


# --- json_to_seqrecords: mol_type ---


def test_json_to_seqrecords_mol_type_containing_dna_sets_molecule_type_dna(tmp_path: Path) -> None:
    data = _build_record_json(common_source={"mol_type": "genomic DNA"})
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["molecule_type"] == "DNA"


def test_json_to_seqrecords_mol_type_containing_rna_sets_molecule_type_rna(tmp_path: Path) -> None:
    data = _build_record_json(common_source={"mol_type": "genomic RNA"})
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["molecule_type"] == "RNA"


def test_json_to_seqrecords_mol_type_without_dna_or_rna_raises_value_error(tmp_path: Path) -> None:
    data = _build_record_json(common_source={"mol_type": "protein"})
    with pytest.raises(ValueError, match="protein"):
        json_to_seqrecords(_write_record_json(tmp_path, data))


# --- json_to_seqrecords: annotations ---


def test_json_to_seqrecords_multiple_entries_each_get_independent_id_and_topology(tmp_path: Path) -> None:
    entries = [
        _build_entry(entry_id="chr", name="chr", topology="circular", features=[]),
        _build_entry(entry_id="plasmid1", name="plasmid1", entry_type="plasmid", topology="linear", features=[]),
    ]
    data = _build_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert [r.id for r in records] == ["chr", "plasmid1"]
    assert records[0].annotations["topology"] == "circular"
    assert records[1].annotations["topology"] == "linear"


def test_json_to_seqrecords_annotations_include_division_and_organism(tmp_path: Path) -> None:
    data = _build_record_json(
        common_source={"organism": "Salmonella enterica"},
        common_meta={"division": "BCT"},
    )
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["data_file_division"] == "BCT"
    assert records[0].annotations["organism"] == "Salmonella enterica"
    assert records[0].annotations["source"] == "Salmonella enterica"


def test_json_to_seqrecords_keywords_present_sets_keywords_annotation(tmp_path: Path) -> None:
    data = _build_record_json(common={"KEYWORD": {"keyword": ["WGS", "STANDARD_DRAFT"]}})
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["keywords"] == ["WGS", "STANDARD_DRAFT"]  # type: ignore[comparison-overlap]


def test_json_to_seqrecords_keywords_absent_does_not_set_keywords_annotation(tmp_path: Path) -> None:
    data = _build_record_json()
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert "keywords" not in records[0].annotations


# --- json_to_seqrecords: CDS translation wiring (実データ, full_record) ---


def test_json_to_seqrecords_cds_feature_gets_translation_other_feature_types_do_not(
    full_record: Dict[str, Path],
) -> None:
    records = json_to_seqrecords(full_record["v1"])
    rec = records[0]
    seq_str = str(rec.seq)

    cds_features = [f for f in rec.features if f.type == "CDS"]
    assert cds_features, "この entry には CDS が含まれているはず"
    first_cds = cds_features[0]
    assert "translation" in first_cds.qualifiers
    assert len(first_cds.qualifiers["translation"]) == 1

    # dr_tools の翻訳関数を経由せず、Bio.Seq.translate を直接使った独立したオラクルで検算する
    start = int(first_cds.location.start)  # type: ignore[union-attr]
    end = int(first_cds.location.end)  # type: ignore[union-attr]
    table = int(first_cds.qualifiers["transl_table"][0])
    oracle = str(Seq(seq_str[start:end]).translate(table=table, cds=True))  # type: ignore[no-untyped-call]
    assert first_cds.qualifiers["translation"][0] == oracle

    non_cds_types = {f.type for f in rec.features if f.type != "CDS"}
    assert non_cds_types, "この entry には CDS 以外の feature も含まれているはず"
    for feature in rec.features:
        if feature.type != "CDS":
            assert "translation" not in feature.qualifiers


# --- json_to_seqrecords: 実データでの sanity check ---


def test_json_to_seqrecords_real_complete_genome_v1_entries_and_description(eg_complete: Dict[str, Path]) -> None:
    records = json_to_seqrecords(eg_complete["v1_json"])
    assert [r.id for r in records] == ["chromosome", "pPLH-1", "pPLH-2"]
    assert records[0].description == "Paucilactobacillus hokkaidonensis LOOC260 DNA, complete genome"


def test_json_to_seqrecords_real_vrl_v1_description_and_rna_molecule_type(eg_vrl: Dict[str, Path]) -> None:
    records = json_to_seqrecords(eg_vrl["v1_json"])
    assert len(records) == 1
    assert records[0].description == (
        "Severe acute respiratory syndrome coronavirus 2 hCov-19/Japan/SZ-NIG-12345/2021 RNA, complete genome"
    )
    assert records[0].annotations["molecule_type"] == "RNA"


def test_json_to_seqrecords_writes_to_genbank_format_without_error(full_record: Dict[str, Path]) -> None:
    records = json_to_seqrecords(full_record["v1"])
    buf = io.StringIO()
    SeqIO.write(records, buf, "genbank")
    output = buf.getvalue()
    assert output.startswith("LOCUS")
    assert "ORIGIN" in output


# --- create_submitter ---


def test_create_submitter_no_ab_name_no_consrtm_returns_none() -> None:
    assert create_submitter({"contact": "X", "institute": "Y"}) is None


def test_create_submitter_single_author_sets_authors_verbatim() -> None:
    ref = create_submitter({"ab_name": ["Smith,J."]})
    assert ref is not None
    assert ref.authors == "Smith,J."


def test_create_submitter_multiple_authors_joins_with_and_before_last() -> None:
    ref = create_submitter({"ab_name": ["Smith,J.", "Doe,A.", "Lee,K."]})
    assert ref is not None
    assert ref.authors == "Smith,J., Doe,A. and Lee,K."


def test_create_submitter_consrtm_only_without_ab_name_returns_reference() -> None:
    ref = create_submitter({"ab_name": [], "consrtm": "Some Consortium"})
    assert ref is not None
    assert ref.consrtm == "Some Consortium"
    assert ref.authors == ""


# --- create_reference ---


def test_create_reference_no_title_returns_none() -> None:
    assert create_reference({}) is None


def test_create_reference_status_published_formats_journal_with_volume_and_page_range() -> None:
    ref = create_reference({
        "title": "T",
        "status": "Published",
        "journal": "J Test",
        "volume": "10",
        "start_page": "1",
        "end_page": "5",
        "year": "2020",
        "ab_name": ["Smith,J."],
    })
    assert ref is not None
    assert ref.journal == "J Test 10:1-5, (2020)"


def test_create_reference_status_published_without_end_page_uses_start_page_only() -> None:
    ref = create_reference({
        "title": "T",
        "status": "Published",
        "journal": "J Test",
        "volume": "10",
        "start_page": "42",
        "year": "2020",
        "ab_name": ["Smith,J."],
    })
    assert ref is not None
    assert ref.journal == "J Test 10:42, (2020)"


def test_create_reference_status_in_press_formats_journal_with_in_press() -> None:
    ref = create_reference({
        "title": "T",
        "status": "In press",
        "journal": "J Test",
        "year": "2021",
        "ab_name": ["Smith,J."],
    })
    assert ref is not None
    assert ref.journal == "J Test (In press, 2021)"


def test_create_reference_status_unpublished_default_formats_journal() -> None:
    ref = create_reference({"title": "T", "year": "2022", "ab_name": ["Smith,J."]})
    assert ref is not None
    assert ref.journal == "Unpublished. (2022)"


def test_create_reference_pubmed_sets_pubmed_id_attribute() -> None:
    ref = create_reference({"title": "T", "pubmed": "30000000", "ab_name": ["Smith,J."]})
    assert ref is not None
    assert ref.pubmed_id == "30000000"


def test_create_reference_no_pubmed_pubmed_id_is_empty_string() -> None:
    ref = create_reference({"title": "T", "ab_name": ["Smith,J."]})
    assert ref is not None
    assert ref.pubmed_id == ""


def test_create_reference_single_author_authors_set_verbatim() -> None:
    ref = create_reference({"title": "T", "ab_name": ["Smith,J."]})
    assert ref is not None
    assert ref.authors == "Smith,J."


def test_create_reference_multiple_authors_joined_with_and_before_last() -> None:
    ref = create_reference({"title": "T", "ab_name": ["Smith,J.", "Doe,A.", "Lee,K."]})
    assert ref is not None
    assert ref.authors == "Smith,J., Doe,A. and Lee,K."


# --- create_definition ---


def test_create_definition_known_placeholder_replaced() -> None:
    assert create_definition("@@[organism]@@ DNA", organism="Escherichia coli") == "Escherichia coli DNA"


def test_create_definition_unknown_placeholder_left_unchanged() -> None:
    assert create_definition("@@[foo]@@ region", organism="Escherichia coli") == "@@[foo]@@ region"


def test_create_definition_multiple_placeholders_all_replaced() -> None:
    result = create_definition(
        "@@[organism]@@ @@[strain]@@, @@[submitter_seqid]@@",
        organism="E. coli",
        strain="K12",
        submitter_seqid="chr1",
    )
    assert result == "E. coli K12, chr1"


def test_create_definition_no_placeholders_returns_text_unchanged() -> None:
    text = "Escherichia coli DNA, complete genome"
    assert create_definition(text) == text


# --- create_seqfeature ---


def test_create_seqfeature_boolean_qualifier_value_converted_to_empty_string() -> None:
    qualifiers = cast(Dict[str, List[str]], {"pseudo": [True]})
    feature = create_seqfeature("gene", "1..5", qualifiers, 5, "linear")
    assert feature.qualifiers["pseudo"] == [""]


def test_create_seqfeature_empty_qualifier_list_converted_to_empty_string() -> None:
    feature = create_seqfeature("gene", "1..5", {"note": []}, 5, "linear")
    assert feature.qualifiers["note"] == [""]


def test_create_seqfeature_normal_string_qualifier_preserved() -> None:
    feature = create_seqfeature("CDS", "1..5", {"product": ["hypothetical protein"]}, 5, "linear")
    assert feature.qualifiers["product"] == ["hypothetical protein"]
    assert feature.type == "CDS"
    assert int(feature.location.start) == 0  # type: ignore[union-attr]
    assert int(feature.location.end) == 5  # type: ignore[union-attr]


def test_create_seqfeature_circular_topology_allows_origin_wrapping_location() -> None:
    feature = create_seqfeature("misc_feature", "9..2", {}, 10, "circular")
    extracted = feature.extract(Seq("ACGTACGTAC"))  # type: ignore[no-untyped-call]
    assert str(extracted) == "ACAC"


def test_create_seqfeature_linear_topology_origin_wrapping_location_raises() -> None:
    with pytest.raises(LocationParserError):
        create_seqfeature("misc_feature", "9..2", {}, 10, "linear")


# --- add_dbxrefs ---


def _new_record() -> SeqRecord:
    return SeqRecord(Seq("ATG"), id="x", name="x", description="")


def test_add_dbxrefs_all_fields_present_formats_each_dbxref() -> None:
    record = _new_record()
    add_dbxrefs(record, {
        "DBLINK": {
            "project": "PRJDB1",
            "biosample": "SAMD1",
            "sequence read archive": "DRR1",
        }
    })
    assert record.dbxrefs == ["BioProject:PRJDB1", "BioSample:SAMD1", "Sequence Read Archive:DRR1"]


def test_add_dbxrefs_list_values_joined_with_comma() -> None:
    record = _new_record()
    add_dbxrefs(record, {
        "DBLINK": {
            "project": "PRJDB1",
            "biosample": "SAMD1",
            "sequence read archive": ["DRR1", "DRR2"],
        }
    })
    assert record.dbxrefs[-1] == "Sequence Read Archive:DRR1, DRR2"


def test_add_dbxrefs_no_dblink_results_in_empty_list() -> None:
    record = _new_record()
    add_dbxrefs(record, {})
    assert record.dbxrefs == []


# --- add_comment ---


def test_add_comment_multiple_comment_blocks_joined_with_newline() -> None:
    record = _new_record()
    add_comment(record, {"COMMENT": [{"line": ["a", "b"]}, {"line": ["c"]}]})
    assert record.annotations["comment"] == "a\nb\nc"


def test_add_comment_no_comments_sets_empty_string() -> None:
    record = _new_record()
    add_comment(record, {})
    assert record.annotations["comment"] == ""


def test_add_comment_st_comment_with_tagset_id_sets_structured_comment() -> None:
    record = _new_record()
    add_comment(record, {
        "ST_COMMENT": {
            "tagset_id": "Genome-Assembly-Data",
            "Assembly Method": "HGAP",
            "Genome Coverage": "60x",
        }
    })
    assert record.annotations["structured_comment"] == {  # type: ignore[comparison-overlap]
        "Genome-Assembly-Data": {"Assembly Method": "HGAP", "Genome Coverage": "60x"}
    }


def test_add_comment_st_comment_without_tagset_id_no_structured_comment_key() -> None:
    record = _new_record()
    add_comment(record, {"ST_COMMENT": {}})
    assert "structured_comment" not in record.annotations


# --- set_date ---


def test_set_date_none_uses_todays_date() -> None:
    record = _new_record()
    expected = datetime.now().strftime("%d-%b-%Y").upper()
    set_date(record, None)
    assert record.annotations["date"] == expected


def test_set_date_explicit_value_used_verbatim() -> None:
    record = _new_record()
    set_date(record, "15-MAR-2024")
    assert record.annotations["date"] == "15-MAR-2024"


# --- create_reference と GenBank export の配線 (PUBMED 行) ---


def test_create_reference_pubmed_appears_in_genbank_pubmed_line() -> None:
    ref = create_reference({"title": "T", "pubmed": "12345678", "ab_name": ["Smith,J."]})
    assert ref is not None
    record = SeqRecord(Seq("ATGAAACGCTAA"), id="test1", name="test1", description="test")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["references"] = [ref]  # type: ignore[assignment]

    buf = io.StringIO()
    SeqIO.write([record], buf, "genbank")
    output = buf.getvalue()
    assert "PUBMED   12345678" in output


def test_create_reference_no_pubmed_no_pubmed_line_in_genbank_output() -> None:
    ref = create_reference({"title": "T", "ab_name": ["Smith,J."]})
    assert ref is not None
    record = SeqRecord(Seq("ATGAAACGCTAA"), id="test1", name="test1", description="test")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["references"] = [ref]  # type: ignore[assignment]

    buf = io.StringIO()
    SeqIO.write([record], buf, "genbank")
    output = buf.getvalue()
    assert "PUBMED" not in output
