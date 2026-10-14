import io
import json
from datetime import datetime
from pathlib import Path
from typing import Any, cast

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import LocationParserError, Reference, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from ddbj_record.schema.v3 import (
    Address,
    Organization,
    Person,
    Publication,
    Relation,
    RelationTarget,
    StructuredComment,
)

from dr_tools.json2biopython import (
    add_comment,
    add_dbxrefs,
    add_translate_qualifier,
    create_definition,
    create_reference,
    create_seqfeature,
    create_submitter,
    json_to_seqrecords,
    set_date,
)

# --- JSON builders (直接 ddbj_record v1 スキーマを満たす最小 JSON を組み立てる) ---
# v1 の JSON は v1 -> v2 -> v3 と変換してから読むので、v1 の形のまま入力にして変換の配線も確かめる。


def _build_source_feature(
    *,
    location: str = "1..12",
    ff_definition: str = "@@[organism]@@",
    extra_qualifiers: dict[str, list[str]] | None = None,
) -> dict[str, Any]:
    qualifiers: dict[str, list[str]] = {"ff_definition": [ff_definition]}
    if extra_qualifiers:
        qualifiers.update(extra_qualifiers)
    return {"id": "feature_source", "type": "source", "location": location, "qualifiers": qualifiers}


def _build_cds_feature(
    location: str,
    *,
    product: str = "hypothetical protein",
    transl_table: str = "11",
    codon_start: str = "1",
) -> dict[str, Any]:
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
    features: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
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
    common: dict[str, Any] | None = None,
    common_source: dict[str, Any] | None = None,
    common_meta: dict[str, Any] | None = None,
    entries: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    base_common: dict[str, Any] = {
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

    base_common_source: dict[str, Any] = {
        "organism": "Escherichia coli",
        "mol_type": "genomic DNA",
    }
    if common_source:
        base_common_source.update(common_source)

    base_common_meta: dict[str, Any] = {"division": "BCT"}
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


def _build_v3_record_json(
    *,
    common_source: dict[str, Any] | None = None,
    entries: list[dict[str, Any]] | None = None,
    features: list[dict[str, Any]] | None = None,
    **extra: Any,
) -> dict[str, Any]:
    """v3 の最小の record。entries を省くと、source を持たない entry を 1 つ置く"""
    if common_source is None:
        common_source = {"organism": {"name": "Escherichia coli"}, "mol_type": "genomic DNA"}
    if entries is None:
        entries = [
            {
                "alias": "entry1",
                "name": "entry1",
                "topology": "linear",
                "sequence": "atgaaacgctaa",
                "source_features": [{"location": "1..12", "definition": ["@@[organism]@@"]}],
            }
        ]
    record: dict[str, Any] = {
        "schema_version": "v3",
        "sequences": {"common_source": common_source, "entries": entries},
    }
    if features is not None:
        record["features"] = features
    record.update(extra)
    return record


def _write_record_json(tmp_path: Path, data: dict[str, Any]) -> Path:
    path = tmp_path / "record.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


# --- add_translate_qualifier ---


def test_add_translate_qualifier_valid_coding_sequence_sets_single_element_translation_list() -> None:
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 12, strand=1),  # type: ignore[no-untyped-call]
        type="CDS",
        qualifiers={"transl_table": ["11"], "codon_start": ["1"]},
    )
    add_translate_qualifier(feature, Seq("atgaaacgctaa"))
    assert feature.qualifiers["translation"] == ["MKR"]


def test_add_translate_qualifier_non_nucleotide_characters_no_exception_and_no_translation_qualifier() -> None:
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 12, strand=1),  # type: ignore[no-untyped-call]
        type="CDS",
        qualifiers={"transl_table": ["11"], "codon_start": ["1"]},
    )
    add_translate_qualifier(feature, Seq("atgaaacgc..."))
    assert "translation" not in feature.qualifiers


@pytest.mark.parametrize(
    "qualifiers",
    [
        {"transl_table": ["99"], "codon_start": ["1"]},
        {"transl_table": ["no such table"], "codon_start": ["1"]},
        {"transl_table": ["11"], "codon_start": ["x"]},
    ],
)
def test_add_translate_qualifier_unusable_transl_table_or_codon_start_logs_warning_and_sets_no_translation(
    qualifiers: dict[str, list[str]],
    caplog: pytest.LogCaptureFixture,
) -> None:
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 12, strand=1),  # type: ignore[no-untyped-call]
        type="CDS",
        qualifiers=dict(qualifiers),
    )
    add_translate_qualifier(feature, Seq("atgaaacgctaa"))
    assert "translation" not in feature.qualifiers
    assert "Cannot translate" in caplog.text


def test_json_to_seqrecords_cds_with_unknown_transl_table_keeps_other_features(tmp_path: Path) -> None:
    sequence = "atgaaacgctaaatgaaacgctaa"
    data = _build_record_json(
        entries=[
            _build_entry(
                sequence=sequence,
                features=[
                    _build_source_feature(location=f"1..{len(sequence)}"),
                    _build_cds_feature("1..12", transl_table="99"),
                    _build_cds_feature("13..24"),
                ],
            )
        ]
    )
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    cds_features = [feature for feature in records[0].features if feature.type == "CDS"]
    assert "translation" not in cds_features[0].qualifiers
    assert cds_features[1].qualifiers["translation"] == ["MKR"]


def test_add_translate_qualifier_translation_reduces_to_empty_string_does_not_set_qualifier() -> None:
    feature = SeqFeature(  # type: ignore[no-untyped-call]
        SimpleLocation(0, 3, strand=1),  # type: ignore[no-untyped-call]
        type="CDS",
        qualifiers={"transl_table": ["11"], "codon_start": ["1"]},
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


def test_json_to_seqrecords_entry_placeholder_uses_entry_name(tmp_path: Path) -> None:
    entries = [
        {
            "alias": "contig_001",
            "name": "contig_001",
            "sequence": "atgaaacgctaa",
            "source_features": [{"location": "1..12", "definition": ["@@[organism]@@ DNA, @@[entry]@@"]}],
        }
    ]
    data = _build_v3_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].description == "Escherichia coli DNA, contig_001"


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


def test_json_to_seqrecords_v3_plasmid_placeholder_uses_entry_source_qualifier(tmp_path: Path) -> None:
    entries = [
        {
            "alias": "pX",
            "name": "pX",
            "sequence": "atgaaacgctaa",
            "source_features": [
                {
                    "location": "1..12",
                    "source": {"qualifiers": {"plasmid": [{"value": "pX1"}]}},
                    "definition": ["@@[organism]@@ plasmid @@[plasmid]@@ DNA"],
                }
            ],
        }
    ]
    data = _build_v3_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].description == "Escherichia coli plasmid pX1 DNA"
    source = records[0].features[0]
    assert source.type == "source"
    assert source.qualifiers["plasmid"] == ["pX1"]
    assert source.qualifiers["organism"] == ["Escherichia coli"]


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


def test_json_to_seqrecords_v3_entry_source_mol_type_overrides_common_source(tmp_path: Path) -> None:
    entries = [
        {
            "alias": "e1",
            "name": "e1",
            "sequence": "atgaaacgctaa",
            "source_features": [{"location": "1..12", "source": {"mol_type": "genomic RNA"}}],
        }
    ]
    data = _build_v3_record_json(entries=entries)
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["molecule_type"] == "RNA"


def test_json_to_seqrecords_v3_no_mol_type_defaults_to_dna(tmp_path: Path) -> None:
    data = _build_v3_record_json(common_source={"organism": {"name": "Escherichia coli"}})
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["molecule_type"] == "DNA"


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


def test_json_to_seqrecords_v3_entry_without_division_uses_unk(tmp_path: Path) -> None:
    data = _build_v3_record_json()
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["data_file_division"] == "UNK"


def test_json_to_seqrecords_keywords_present_sets_keywords_annotation(tmp_path: Path) -> None:
    data = _build_record_json(common={"KEYWORD": {"keyword": ["WGS", "STANDARD_DRAFT"]}})
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].annotations["keywords"] == ["WGS", "STANDARD_DRAFT"]  # type: ignore[comparison-overlap]


def test_json_to_seqrecords_keywords_absent_does_not_set_keywords_annotation(tmp_path: Path) -> None:
    data = _build_record_json()
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert "keywords" not in records[0].annotations


def test_json_to_seqrecords_v3_relations_become_dbxrefs(tmp_path: Path) -> None:
    data = _build_v3_record_json(
        relations=[
            {"type": "xref", "target": {"db": "bioproject", "id": "PRJDB1"}},
            {"type": "xref", "target": {"db": "biosample", "id": "SAMD1"}},
        ]
    )
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    assert records[0].dbxrefs == ["BioProject:PRJDB1", "BioSample:SAMD1"]


def test_json_to_seqrecords_v3_publication_journal_fields_reach_reference(tmp_path: Path) -> None:
    data = _build_v3_record_json(
        projects=[
            {
                "publications": [
                    {
                        "title": "T",
                        "status": "published",
                        "date": "2024-05-01",
                        "journal": "J Test",
                        "volume": "3",
                        "pages_from": "10",
                        "pages_to": "12",
                        "pubmed_id": "99999999",
                        "authors": [{"abbreviation": "Smith,J."}],
                    }
                ]
            }
        ]
    )
    records = json_to_seqrecords(_write_record_json(tmp_path, data))
    references = cast(list[Reference], records[0].annotations["references"])
    assert len(references) == 1  # submitters が無いので Direct Submission は無い
    publication = references[0]
    assert publication.journal == "J Test 3:10-12, (2024)"
    assert publication.pubmed_id == "99999999"


# --- json_to_seqrecords: CDS translation wiring (実データ, full_record) ---


def test_json_to_seqrecords_cds_feature_gets_translation_other_feature_types_do_not(
    full_record: dict[str, Path],
) -> None:
    records = json_to_seqrecords(full_record["v3"])
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


def test_json_to_seqrecords_real_complete_genome_v1_entries_and_description(eg_complete: dict[str, Path]) -> None:
    records = json_to_seqrecords(eg_complete["v1_json"])
    assert [r.id for r in records] == ["chromosome", "pPLH-1", "pPLH-2"]
    assert records[0].description == "Paucilactobacillus hokkaidonensis LOOC260 DNA, complete genome"


def test_json_to_seqrecords_real_complete_genome_v3_plasmid_description_and_source_qualifier(
    eg_complete: dict[str, Path],
) -> None:
    records = json_to_seqrecords(eg_complete["v3_json"])
    assert [r.id for r in records] == ["chromosome", "pPLH-1", "pPLH-2"]
    assert records[1].description == "Paucilactobacillus hokkaidonensis LOOC260 plasmid pPLH-1 DNA, complete sequence"
    assert records[1].features[0].qualifiers["plasmid"] == ["pPLH-1"]
    assert records[1].annotations["data_file_division"] == "BCT"
    assert records[1].dbxrefs == ["BioProject:PRJDB99999", "BioSample:SAMD999999", "Sequence Read Archive:DRR999990"]


def test_json_to_seqrecords_real_vrl_v1_description_and_rna_molecule_type(eg_vrl: dict[str, Path]) -> None:
    records = json_to_seqrecords(eg_vrl["v1_json"])
    assert len(records) == 1
    assert records[0].description == (
        "Severe acute respiratory syndrome coronavirus 2 hCov-19/Japan/SZ-NIG-12345/2021 RNA, complete genome"
    )
    assert records[0].annotations["molecule_type"] == "RNA"


def test_json_to_seqrecords_real_vrl_v3_description_division_and_rna_molecule_type(eg_vrl: dict[str, Path]) -> None:
    records = json_to_seqrecords(eg_vrl["v3_json"])
    assert len(records) == 1
    assert records[0].description == (
        "Severe acute respiratory syndrome coronavirus 2 hCov-19/Japan/SZ-NIG-12345/2021 RNA, complete genome"
    )
    assert records[0].annotations["molecule_type"] == "RNA"
    assert records[0].annotations["data_file_division"] == "VRL"


def test_json_to_seqrecords_same_genbank_text_from_v1_v2_v3_of_same_ann(full_record: dict[str, Path]) -> None:
    outputs = []
    for version in ("v1", "v2", "v3"):
        buf = io.StringIO()
        SeqIO.write(json_to_seqrecords(full_record[version]), buf, "genbank")
        outputs.append(buf.getvalue())
    assert outputs[0] == outputs[1] == outputs[2]


def test_json_to_seqrecords_writes_to_genbank_format_without_error(full_record: dict[str, Path]) -> None:
    records = json_to_seqrecords(full_record["v3"])
    buf = io.StringIO()
    SeqIO.write(records, buf, "genbank")
    output = buf.getvalue()
    assert output.startswith("LOCUS")
    assert "ORIGIN" in output


# --- create_submitter ---


def _institution(**kwargs: Any) -> Organization:
    return Organization(type="institution", **kwargs)


def test_create_submitter_none_returns_none() -> None:
    assert create_submitter(None) is None


def test_create_submitter_no_abbreviation_no_consortium_returns_none() -> None:
    assert create_submitter([Person(name="X", organizations=[_institution(name="Y")])]) is None


def test_create_submitter_single_author_sets_authors_verbatim() -> None:
    ref = create_submitter([Person(abbreviation="Smith,J.")])
    assert ref is not None
    assert ref.authors == "Smith,J."


def test_create_submitter_multiple_authors_joins_with_and_before_last() -> None:
    ref = create_submitter(
        [Person(abbreviation="Smith,J."), Person(abbreviation="Doe,A."), Person(abbreviation="Lee,K.")]
    )
    assert ref is not None
    assert ref.authors == "Smith,J., Doe,A. and Lee,K."


def test_create_submitter_consortium_only_without_abbreviation_returns_reference() -> None:
    ref = create_submitter([Person(organizations=[Organization(name="Some Consortium", type="consortium")])])
    assert ref is not None
    assert ref.consrtm == "Some Consortium"
    assert ref.authors == ""


def test_create_submitter_journal_has_contact_affiliation_address_and_url() -> None:
    contact = Person(
        name="John Smith",
        abbreviation="Smith,J.",
        organizations=[
            _institution(
                name="Example Institute",
                department="Dept",
                url="https://example.com",
                address=Address(country="Japan", state="Tokyo", city="Bunkyo", street="1-1-1", postal_code="100-0001"),
            )
        ],
    )
    ref = create_submitter([contact])
    assert ref is not None
    assert (
        ref.journal
        == "Contact:John Smith, Example Institute, Dept;, 1-1-1, Bunkyo, Tokyo 100-0001, Japan URL:https://example.com"
    )


def test_create_submitter_missing_affiliation_uses_placeholders() -> None:
    ref = create_submitter([Person(abbreviation="Smith,J.")])
    assert ref is not None
    assert ref.journal == "Contact:Contact Person, Institute Name, Street Name, City ###-####, Country"


# --- create_reference ---


def test_create_reference_no_title_returns_none() -> None:
    assert create_reference(Publication()) is None


def test_create_reference_status_published_formats_journal_with_volume_and_page_range() -> None:
    ref = create_reference(
        Publication(
            title="T",
            status="published",
            journal="J Test",
            volume="10",
            pages_from="1",
            pages_to="5",
            date="2020",
            authors=[Person(abbreviation="Smith,J.")],
        )
    )
    assert ref is not None
    assert ref.journal == "J Test 10:1-5, (2020)"


def test_create_reference_status_published_without_pages_to_uses_pages_from_only() -> None:
    ref = create_reference(
        Publication(
            title="T",
            status="published",
            journal="J Test",
            volume="10",
            pages_from="42",
            date="2020",
            authors=[Person(abbreviation="Smith,J.")],
        )
    )
    assert ref is not None
    assert ref.journal == "J Test 10:42, (2020)"


def test_create_reference_status_in_press_formats_journal_with_in_press() -> None:
    ref = create_reference(
        Publication(
            title="T", status="in-press", journal="J Test", date="2021", authors=[Person(abbreviation="Smith,J.")]
        )
    )
    assert ref is not None
    assert ref.journal == "J Test (In press, 2021)"


def test_create_reference_status_unpublished_default_formats_journal() -> None:
    ref = create_reference(Publication(title="T", date="2022", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    assert ref.journal == "Unpublished. (2022)"


def test_create_reference_year_is_first_four_characters_of_date() -> None:
    ref = create_reference(Publication(title="T", date="2022-11-30", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    assert ref.journal == "Unpublished. (2022)"


def test_create_reference_pubmed_id_sets_pubmed_id_attribute() -> None:
    ref = create_reference(Publication(title="T", pubmed_id="30000000", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    assert ref.pubmed_id == "30000000"


def test_create_reference_no_pubmed_id_pubmed_id_is_empty_string() -> None:
    ref = create_reference(Publication(title="T", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    assert ref.pubmed_id == ""


def test_create_reference_single_author_authors_set_verbatim() -> None:
    ref = create_reference(Publication(title="T", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    assert ref.authors == "Smith,J."


def test_create_reference_multiple_authors_joined_with_and_before_last() -> None:
    ref = create_reference(
        Publication(
            title="T",
            authors=[Person(abbreviation="Smith,J."), Person(abbreviation="Doe,A."), Person(abbreviation="Lee,K.")],
        )
    )
    assert ref is not None
    assert ref.authors == "Smith,J., Doe,A. and Lee,K."


def test_create_reference_consortiums_joined_into_consrtm() -> None:
    ref = create_reference(Publication(title="T", consortiums=["C1", "C2"]))
    assert ref is not None
    assert ref.consrtm == "C1, C2"


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
    qualifiers = cast(dict[str, list[str | bool]], {"pseudo": [True]})
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


def _xref(db: str, id_: str) -> Relation:
    return Relation(type="xref", target=RelationTarget(db=db, id=id_))


def test_add_dbxrefs_all_dbs_present_formats_each_dbxref() -> None:
    record = _new_record()
    add_dbxrefs(record, [_xref("bioproject", "PRJDB1"), _xref("biosample", "SAMD1"), _xref("insdc.sra", "DRR1")])
    assert record.dbxrefs == ["BioProject:PRJDB1", "BioSample:SAMD1", "Sequence Read Archive:DRR1"]


def test_add_dbxrefs_multiple_ids_of_same_db_joined_with_comma() -> None:
    record = _new_record()
    add_dbxrefs(record, [_xref("bioproject", "PRJDB1"), _xref("insdc.sra", "DRR1"), _xref("insdc.sra", "DRR2")])
    assert record.dbxrefs[-1] == "Sequence Read Archive:DRR1, DRR2"


def test_add_dbxrefs_order_is_fixed_regardless_of_relations_order() -> None:
    record = _new_record()
    add_dbxrefs(record, [_xref("biosample", "SAMD1"), _xref("bioproject", "PRJDB1")])
    assert record.dbxrefs == ["BioProject:PRJDB1", "BioSample:SAMD1"]


def test_add_dbxrefs_no_relations_results_in_empty_list() -> None:
    record = _new_record()
    add_dbxrefs(record, None)
    assert record.dbxrefs == []


def test_add_dbxrefs_unknown_db_and_non_xref_relations_are_ignored() -> None:
    record = _new_record()
    add_dbxrefs(
        record,
        [_xref("assembly", "GCA_1"), Relation(type="part_of", target=RelationTarget(db="bioproject", id="PRJDB1"))],
    )
    assert record.dbxrefs == []


# --- add_comment ---


def test_add_comment_multiple_comments_joined_with_newline() -> None:
    record = _new_record()
    add_comment(record, ["a\nb", "c"], None)
    assert record.annotations["comment"] == "a\nb\nc"


def test_add_comment_no_comments_sets_empty_string() -> None:
    record = _new_record()
    add_comment(record, None, None)
    assert record.annotations["comment"] == ""


def test_add_comment_structured_comment_with_tagset_id_sets_structured_comment() -> None:
    record = _new_record()
    add_comment(
        record,
        None,
        [
            StructuredComment(
                tagset_id="Genome-Assembly-Data",
                fields={"Assembly Method": "HGAP", "Genome Coverage": "60x"},
            )
        ],
    )
    assert record.annotations["structured_comment"] == {  # type: ignore[comparison-overlap]
        "Genome-Assembly-Data": {"Assembly Method": "HGAP", "Genome Coverage": "60x"}
    }


def test_add_comment_multiple_structured_comments_keyed_by_tagset_id() -> None:
    record = _new_record()
    add_comment(
        record,
        None,
        [StructuredComment(tagset_id="A", fields={"x": "1"}), StructuredComment(tagset_id="B", fields={"y": "2"})],
    )
    assert record.annotations["structured_comment"] == {"A": {"x": "1"}, "B": {"y": "2"}}  # type: ignore[comparison-overlap]


def test_add_comment_structured_comment_without_tagset_id_no_structured_comment_key() -> None:
    record = _new_record()
    add_comment(record, None, [StructuredComment(fields={"x": "1"})])
    assert "structured_comment" not in record.annotations


# --- set_date ---


def test_set_date_none_uses_todays_date() -> None:
    record = _new_record()
    expected = datetime.now().astimezone().strftime("%d-%b-%Y").upper()
    set_date(record, None)
    assert record.annotations["date"] == expected


def test_set_date_explicit_value_used_verbatim() -> None:
    record = _new_record()
    set_date(record, "15-MAR-2024")
    assert record.annotations["date"] == "15-MAR-2024"


# --- create_reference と GenBank export の配線 (PUBMED 行) ---


def test_create_reference_pubmed_appears_in_genbank_pubmed_line() -> None:
    ref = create_reference(Publication(title="T", pubmed_id="12345678", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    record = SeqRecord(Seq("ATGAAACGCTAA"), id="test1", name="test1", description="test")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["references"] = [ref]  # type: ignore[assignment]

    buf = io.StringIO()
    SeqIO.write([record], buf, "genbank")
    output = buf.getvalue()
    assert "PUBMED   12345678" in output


def test_create_reference_no_pubmed_no_pubmed_line_in_genbank_output() -> None:
    ref = create_reference(Publication(title="T", authors=[Person(abbreviation="Smith,J.")]))
    assert ref is not None
    record = SeqRecord(Seq("ATGAAACGCTAA"), id="test1", name="test1", description="test")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["references"] = [ref]  # type: ignore[assignment]

    buf = io.StringIO()
    SeqIO.write([record], buf, "genbank")
    output = buf.getvalue()
    assert "PUBMED" not in output
