import json
import tempfile
import warnings
from pathlib import Path
from typing import Any

import pytest
from ddbj_record.schema.v1 import Common, CommonMeta, CommonSource
from ddbj_record.schema.v2 import DdbjRecord as V2DdbjRecord
from ddbj_record.schema.v3 import DdbjRecord as V3DdbjRecord
from hypothesis import given, settings
from hypothesis import strategies as st

from dr_tools.ann2json import (
    ann2json_for_dfast,
    build_submitters,
    common_source_qualifiers,
    common_to_dict,
    find_contact_index,
    get_hold_date,
    infer_completeness,
    infer_entry_type,
    infer_meta_from_ann,
    make_entry_dict,
    mss_to_record_v3,
    source_to_dict,
)
from dr_tools.json2ann import json2ann
from dr_tools.MSS import MSS, Entry, Feature, Sequence

Qualifiers = dict[str, list[str | bool]]


def test_common_to_dict_complete(mss_complete: MSS) -> None:
    common_dict = common_to_dict(mss_complete)
    assert "COMMON" in common_dict
    common_obj = common_dict["COMMON"]
    common_obj["trad_submission_category"] = "GNM"  # for test
    Common.model_validate(common_obj)


def test_common_to_dict_vrl(mss_vrl: MSS) -> None:
    common_dict = common_to_dict(mss_vrl)
    assert "COMMON" in common_dict
    common_obj = common_dict["COMMON"]
    common_obj["trad_submission_category"] = "WGS"  # for test
    Common.model_validate(common_obj)


def test_source_to_dict_complete(mss_complete: MSS) -> None:
    source_dict = source_to_dict(mss_complete)
    assert "COMMON_SOURCE" in source_dict
    source_obj = source_dict["COMMON_SOURCE"]
    CommonSource.model_validate(source_obj)


def test_source_to_dict_vrl(mss_vrl: MSS) -> None:
    source_dict = source_to_dict(mss_vrl)
    assert "COMMON_SOURCE" in source_dict
    source_obj = source_dict["COMMON_SOURCE"]
    CommonSource.model_validate(source_obj)


def test_source_to_dict_common_note_becomes_list_and_other_qualifiers_stay_strings() -> None:
    # set の並びによらず、note だけが list になり、organism などは文字列のままであること
    shared: Qualifiers = {
        "organism": ["Escherichia coli"],
        "mol_type": ["genomic DNA"],
        "strain": ["K-12"],
        "note": ["note 1", "note 2"],
    }
    mss = _mss_with_sources(dict(shared), dict(shared))
    common_source = source_to_dict(mss)["COMMON_SOURCE"]
    assert common_source == {
        "organism": "Escherichia coli",
        "mol_type": "genomic DNA",
        "strain": "K-12",
        "note": ["note 1", "note 2"],
    }
    CommonSource.model_validate(common_source)


def test_ann2json_for_dfast_v1_common_note_is_written_as_list(eg_complete: dict[str, Path], tmp_path: Path) -> None:
    lines = []
    for line in eg_complete["ann"].read_text(encoding="utf-8").splitlines():
        lines.append(line)
        if line.split("\t")[1] == "source":
            lines.append("\t\t\tnote\tcommon note")
    ann_file = tmp_path.joinpath("with_note.ann")
    ann_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out_json_file = tmp_path.joinpath("out.json")
    ann2json_for_dfast(ann_file, eg_complete["seq"], out_json_file, record_version="v1")
    common_source = json.loads(out_json_file.read_text(encoding="utf-8"))["COMMON_SOURCE"]
    assert common_source["organism"] == "Paucilactobacillus hokkaidonensis"
    assert common_source["note"] == ["common note"]


def test_infer_meta_from_ann_complete(mss_complete: MSS) -> None:
    common_meta, _seq_info = infer_meta_from_ann(mss_complete)
    assert "COMMON_META" in common_meta
    common_meta_obj = common_meta["COMMON_META"]
    common_meta_obj.pop("trad_submission_category", None)  # not in ANN
    CommonMeta.model_validate(common_meta_obj)


def test_ann2json_for_dfast_complete(eg_complete: dict[str, Path]) -> None:
    for record_version in ("v1", "v2", "v3"):
        ann2json_for_dfast(
            ann_file=eg_complete["ann"],
            seq_file=eg_complete["seq"],
            out_json_file=Path("/dev/null"),
            record_version=record_version,
        )


def test_ann2json_for_dfast_vrl(eg_vrl: dict[str, Path]) -> None:
    for record_version in ("v1", "v2", "v3"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # VRL の contact はどの ab_name にも当たらない
            ann2json_for_dfast(
                ann_file=eg_vrl["ann"],
                seq_file=eg_vrl["seq"],
                out_json_file=Path("/dev/null"),
                record_version=record_version,
            )


def test_ann2_json_with_entry_comment(eg_complete: dict[str, Path]) -> None:
    ann_file = eg_complete["ann"].parent.joinpath("complete_genome_with_comment.ann")
    seq_file = eg_complete["seq"]

    with tempfile.NamedTemporaryFile() as tmp:
        tmp_path = Path(tmp.name)

        ann2json_for_dfast(
            ann_file=ann_file,
            seq_file=seq_file,
            out_json_file=tmp_path,
            record_version="v2",
        )

        json_str = tmp_path.read_text(encoding="utf-8")

    record = V2DdbjRecord.model_validate_json(json_str)
    has_entry_comment = False
    for entry in record.sequences.entries:
        if entry.comments:
            has_entry_comment = True
            break

    assert has_entry_comment


def _count_value_lines(ann_text: str) -> int:
    """5列TSVの5列目 (qualifierの値) が非空の行数を数える"""
    count = 0
    for line in ann_text.splitlines():
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) >= 5 and cols[4] != "":
            count += 1
    return count


def _wgs_mss_with_ff_definition(ff_definition: str) -> MSS:
    """
    WGS の COMMON entry と、ff_definition のみを持つ source feature の entry からなる
    最小の MSS を組み立てる。seq_type の ff_definition 判定が GNM に限定されることを
    確認するために使う (examples の WGS 相当データには plasmid/complete genome 等の
    文言を含む ff_definition が無いため、専用に組み立てる)。
    """
    datatype_qualifiers: Qualifiers = {"type": ["WGS"]}
    common_entry = Entry(
        id="COMMON",
        name="COMMON",
        features=[Feature(type="DATATYPE", id="feature_1", qualifiers=datatype_qualifiers)],
    )
    source_qualifiers: Qualifiers = {"ff_definition": [ff_definition]}
    contig_entry = Entry(
        id="contig01",
        name="contig01",
        features=[Feature(type="source", id="feature_2", location="1..30", qualifiers=source_qualifiers)],
    )
    return MSS(ann_file=None, seq_file=None, entries=[common_entry, contig_entry], sequences=[])


def test_infer_meta_from_ann_wgs_entry_with_topology_feature_sets_seq_topology_circular(
    wgs_mss_files: tuple[Path, Path],
) -> None:
    ann_file, seq_file = wgs_mss_files
    mss = MSS.parse(ann_file, seq_file)
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig01"]["seq_topology"] == "circular"


def test_infer_meta_from_ann_wgs_entry_without_topology_feature_sets_seq_topology_linear(
    wgs_mss_files: tuple[Path, Path],
) -> None:
    ann_file, seq_file = wgs_mss_files
    mss = MSS.parse(ann_file, seq_file)
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig02"]["seq_topology"] == "linear"


def test_infer_meta_from_ann_wgs_ff_definition_plasmid_keyword_seq_type_stays_other() -> None:
    mss = _wgs_mss_with_ff_definition("Escherichia coli DNA, plasmid pXYZ")
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig01"]["seq_type"] == "other"


def test_infer_meta_from_ann_wgs_ff_definition_complete_genome_keyword_seq_type_stays_other() -> None:
    mss = _wgs_mss_with_ff_definition("Escherichia coli DNA, complete genome")
    _common_meta, seq_info = infer_meta_from_ann(mss)
    assert seq_info["contig01"]["seq_type"] == "other"


def test_infer_meta_from_ann_gnm_ff_definition_plasmid_keyword_sets_seq_type_plasmid(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    assert seq_info["pPLH-1"]["seq_type"] == "plasmid"


def test_infer_meta_from_ann_gnm_ff_definition_complete_genome_keyword_sets_seq_type_chromosome(
    mss_complete: MSS,
) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    assert seq_info["chromosome"]["seq_type"] == "chromosome"


def test_infer_meta_from_ann_common_entry_included_in_seq_info(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    assert "COMMON" in seq_info


def test_make_entry_dict_excludes_common_entry_from_output(mss_complete: MSS) -> None:
    _common_meta, seq_info = infer_meta_from_ann(mss_complete)
    entry_dict = make_entry_dict(mss_complete, seq_info)
    entry_ids = [entry["id"] for entry in entry_dict["ENTRIES"]]
    assert "COMMON" not in entry_ids


def test_ann2json_json2ann_roundtrip_wgs_topology_line_count_matches_circular_entries(
    wgs_mss_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("wgs.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    topology_lines = [line for line in out_ann_text.splitlines() if "\tTOPOLOGY\t" in line]
    # wgs_mss_files では contig01 のみ circular (TOPOLOGY 行を持つ)
    assert len(topology_lines) == 1


def test_ann2json_json2ann_roundtrip_wgs_locus_tag_lines_preserved(
    wgs_mss_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("wgs.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    assert "\t\t\tlocus_tag\tTEST_00010" in out_ann_text
    assert "\t\t\tlocus_tag\tTEST_00020" in out_ann_text


def test_ann2json_json2ann_roundtrip_wgs_value_line_count_not_decreased(
    wgs_mss_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("wgs.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    orig_count = _count_value_lines(ann_file.read_text(encoding="utf-8"))
    new_count = _count_value_lines(tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8"))
    assert new_count >= orig_count


def test_ann2json_json2ann_roundtrip_vrl_empty_qualifier_value_line_dropped(
    eg_vrl: dict[str, Path],
    tmp_path: Path,
) -> None:
    json_file = tmp_path.joinpath("vrl.json")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ann2json_for_dfast(eg_vrl["ann"], eg_vrl["seq"], json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    # examples/vrl_result.ann の geo_loc_name は値を持たない (bool qualifierでもない) 唯一の qualifier
    assert "geo_loc_name" not in out_ann_text


def test_ann2json_json2ann_roundtrip_vrl_value_line_count_not_decreased(
    eg_vrl: dict[str, Path],
    tmp_path: Path,
) -> None:
    json_file = tmp_path.joinpath("vrl.json")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ann2json_for_dfast(eg_vrl["ann"], eg_vrl["seq"], json_file, record_version="v2")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    orig_count = _count_value_lines(eg_vrl["ann"].read_text(encoding="utf-8"))
    new_count = _count_value_lines(tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8"))
    assert new_count >= orig_count


# === record v3 ===


def _ann2json_v3(ann_file: Path, seq_file: Path, tmp_path: Path, division: str | None = None) -> dict[str, Any]:
    """ann / seq を v3 の JSON にして dict で返す"""
    json_file = tmp_path.joinpath("record.v3.json")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ann2json_for_dfast(ann_file, seq_file, json_file, division=division, record_version="v3")
    data: dict[str, Any] = json.loads(json_file.read_text(encoding="utf-8"))
    return data


def _mss_signature(mss: MSS) -> list[tuple[str, list[tuple[str, str, tuple[tuple[str, tuple[str | bool, ...]], ...]]]]]:
    """
    MSS を、qualifier の並びに依らずに比べられる形にする
    entry ごとに (name, [(feature type, location, qualifier の (名前, 値の組) を名前順に並べたもの)])
    COMMON entry の feature は record に置き場ごとに分かれて並びが変わるので、feature も並べ直す
    空文字列の値は json2ann が行を書かないので落とす
    """
    signature = []
    for entry in mss.entries:
        features = []
        for feature in entry.features:
            qualifiers = tuple(
                sorted(
                    (key, tuple(value for value in values if value != ""))
                    for key, values in feature.qualifiers.items()
                    if any(value != "" for value in values)
                )
            )
            features.append((feature.type, feature.location, qualifiers))
        if entry.id == "COMMON":
            features.sort()
        signature.append((entry.name, features))
    return signature


def _comparable_with_trad_gnm_fixture(record: dict[str, Any], feature_aliases: set[str]) -> dict[str, Any]:
    """
    record の tests/fixtures/v3/records/trad_gnm.json と比べられる形にする。fixture は手で整えたものなので、次を外す
    - provenance.source_format (fixture は "DFAST"、ann2json は "MSS") と dfast_version (fixture に無い)
    - submitters の first_name / last_name (fixture にだけある。ann2json は名前を分けない)
    - assembly (fixture にだけある。ann からは決められない)
    - entry の sequence (fixture は切り詰めてある)、division と completeness (fixture に無い)
    - entry の source_features[].source (fixture は /plasmid を落とした v2 経路で作られている)
    - features は fixture が 4 件に絞ってあるので、fixture にある alias のものだけ
    """
    data: dict[str, Any] = json.loads(json.dumps(record))
    data["provenance"].pop("source_format", None)
    data["provenance"].pop("dfast_version", None)
    data.pop("assembly", None)
    for person in data["submission"]["submitters"]:
        person.pop("first_name", None)
        person.pop("last_name", None)
    for entry in data["sequences"]["entries"]:
        entry.pop("sequence", None)
        entry.pop("division", None)
        entry.pop("completeness", None)
        for source_feature in entry.get("source_features", []):
            source_feature.pop("source", None)
    data["features"] = {
        feature["alias"]: feature for feature in data["features"] if feature["alias"] in feature_aliases
    }
    return data


def test_ann2json_v3_complete_genome_matches_record_fixture_trad_gnm(
    eg_complete: dict[str, Path],
    trad_gnm_v3_json: Path,
    tmp_path: Path,
) -> None:
    actual = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path)
    expected: dict[str, Any] = json.loads(trad_gnm_v3_json.read_text(encoding="utf-8"))
    feature_aliases = {feature["alias"] for feature in expected["features"]}
    assert _comparable_with_trad_gnm_fixture(actual, feature_aliases) == _comparable_with_trad_gnm_fixture(
        expected, feature_aliases
    )


def test_ann2json_v3_output_is_valid_v3_record_with_source_format_mss(
    eg_complete: dict[str, Path], tmp_path: Path
) -> None:
    data = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path)
    record = V3DdbjRecord.model_validate(data)
    assert record.schema_version == "v3"
    assert record.provenance is not None
    assert record.provenance.source_format == "MSS"
    assert record.provenance.submission_category == "GNM"
    assert record.provenance.model_extra == {"dfast_version": "1.3.4"}


def test_ann2json_v3_plasmid_qualifier_stays_in_entry_source_and_common_qualifiers_do_not(
    eg_complete: dict[str, Path],
    tmp_path: Path,
) -> None:
    data = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path)
    entries = {entry["alias"]: entry for entry in data["sequences"]["entries"]}
    assert entries["pPLH-1"]["source_features"][0]["source"] == {"qualifiers": {"plasmid": [{"value": "pPLH-1"}]}}
    assert "source" not in entries["chromosome"]["source_features"][0]
    assert "plasmid" not in data["sequences"]["common_source"]["qualifiers"]
    assert "strain" in data["sequences"]["common_source"]["qualifiers"]


def test_ann2json_v3_complete_genome_entry_type_completeness_and_division(
    eg_complete: dict[str, Path], tmp_path: Path
) -> None:
    data = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path)
    assert [(e["alias"], e["type"], e["completeness"], e["division"]) for e in data["sequences"]["entries"]] == [
        ("chromosome", "chromosome", "complete", "BCT"),
        ("pPLH-1", "plasmid", "complete", "BCT"),
        ("pPLH-2", "plasmid", "complete", "BCT"),
    ]


def test_ann2json_v3_explicit_division_overrides_default(eg_complete: dict[str, Path], tmp_path: Path) -> None:
    data = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path, division="ENV")
    assert {entry["division"] for entry in data["sequences"]["entries"]} == {"ENV"}


def test_ann2json_v3_without_dfast_comment_division_defaults_to_unk(eg_vrl: dict[str, Path], tmp_path: Path) -> None:
    data = _ann2json_v3(eg_vrl["ann"], eg_vrl["seq"], tmp_path)
    assert data["sequences"]["entries"][0]["division"] == "UNK"
    assert "dfast_version" not in data["provenance"]


def test_ann2json_v3_vrl_value_less_qualifier_has_no_value_and_empty_value_is_kept(
    eg_vrl: dict[str, Path], tmp_path: Path
) -> None:
    data = _ann2json_v3(eg_vrl["ann"], eg_vrl["seq"], tmp_path, division="VRL")
    assert data["sequences"]["common_source"]["mol_type"] == "genomic RNA"
    assert data["sequences"]["common_source"]["qualifiers"]["geo_loc_name"] == [{"value": ""}]
    orf1ab = next(feature for feature in data["features"] if feature["location"].startswith("join(251"))
    assert orf1ab["qualifiers"]["ribosomal_slippage"] == [{}]
    assert data["sequences"]["entries"][0]["division"] == "VRL"
    assert data["sequences"]["entries"][0]["type"] == "segment"


def test_ann2json_v3_vrl_without_division_argument_entry_type_is_chromosome(
    eg_vrl: dict[str, Path], tmp_path: Path
) -> None:
    # division を省くと UNK になり、VRL の決まり (segment) は効かない
    data = _ann2json_v3(eg_vrl["ann"], eg_vrl["seq"], tmp_path)
    assert data["sequences"]["entries"][0]["division"] == "UNK"
    assert data["sequences"]["entries"][0]["type"] == "chromosome"


def test_ann2json_v3_wgs_category_seq_prefix_keywords_and_unplaced_entries(
    wgs_mss_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    data = _ann2json_v3(ann_file, seq_file, tmp_path)
    assert data["provenance"]["submission_category"] == "WGS"
    assert data["sequences"]["seq_prefix"] == "contig"
    assert data["sequences"]["keywords"] == ["WGS"]
    assert [(e["alias"], e["type"], e["topology"]) for e in data["sequences"]["entries"]] == [
        ("contig01", "unplaced", "circular"),
        ("contig02", "unplaced", "linear"),
    ]
    assert "completeness" not in data["sequences"]["entries"][0]
    assert data["projects"][0]["locus_tag_prefix"] == [{"prefix": "TEST"}]
    assert [feature["locus_tag_id"] for feature in data["features"]] == ["00010", "00020"]


def test_ann2json_v3_submitters_contact_first_with_organization_and_others_with_abbreviation_only(
    eg_complete: dict[str, Path],
    tmp_path: Path,
) -> None:
    data = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path)
    submitters = data["submission"]["submitters"]
    assert [(p.get("abbreviation"), p["role"]) for p in submitters] == [
        ("Mishima,H.", "contact"),
        ("Tanizawa,Y.", "submitter"),
        ("Nakamura,Y.", "submitter"),
    ]
    assert submitters[0]["name"] == "Hanako Mishima"
    assert submitters[0]["organizations"][0]["type"] == "institution"
    assert submitters[0]["organizations"][0]["address"]["postal_code"] == "411-8540"
    assert "organizations" not in submitters[1]


def test_ann2json_v3_hold_date_and_publication_date(eg_complete: dict[str, Path], tmp_path: Path) -> None:
    data = _ann2json_v3(eg_complete["ann"], eg_complete["seq"], tmp_path)
    assert data["submission"]["hold_date"] == "2025-03-31"
    assert data["projects"][0]["publications"][0]["date"] == "2023"
    assert data["projects"][0]["publications"][0]["status"] == "unpublished"


def test_ann2json_v3_entry_comment_is_kept_in_entry_comments(eg_complete: dict[str, Path], tmp_path: Path) -> None:
    ann_file = eg_complete["ann"].parent.joinpath("complete_genome_with_comment.ann")
    data = _ann2json_v3(ann_file, eg_complete["seq"], tmp_path)
    assert data["sequences"]["entries"][0]["comments"] == ["comment at entry for test"]
    assert all(feature["type"] != "COMMENT" for feature in data["features"])


# --- ann -> v3 -> ann の往復 ---


def test_ann2json_json2ann_roundtrip_v3_complete_genome_keeps_plasmid_line(
    eg_complete: dict[str, Path], tmp_path: Path
) -> None:
    json_file = tmp_path.joinpath("v3.json")
    ann2json_for_dfast(eg_complete["ann"], eg_complete["seq"], json_file, record_version="v3")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    out_ann_text = tmp_path.joinpath("roundtrip.ann").read_text(encoding="utf-8")
    assert "\t\t\tplasmid\tpPLH-1\n" in out_ann_text
    assert "\t\t\tplasmid\tpPLH-2\n" in out_ann_text


@pytest.mark.parametrize("example", ["complete", "vrl"])
def test_ann2json_json2ann_roundtrip_v3_keeps_all_entries_features_and_qualifiers(
    example: str,
    eg_complete: dict[str, Path],
    eg_vrl: dict[str, Path],
    tmp_path: Path,
) -> None:
    files = eg_complete if example == "complete" else eg_vrl
    json_file = tmp_path.joinpath("v3.json")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ann2json_for_dfast(files["ann"], files["seq"], json_file, record_version="v3")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    original = MSS.parse(files["ann"], files["seq"])
    roundtrip = MSS.parse(tmp_path.joinpath("roundtrip.ann"), tmp_path.joinpath("roundtrip.fasta"))
    assert _mss_signature(roundtrip) == _mss_signature(original)
    assert [(s.id, s.seq) for s in roundtrip.sequences] == [(s.id, s.seq) for s in original.sequences]


def test_ann2json_json2ann_roundtrip_v3_wgs_keeps_datatype_keyword_and_topology(
    wgs_mss_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    ann_file, seq_file = wgs_mss_files
    json_file = tmp_path.joinpath("v3.json")
    ann2json_for_dfast(ann_file, seq_file, json_file, record_version="v3")
    json2ann(json_file, out_dir=tmp_path, out_prefix="roundtrip")

    roundtrip = MSS.parse(tmp_path.joinpath("roundtrip.ann"), tmp_path.joinpath("roundtrip.fasta"))
    assert _mss_signature(roundtrip) == _mss_signature(MSS.parse(ann_file, seq_file))


def test_json2ann_v1_v2_v3_of_same_ann_write_same_ann_and_fasta(full_record: dict[str, Path], tmp_path: Path) -> None:
    outputs = []
    for version in ("v1", "v2", "v3"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            json2ann(full_record[version], out_dir=tmp_path, out_prefix=version)
        outputs.append(
            (
                tmp_path.joinpath(f"{version}.ann").read_text(encoding="utf-8"),
                tmp_path.joinpath(f"{version}.fasta").read_text(encoding="utf-8"),
            )
        )
    assert outputs[0] == outputs[1] == outputs[2]


# --- infer_entry_type ---


@pytest.mark.parametrize(
    ("qualifiers", "category", "division", "expected"),
    [
        ({"plasmid": ["pX"]}, "GNM", "BCT", "plasmid"),
        ({"chromosome": ["1"]}, "GNM", "BCT", "chromosome"),
        ({"organelle": ["mitochondrion"]}, "GNM", "BCT", "organelle"),
        ({"segment": ["4"]}, "GNM", "VRL", "segment"),
        ({"plasmid": ["pX"], "chromosome": ["1"]}, "GNM", "BCT", "plasmid"),
        ({"organelle": ["plastid"], "segment": ["1"]}, "WGS", "BCT", "organelle"),
        ({"submitter_seqid": ["contig01"]}, "WGS", "BCT", "unplaced"),
        ({}, "WGS", "BCT", "unplaced"),
        ({}, "WGS", "VRL", "unplaced"),
        ({}, "GNM", "BCT", "chromosome"),
        ({}, "GNM", "UNK", "chromosome"),
        ({}, "GNM", None, "chromosome"),
        ({}, "GNM", "VRL", "segment"),
        ({"isolate": ["x"]}, "GNM", "VRL", "segment"),
        ({"plasmid": ["pX"]}, "GNM", "VRL", "plasmid"),
        ({"plasmid": ["pX"]}, "WGS", "BCT", "plasmid"),
    ],
)
def test_infer_entry_type(qualifiers: dict[str, list[str]], category: str, division: str | None, expected: str) -> None:
    assert infer_entry_type(qualifiers, category, division) == expected


def test_infer_entry_type_division_defaults_to_none() -> None:
    assert infer_entry_type({}, "GNM") == "chromosome"


# --- infer_completeness ---


@pytest.mark.parametrize(
    ("ff_definition", "expected"),
    [
        ("@@[organism]@@ @@[strain]@@ DNA, complete genome", "complete"),
        ("@@[organism]@@ plasmid @@[plasmid]@@ DNA, complete sequence", "complete"),
        ("@@[organism]@@ DNA, nearly complete genome", "nearly complete"),
        ("@@[organism]@@ DNA, partial genome", "partial"),
        ("@@[organism]@@ DNA, partial sequence", "partial"),
        ("@@[organism]@@ DNA, draft genome", "partial"),
        ("@@[organism]@@ DNA, draft sequence", "partial"),
        ("@@[organism]@@ DNA, @@[submitter_seqid]@@", None),
        ("Escherichia coli DNA, Complete Genome", "complete"),
        ("", None),
        (None, None),
    ],
)
def test_infer_completeness(ff_definition: str | None, expected: str | None) -> None:
    assert infer_completeness(ff_definition) == expected


@given(text=st.text(max_size=60))
def test_infer_completeness_pbt_returns_known_value_or_none(text: str) -> None:
    assert infer_completeness(text) in (None, "complete", "nearly complete", "partial")


# --- build_submitters / find_contact_index ---


def _common_entry_with_submitter(qualifiers: Qualifiers) -> Entry:
    return Entry(
        id="COMMON", name="COMMON", features=[Feature(type="SUBMITTER", id="feature_1", qualifiers=qualifiers)]
    )


def test_build_submitters_contact_matching_later_ab_name_is_moved_to_front() -> None:
    entry = _common_entry_with_submitter(
        {
            "ab_name": ["Doe,A.", "Smith,J."],
            "contact": ["John Smith"],
            "email": ["smith@example.com"],
            "institute": ["Inst"],
        }
    )
    submitters = build_submitters(entry)
    assert submitters is not None
    assert [(p.abbreviation, p.role, p.name) for p in submitters] == [
        ("Smith,J.", "contact", "John Smith"),
        ("Doe,A.", "submitter", None),
    ]
    assert submitters[0].email == "smith@example.com"
    assert submitters[0].organizations is not None
    assert submitters[0].organizations[0].name == "Inst"
    assert submitters[1].organizations is None


def test_build_submitters_contact_without_matching_ab_name_adds_person_without_abbreviation_and_warns() -> None:
    entry = _common_entry_with_submitter({"ab_name": ["Doe,A."], "contact": ["John Smith"], "email": ["s@example.com"]})
    with pytest.warns(UserWarning, match="John Smith"):
        submitters = build_submitters(entry)
    assert submitters is not None
    assert [(p.abbreviation, p.name, p.role) for p in submitters] == [
        (None, "John Smith", "contact"),
        ("Doe,A.", None, "submitter"),
    ]


def test_build_submitters_no_contact_assigns_email_to_first_ab_name() -> None:
    entry = _common_entry_with_submitter({"ab_name": ["Doe,A.", "Lee,K."], "email": ["doe@example.com"]})
    submitters = build_submitters(entry)
    assert submitters is not None
    assert submitters[0].abbreviation == "Doe,A."
    assert submitters[0].email == "doe@example.com"
    assert submitters[0].name is None


def test_build_submitters_no_ab_name_no_contact_with_email_returns_single_person() -> None:
    entry = _common_entry_with_submitter({"email": ["x@example.com"]})
    submitters = build_submitters(entry)
    assert submitters is not None
    assert len(submitters) == 1
    assert submitters[0].email == "x@example.com"
    assert submitters[0].role == "contact"


def test_build_submitters_empty_submitter_feature_returns_none() -> None:
    assert build_submitters(_common_entry_with_submitter({})) is None


def test_build_submitters_no_submitter_feature_returns_none() -> None:
    assert build_submitters(Entry(id="COMMON", name="COMMON")) is None


def test_build_submitters_consrtm_becomes_consortium_organization_after_institution() -> None:
    entry = _common_entry_with_submitter(
        {
            "ab_name": ["Smith,J."],
            "contact": ["John Smith"],
            "institute": ["Inst"],
            "consrtm": ["Big Consortium"],
        }
    )
    submitters = build_submitters(entry)
    assert submitters is not None
    assert submitters[0].organizations is not None
    assert [(o.name, o.type) for o in submitters[0].organizations] == [
        ("Inst", "institution"),
        ("Big Consortium", "consortium"),
    ]


def test_build_submitters_address_only_institution_is_kept() -> None:
    entry = _common_entry_with_submitter({"ab_name": ["Smith,J."], "contact": ["John Smith"], "country": ["Japan"]})
    submitters = build_submitters(entry)
    assert submitters is not None
    assert submitters[0].organizations is not None
    assert submitters[0].organizations[0].address is not None
    assert submitters[0].organizations[0].address.country == "Japan"
    assert submitters[0].organizations[0].name is None


@pytest.mark.parametrize(
    ("contact", "ab_names", "expected"),
    [
        ("Hanako Mishima", ["Mishima,H.", "Tanizawa,Y."], 0),
        ("Hanako Mishima", ["Tanizawa,Y.", "Mishima,H."], 1),
        ("Mishima, Hanako", ["Mishima,H."], 0),
        ("Yamada Taro", ["Yamada,T."], 0),
        ("John A. Doe", ["Doe,J.A."], 0),
        ("John A. Doe", ["Doe,JA"], 0),
        ("Masanori Arita", ["Fujisawa,T.", "Tanizawa,Y."], None),
        (None, ["Mishima,H."], None),
        ("", ["Mishima,H."], None),
        ("Hanako Mishima", [], None),
    ],
)
def test_find_contact_index(contact: str | None, ab_names: list[str], expected: int | None) -> None:
    assert find_contact_index(contact, ab_names) == expected


# --- get_hold_date ---


@pytest.mark.parametrize(
    ("hold_date", "expected"), [("20250331", "2025-03-31"), ("2025-03-31", "2025-03-31"), ("202503", "202503")]
)
def test_get_hold_date(hold_date: str, expected: str) -> None:
    entry = Entry(
        id="COMMON",
        name="COMMON",
        features=[Feature(type="DATE", id="feature_1", qualifiers={"hold_date": [hold_date]})],
    )
    assert get_hold_date(entry) == expected


def test_get_hold_date_no_date_feature_returns_none() -> None:
    assert get_hold_date(Entry(id="COMMON", name="COMMON")) is None


# --- common_source_qualifiers ---


def _mss_with_sources(*qualifiers: Qualifiers) -> MSS:
    entries = [Entry(id="COMMON", name="COMMON")]
    sequences = []
    for i, q in enumerate(qualifiers):
        entry_id = f"entry{i + 1}"
        entries.append(
            Entry(
                id=entry_id,
                name=entry_id,
                features=[
                    Feature(type="source", id=f"feature_{i + 1}", location="1..10", qualifiers=q),
                ],
            )
        )
        sequences.append(Sequence(id=entry_id, seq="acgtacgtac"))
    return MSS(ann_file=None, seq_file=None, entries=entries, sequences=sequences)


def test_common_source_qualifiers_only_qualifiers_equal_in_all_entries_are_common() -> None:
    mss = _mss_with_sources(
        {"organism": ["E. coli"], "strain": ["A"], "note": ["n1", "n2"], "plasmid": ["pX"]},
        {"organism": ["E. coli"], "strain": ["B"], "note": ["n1", "n2"]},
    )
    assert common_source_qualifiers(mss) == {"organism": ["E. coli"], "note": ["n1", "n2"]}


def test_common_source_qualifiers_submitter_seqid_and_ff_definition_are_never_common() -> None:
    mss = _mss_with_sources(
        {"organism": ["E. coli"], "submitter_seqid": ["x"], "ff_definition": ["@@[organism]@@"]},
        {"organism": ["E. coli"], "submitter_seqid": ["x"], "ff_definition": ["@@[organism]@@"]},
    )
    assert common_source_qualifiers(mss) == {"organism": ["E. coli"]}


def test_common_source_qualifiers_does_not_modify_the_mss() -> None:
    mss = _mss_with_sources({"organism": ["E. coli"], "strain": ["A"]}, {"organism": ["E. coli"], "strain": ["A"]})
    common_source_qualifiers(mss)
    assert mss.entries[1].features[0].qualifiers == {"organism": ["E. coli"], "strain": ["A"]}


def test_common_source_qualifiers_no_source_features_returns_empty_dict() -> None:
    assert common_source_qualifiers(_mss_with_sources()) == {}


def test_mss_to_record_v3_entry_without_sequence_raises_value_error() -> None:
    mss = _mss_with_sources({"organism": ["E. coli"]})
    mss.sequences = []
    with pytest.raises(ValueError, match="entry1"):
        mss_to_record_v3(mss)


def test_mss_to_record_v3_locus_tag_without_underscore_raises_value_error() -> None:
    mss = _mss_with_sources({"organism": ["E. coli"]})
    mss.entries[1].features.append(
        Feature(type="CDS", id="feature_9", location="1..9", qualifiers={"locus_tag": ["ABC123"]})
    )
    with pytest.raises(ValueError, match="ABC123"):
        mss_to_record_v3(mss)


def test_mss_to_record_v3_without_common_entry_raises_value_error() -> None:
    mss = _mss_with_sources({"organism": ["E. coli"]})
    mss.entries = mss.entries[1:]
    with pytest.raises(ValueError, match="COMMON entry not found"):
        mss_to_record_v3(mss)


# --- PBT: MSS -> v3 -> MSS の往復 ---

_text = st.text(alphabet=st.characters(blacklist_categories=("Cs", "Cc")), min_size=1, max_size=6)
_qualifier_value = st.one_of(st.just(True), _text)
_source_qualifiers = st.dictionaries(
    st.sampled_from(["strain", "note", "isolate", "plasmid", "submitter_seqid", "pseudo", "collection_date"]),
    st.lists(_qualifier_value, min_size=1, max_size=2),
    max_size=5,
)
_feature_qualifiers = st.dictionaries(
    st.sampled_from(["product", "note", "pseudo", "gene"]),
    st.lists(_qualifier_value, min_size=1, max_size=2),
    max_size=3,
)
_entry_name = st.from_regex(r"[A-Za-z0-9_.-]{1,8}", fullmatch=True).filter(lambda name: name != "COMMON")


@st.composite
def _mss_strategy(draw: st.DrawFn) -> MSS:
    entries = [
        Entry(
            id="COMMON",
            name="COMMON",
            features=[
                Feature(type="DATATYPE", id="feature_0", qualifiers={"type": ["WGS"]}),
            ]
            if draw(st.booleans())
            else [],
        )
    ]
    sequences = []
    names = draw(st.lists(_entry_name, min_size=1, max_size=3, unique=True))
    counter = 1
    for name in names:
        source_qualifiers: Qualifiers = {}
        if draw(st.booleans()):
            source_qualifiers["organism"] = [draw(_text)]
        if draw(st.booleans()):
            source_qualifiers["mol_type"] = [draw(_text)]
        source_qualifiers.update(draw(_source_qualifiers))
        if draw(st.booleans()):
            source_qualifiers["ff_definition"] = draw(st.lists(_text, min_size=1, max_size=2))
        features = [Feature(type="source", id=f"feature_{counter}", location="1..10", qualifiers=source_qualifiers)]
        counter += 1
        if draw(st.booleans()):
            features.append(Feature(type="TOPOLOGY", id=f"feature_{counter}", qualifiers={"circular": [True]}))
            counter += 1
        for feature_type in draw(st.lists(st.sampled_from(["CDS", "gene", "rRNA"]), max_size=2)):
            qualifiers: Qualifiers = draw(_feature_qualifiers)
            if draw(st.booleans()):
                qualifiers["locus_tag"] = [f"PFX_{counter:05d}"]
            features.append(
                Feature(type=feature_type, id=f"feature_{counter}", location="1..10", qualifiers=qualifiers)
            )
            counter += 1
        entries.append(Entry(id=name, name=name, features=features))
        sequences.append(Sequence(id=name, seq=draw(st.text(alphabet="acgtn", max_size=12))))
    return MSS(ann_file=None, seq_file=None, entries=entries, sequences=sequences)


@settings(max_examples=60, deadline=None)
@given(mss=_mss_strategy())
def test_mss_to_record_v3_and_back_keeps_entries_source_qualifiers_features_and_sequences(mss: MSS) -> None:
    record = mss_to_record_v3(mss)
    # JSON に書いて読み直しても同じ record になる
    reloaded = V3DdbjRecord.model_validate_json(record.model_dump_json(exclude_none=True))
    assert reloaded.model_dump(exclude_none=True) == record.model_dump(exclude_none=True)

    roundtrip = MSS.from_record(reloaded)
    assert [entry.name for entry in roundtrip.entries] == [entry.name for entry in mss.entries]
    assert [(s.id, s.seq) for s in roundtrip.sequences] == [(s.id, s.seq) for s in mss.sequences]
    for original_entry, roundtrip_entry in zip(mss.entries[1:], roundtrip.entries[1:], strict=True):
        original_features = [(f.type, f.location, f.qualifiers) for f in original_entry.features]
        roundtrip_features = [(f.type, f.location, f.qualifiers) for f in roundtrip_entry.features]
        # source と TOPOLOGY は並びが決まっている。source の qualifier は名前ごとの比較 (dict) なので並びを問わない
        assert roundtrip_features == original_features
