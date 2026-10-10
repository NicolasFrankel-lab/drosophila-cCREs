from collections import defaultdict
from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.sample_sheet import read_sample_sheet


def test_four_column_csv_is_accepted_and_defaults_are_resolved(tmp_path: Path):
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "accession,library_id,assay,context\n"
        "SRR123,atac_rep1,atac,eye\n"
    )
    row = read_sample_sheet(sheet)[0]

    assert row["accession"] == "SRR123"
    assert row["library_id"] == "atac_rep1"
    assert row["role"] == "treatment"
    assert row["peak_caller"] == "callpeak"
    assert row["macs3_qvalue"] == 0.1
    assert row["macs3_shift"] == -75
    assert row["macs3_extsize"] == 150
    assert row["adapter_preset"] == "nextera"
    assert row["mapq_minimum"] == 30


def test_h3k27ac_alias_resolves_histone_defaults(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR123\th3_rep1\th3k27ac\teye\n"
    )
    row = read_sample_sheet(sheet)[0]

    assert row["assay"] == "chip_histone"
    assert row["macs3_broad"] is True
    assert row["macs3_broad_cutoff"] == 0.1
    assert row["adapter_preset"] == "truseq"


def test_missing_required_column_fails(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\n"
        "SRR123\tatac_rep1\tatac\n"
    )
    with pytest.raises(AcquisitionError, match="missing required columns: context"):
        read_sample_sheet(sheet)


def test_technical_runs_must_have_identical_contexts(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR123\tatac_rep1\tatac\teye\n"
        "SRR124\tatac_rep1\tatac\twing\n"
    )
    with pytest.raises(AcquisitionError, match="disagree on: context"):
        read_sample_sheet(sheet)


def test_h3k27ac_callpeak_shift_requires_bam_format(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\tmacs3_shift\tmacs3_extsize\n"
        "SRR123\th3_rep1\th3k27ac\teye\t-75\t150\n"
    )
    with pytest.raises(AcquisitionError, match="requires macs3_format=BAM"):
        read_sample_sheet(sheet)


@pytest.mark.parametrize("assay", ["chip_tf", "chip_histone"])
def test_only_atac_and_h3k27ac_assays_are_accepted(tmp_path: Path, assay):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        f"SRR123\tlib_rep1\t{assay}\teye\n"
    )
    with pytest.raises(AcquisitionError, match="assay must be one of atac, h3k27ac"):
        read_sample_sheet(sheet)


def test_h3k27ac_requires_broad_peaks(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\tmacs3_broad\n"
        "SRR123\th3_rep1\th3k27ac\teye\tfalse\n"
    )
    with pytest.raises(AcquisitionError, match="macs3_broad=true"):
        read_sample_sheet(sheet)


def test_atlas_sample_sheet_is_the_eight_context_selection():
    resources = Path(__file__).parents[1] / "resources"
    rows = read_sample_sheet(resources / "atlas_samples.tsv")
    contexts = {"ab", "e13", "e5", "ead", "hid", "lb", "o", "wid"}

    assert len(rows) == 35
    assert len({row["accession"] for row in rows}) == 35
    assert [row["library_id"] for row in rows].count("e5_atac_rep1") == 2
    libraries: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        libraries[(str(row["assay"]), str(row["context"]))].add(str(row["library_id"]))
    assert {context for _assay, context in libraries} == contexts
    assert sum(len(values) for (assay, _), values in libraries.items() if assay == "atac") == 20
    assert (
        sum(len(values) for (assay, _), values in libraries.items() if assay == "chip_histone")
        == 14
    )
    assert all(len(libraries[("atac", context)]) >= 2 for context in contexts)
    assert all(libraries[("chip_histone", context)] for context in contexts)


def test_atlas_qc_decisions_cover_every_atlas_library_once():
    resources = Path(__file__).parents[1] / "resources"
    rows = read_sample_sheet(resources / "atlas_samples.tsv")
    expected = {
        str(row["library_id"]): (
            "h3k27ac" if row["assay"] == "chip_histone" else str(row["assay"]),
            str(row["context"]),
        )
        for row in rows
    }
    lines = (resources / "atlas_qc_decisions.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    decisions = [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    assert header == [
        "library_id",
        "assay",
        "context",
        "qc_decision",
        "estimated_fragment_length_bp",
        "notes",
    ]
    assert len(decisions) == len({row["library_id"] for row in decisions})
    assert {
        row["library_id"]: (row["assay"], row["context"]) for row in decisions
    } == expected
    assert {row["qc_decision"] for row in decisions} == {"pass"}
    assert all(
        not row["estimated_fragment_length_bp"]
        for row in decisions
        if row["assay"] == "atac"
    )
