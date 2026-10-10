import sys
from pathlib import Path

import pytest

from review_final_bam_manifest import main
from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import (
    FINAL_BAM_FILTERING_CONTRACT,
    read_final_bam_manifest,
    sha256_file,
)


def _source_manifest(tmp_path: Path) -> Path:
    bam = tmp_path / "input" / "h3.final.bam"
    bai = tmp_path / "input" / "h3.final.bam.bai"
    bam.parent.mkdir()
    bam.write_bytes(b"bam")
    bai.write_bytes(b"bai")
    manifest = tmp_path / "pending.tsv"
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n"
        f"h3\tchip_histone\teye\ttreatment\tsingle\t{bam}\t{bai}\tdm6\t"
        f"{FINAL_BAM_FILTERING_CONTRACT}\t{sha256_file(bam)}\t"
        f"{sha256_file(bai)}\tpending_review\n"
    )
    return manifest


def test_review_manifest_records_single_end_fragment_length(tmp_path, monkeypatch):
    source = _source_manifest(tmp_path)
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "library_id\tqc_status\testimated_fragment_length_bp\tnotes\n"
        "h3\taccepted\t165\taccepted after QC\n"
    )
    output = tmp_path / "reviewed.tsv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "review_final_bam_manifest.py",
            str(source),
            "--decisions",
            str(decisions),
            "--output",
            str(output),
        ],
    )

    assert main() == 0
    row = read_final_bam_manifest(output)["h3"]
    assert row["qc_status"] == "accepted"
    assert row["estimated_fragment_length_bp"] == "165"


def test_review_manifest_requires_single_end_fragment_length(tmp_path, monkeypatch):
    source = _source_manifest(tmp_path)
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "library_id\tqc_status\testimated_fragment_length_bp\tnotes\n"
        "h3\taccepted\t\taccepted after QC\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "review_final_bam_manifest.py",
            str(source),
            "--decisions",
            str(decisions),
            "--output",
            str(tmp_path / "reviewed.tsv"),
        ],
    )

    with pytest.raises(AcquisitionError, match="requires a positive"):
        main()


def test_review_manifest_accepts_rich_pass_fail_review_table(tmp_path, monkeypatch):
    source = _source_manifest(tmp_path)
    review = tmp_path / "library-review.tsv"
    review.write_text(
        "library_id\tcontext\tfrip\tqc_decision\t"
        "estimated_fragment_length_bp\tnotes\n"
        "h3\teye\t0.12\tPASS\t165\tlooks good\n"
    )
    output = tmp_path / "reviewed.tsv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "review_final_bam_manifest.py",
            str(source),
            "--review-table",
            str(review),
            "--output",
            str(output),
        ],
    )

    assert main() == 0
    row = read_final_bam_manifest(output)["h3"]
    assert row["qc_status"] == "accepted"
    assert row["estimated_fragment_length_bp"] == "165"


def test_review_manifest_selects_one_assay_from_atlas_decision_table(
    tmp_path, monkeypatch
):
    source = _source_manifest(tmp_path)
    decisions = tmp_path / "decisions.tsv"
    decisions.write_text(
        "library_id\tassay\tcontext\tqc_decision\testimated_fragment_length_bp\tnotes\n"
        "atac_rep1\tatac\teye\tpass\t\t\n"
        "h3\th3k27ac\teye\tpass\t170\t\n"
    )
    output = tmp_path / "reviewed.tsv"
    arguments = [
        "review_final_bam_manifest.py",
        str(source),
        "--review-table",
        str(decisions),
        "--output",
        str(output),
    ]
    monkeypatch.setattr(sys, "argv", arguments)
    with pytest.raises(AcquisitionError, match="unknown libraries: atac_rep1"):
        main()

    monkeypatch.setattr(sys, "argv", [*arguments, "--assay", "h3k27ac"])
    assert main() == 0
    row = read_final_bam_manifest(output)["h3"]
    assert row["qc_status"] == "accepted"
    assert row["estimated_fragment_length_bp"] == "170"
