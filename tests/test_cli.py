import argparse
import json
from pathlib import Path
import sys

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import sha256_file
from short_read_processing.cli import read_accession_column
import run_pipeline
from run_pipeline import main as run_pipeline_main, rule_threads


def test_read_configurable_accession_column(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("run_id\tlabel\nSRR123\ta\nERX456\tb\nSRR123\tduplicate\n")
    assert read_accession_column(metadata, "run_id") == ["SRR123", "ERX456"]


def test_read_accession_column_from_csv(tmp_path):
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("run_id,label\nSRR123,a\nERR456,b\n")
    assert read_accession_column(metadata, "run_id") == ["SRR123", "ERR456"]


def test_missing_accession_column_is_clear(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("wrong\nSRR123\n")
    with pytest.raises(AcquisitionError, match="available columns: wrong"):
        read_accession_column(metadata, "run_id")


@pytest.mark.parametrize("value", ["align_lane", "align_lane=0", "=4", "align_lane=x"])
def test_rule_thread_override_rejects_invalid_values(value):
    with pytest.raises(argparse.ArgumentTypeError, match="RULE=THREADS"):
        rule_threads(value)


def test_rule_thread_override_accepts_positive_count():
    assert rule_threads("align_lane=16") == "align_lane=16"


def test_alignment_chunk_size_must_be_positive(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--alignment-chunk-pairs",
            "0",
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_final_bam_stage_requires_explicit_manifest(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_pipeline.py", str(sheet), "--from-stage", "final-bam"],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_alignment_boundary_requires_one_reuse_manifest(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "alignment",
            "--until-stage",
            "qc",
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_alignment_checkpoint_supplies_final_bam_manifest(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    final_bams = tmp_path / "final-bams.tsv"
    final_bams.write_text("manifest\n")
    checkpoint = tmp_path / "alignment.checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "stage": "alignment",
                "source_project": "atlas",
                "source_run_id": "alignment-v1",
                "semantic_sha256": "a" * 64,
                "parameters": {},
                "artifacts": {
                    "final_bam_manifest": {
                        "path": final_bams.name,
                        "sha256": sha256_file(final_bams),
                    }
                },
            }
        )
        + "\n"
    )
    captured = {}

    def fake_generate_configs(**kwargs):
        captured.update(kwargs)
        output = tmp_path / "resolved.yaml"
        output.write_text("project: atlas\n")
        return [output]

    monkeypatch.setattr(run_pipeline, "generate_configs", fake_generate_configs)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "alignment",
            "--checkpoint-manifest",
            str(checkpoint),
            "--until-stage",
            "qc",
            "--config-only",
        ],
    )

    assert run_pipeline_main() == 0
    assert captured["input_stage"] == "final-bam"
    assert captured["start_stage"] == "alignment"
    assert captured["final_bam_manifest_path"] == final_bams.resolve()


def test_reuse_stage_rejects_download_fallback_flags(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    manifest = tmp_path / "final-bams.tsv"
    manifest.write_text("")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "final-bam",
            "--final-bam-manifest",
            str(manifest),
            "--skip-download",
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_quantification_stage_requires_all_explicit_artifact_inputs(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "quantification",
            "--master-manifest",
            str(tmp_path / "master.tsv"),
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_accession_stage_rejects_activity_artifacts(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--activity-bam-manifest",
            str(tmp_path / "activity.tsv"),
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_accession_stage_requires_qc_review_before_master(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_pipeline.py", str(sheet), "--until-stage", "master"],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_nonactivity_reuse_stage_rejects_activity_options(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "master",
            "--master-manifest",
            str(tmp_path / "master.tsv"),
            "--activity-bam-manifest",
            str(tmp_path / "activity.tsv"),
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_final_bam_input_rejects_trimming_output_stage(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "final-bam",
            "--until-stage",
            "trimming",
            "--final-bam-manifest",
            str(tmp_path / "final-bams.tsv"),
        ],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_activity_threshold_reaches_activity_configuration(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    master = tmp_path / "master.tsv"
    master.write_text("placeholder\n")
    activity = tmp_path / "activity.tsv"
    activity.write_text("placeholder\n")
    captured = {}

    def fake_generate_activity_config(**kwargs):
        captured.update(kwargs)
        output = tmp_path / "resolved.yaml"
        output.write_text("project: atlas\n")
        return output

    monkeypatch.setattr(
        run_pipeline, "generate_activity_config", fake_generate_activity_config
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            str(sheet),
            "--from-stage",
            "master",
            "--until-stage",
            "report",
            "--master-manifest",
            str(master),
            "--activity-bam-manifest",
            str(activity),
            "--h3k27ac-activity-threshold",
            "0.75",
            "--config-only",
        ],
    )

    assert run_pipeline_main() == 0
    assert captured["h3k27ac_activity_threshold"] == 0.75


@pytest.mark.parametrize(
    "arguments",
    [
        ["--h3k27ac-activity-threshold", "1.0"],
        ["--h3k27ac-activity-threshold", "0.7"],
    ],
)
def test_activity_threshold_is_rejected_when_invalid_or_unused(
    monkeypatch, tmp_path, arguments
):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    monkeypatch.setattr(sys, "argv", ["run_pipeline.py", str(sheet), *arguments])

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2


def test_footprints_require_master_and_reviewed_bam_manifests(monkeypatch, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    master = tmp_path / "master.tsv"
    master.write_text("")
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_pipeline.py", str(sheet), "--from-stage", "master", "--until-stage", "footprints",
         "--master-manifest", str(master)],
    )

    with pytest.raises(SystemExit) as error:
        run_pipeline_main()

    assert error.value.code == 2
