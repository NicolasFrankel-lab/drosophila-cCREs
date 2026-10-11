import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import sha256_file
from short_read_processing.configuration import ATAC_QPOIS_DEFAULTS, REFERENCE_SOURCES
from short_read_processing.workflow_config import (
    guard_result_namespace,
    validate_stage_selection,
    validate_workflow_config,
    workflow_semantic_sha256,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = yaml.safe_load((REPO_ROOT / "tests/fixtures/workflow_config.yaml").read_text())


def _h3k27ac_callpeak():
    return {
        "command": "callpeak",
        "format": "BAMPE",
        "qvalue": 0.01,
        "broad": True,
        "broad_cutoff": 0.1,
        "nomodel": False,
        "shift": None,
        "extsize": None,
        "write_bedgraph": True,
        "spmr": True,
    }


def _dry_run(tmp_path: Path, config: dict, name: str = "workflow") -> str:
    config["output_dir"] = str(tmp_path / "results")
    config_path = tmp_path / f"{name}.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    snakemake = Path(sys.executable).with_name("snakemake")
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    result = subprocess.run(
        [
            str(snakemake),
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "8",
            "--dry-run",
            "--printshellcmds",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    return output


def _activity_config(tmp_path: Path, contexts=("ctx",)) -> dict:
    config = copy.deepcopy(BASE_CONFIG)
    config["assay"] = "activity"
    config["input_stage"] = "quantification"
    config["output_stage"] = "catalog"
    config["samples"] = []
    config.pop("atac_qpois")
    master = {}
    for field in (
        "master_bed",
        "summits_bed",
        "membership_tsv",
        "context_matrix_tsv",
        "stats_json",
    ):
        path = tmp_path / "master" / field
        path.parent.mkdir(exist_ok=True)
        path.write_text(field)
        master[field] = str(path)
        master[f"{field}_sha256"] = "a" * 64
    master.update(
        {
            "genome": "dm6",
            "method": "reciprocal_summit_complete_linkage_v2",
            "input_filtering_contract": "short-read-processing-final-v2",
            "source_project": "atlas",
            "source_run_id": "master-v1",
        }
    )
    libraries = []
    for context in contexts:
        for replicate in (1, 2):
            for assay in ("atac", "h3k27ac"):
                library_id = f"atlas_{context}_{assay}_{replicate}"
                bam = tmp_path / "bams" / f"{library_id}.bam"
                bai = tmp_path / "bams" / f"{library_id}.bam.bai"
                bam.parent.mkdir(exist_ok=True)
                bam.write_bytes(b"bam")
                bai.write_bytes(b"bai")
                libraries.append(
                    {
                        "id": library_id,
                        "assay": assay,
                        "cohort": "atlas",
                        "context": context,
                        "layout": (
                            "single"
                            if assay == "h3k27ac" and replicate == 1
                            else "paired"
                        ),
                        "genome": "dm6",
                        "bam": str(bam),
                        "bai": str(bai),
                        "bam_sha256": "b" * 64,
                        "bai_sha256": "c" * 64,
                        "filtering_contract": "short-read-processing-final-v2",
                        "qc_status": "accepted",
                        **(
                            {"estimated_fragment_length_bp": 165}
                            if assay == "h3k27ac" and replicate == 1
                            else {}
                        ),
                    }
                )
    config["activity"] = {
        "schema_version": 6,
        "master": master,
        "contexts": list(contexts),
        "libraries": libraries,
        "atac_fragment_maximum": 150,
        "atac_browser_extension_bp": 150,
        "normalization": "background_tmm_10kb_autosomes_v1",
        "h3k27ac_signal": "summit_max3_background_tmm_v3",
        "h3k27ac_activity_percentile_threshold": 0.60,
    }
    config["report"] = {
        "schema_version": 1,
        "source_roots": [],
        "source_files": [
            {
                "path": str(tmp_path / "master" / "stats_json"),
                "sha256": "a" * 64,
                "kind": "master_metrics",
                "source_root": str(tmp_path / "master"),
            }
        ],
    }
    nearest_annotation = tmp_path / "reference" / "nearest.gtf.gz"
    nearest_annotation.parent.mkdir(exist_ok=True)
    nearest_annotation.write_bytes(b"annotation")
    config["nearest_tss_links"] = {
        "schema_version": 2,
        "promoter_annotation": str(nearest_annotation),
        "promoter_annotation_checksum": "sha256:" + "e" * 64,
        "chromosome_scope": "all_annotation_chromosomes_in_chrom_sizes",
        "promoter_width_bp": 1000,
        "enhancer_classes": [
            "proximal_enhancer_like",
            "distal_enhancer_like",
        ],
        "promoter_activity": "summit_within_500bp_catalog_activity_v3",
        "candidate_assignment": "nearest_annotated_tss_all_ties_v1",
    }
    return config


def test_result_namespace_accepts_same_semantic_configuration(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["provenance"] = {}
    config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
    resolved = tmp_path / "results" / "project" / "run" / "provenance" / (
        "resolved_config.json"
    )
    resolved.parent.mkdir(parents=True)
    resolved.write_text(json.dumps(config), encoding="utf-8")

    guard_result_namespace(config, resolved)


def test_output_stage_does_not_change_scientific_semantic_digest():
    config = copy.deepcopy(BASE_CONFIG)
    config["output_stage"] = "trimming"
    trimming = workflow_semantic_sha256(config)
    config["output_stage"] = "master"

    assert workflow_semantic_sha256(config) == trimming


@pytest.mark.parametrize(
    ("start", "stop"),
    [
        ("trimming", "alignment"),
        ("alignment", "qc"),
        ("qc", "master"),
        ("master", "quantification"),
        ("quantification", "catalog"),
        ("catalog", "links"),
        ("links", "report"),
        ("catalog", "report"),
        ("report", "report"),
        ("master", "footprints"),
    ],
)
def test_each_logical_boundary_can_continue_forward(start, stop):
    assert validate_stage_selection(start, stop) == stop


def test_logical_boundary_cannot_move_backward():
    with pytest.raises(AcquisitionError, match="Cannot stop"):
        validate_stage_selection("qc", "alignment")


def test_report_sources_do_not_change_scientific_semantic_digest():
    config = copy.deepcopy(BASE_CONFIG)
    original = workflow_semantic_sha256(config)
    config["report"] = {
        "schema_version": 1,
        "source_roots": ["upstream"],
        "source_files": [
            {
                "path": "upstream/qc/metrics.json",
                "sha256": "a" * 64,
                "kind": "qc_metrics",
                "source_root": "upstream",
            }
        ],
    }

    assert workflow_semantic_sha256(config) == original


def test_result_namespace_rejects_changed_semantic_configuration(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["provenance"] = {}
    config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
    resolved = tmp_path / "results" / "project" / "run" / "provenance" / (
        "resolved_config.json"
    )
    resolved.parent.mkdir(parents=True)
    existing = copy.deepcopy(config)
    existing["provenance"]["semantic_sha256"] = "0" * 64
    resolved.write_text(json.dumps(existing), encoding="utf-8")

    with pytest.raises(AcquisitionError, match="Use a new run_id"):
        guard_result_namespace(config, resolved)


def test_workflow_config_rejects_stale_semantic_digest():
    config = copy.deepcopy(BASE_CONFIG)
    config["provenance"] = {"semantic_sha256": "0" * 64}

    with pytest.raises(AcquisitionError, match="does not match"):
        validate_workflow_config(config)


@pytest.mark.parametrize("branch", ["atac_qpois", "atac_se", "h3k27ac"])
def test_workflow_branches_dry_run(tmp_path, branch):
    config = copy.deepcopy(BASE_CONFIG)
    treatment = config["samples"][0]
    if branch == "atac_se":
        treatment["layout"] = "single"
        treatment.pop("r2")
    elif branch == "h3k27ac":
        config["assay"] = "chip_histone"
        config.pop("atac_qpois")
        treatment["peak_caller"] = _h3k27ac_callpeak()

    output = _dry_run(tmp_path, config, branch)
    assert "export_final_bam_manifest" in output
    if branch in {"atac_qpois", "atac_se"}:
        assert "prepare_atac_tn5_insertions" in output
        assert "call_atac_replicate_qpois" in output
        assert "refine_atac_replicate_qpois" in output
    if branch == "atac_se":
        assert "filter_atac_short_fragments" not in output
    if branch == "h3k27ac":
        assert "callpeak_h3k27ac" in output
        assert "chip_cross_correlation" in output


@pytest.mark.parametrize(
    "change, message",
    [
        ({"broad": False}, "must use broad mode"),
        ({"spmr": False}, "-B --SPMR"),
    ],
)
def test_h3k27ac_callpeak_contract_is_enforced(change, message):
    config = copy.deepcopy(BASE_CONFIG)
    config["assay"] = "chip_histone"
    config.pop("atac_qpois")
    config["samples"][0]["peak_caller"] = {**_h3k27ac_callpeak(), **change}

    with pytest.raises(AcquisitionError, match=message):
        validate_workflow_config(config)


def test_workflow_config_rejects_input_control_samples():
    config = copy.deepcopy(BASE_CONFIG)
    config["assay"] = "chip_histone"
    config.pop("atac_qpois")
    config["samples"][0]["peak_caller"] = _h3k27ac_callpeak()
    config["samples"][0]["control"] = "input_rep1"

    with pytest.raises(AcquisitionError, match="input controls are not supported"):
        validate_workflow_config(config)


def test_technical_lanes_align_separately_then_merge(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    sample = config["samples"][0]
    sample["accessions"].append("SRR123457")
    sample["r1"].append(sample["r1"][0])
    sample["r2"].append(sample["r2"][0])

    output = _dry_run(tmp_path, config, "two-lanes")

    assert re.search(r"align_lane\s+2", output)
    assert re.search(r"merge_and_mark_duplicates\s+1", output)


def test_restartable_alignment_splits_lanes_before_chunk_jobs(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["execution"] = {"alignment_chunk_pairs": 1}

    output = _dry_run(tmp_path, config, "restartable-alignment")

    assert re.search(r"split_lane_fastq\s+1", output)
    assert "--records-per-chunk 1" in output


def test_workflow_config_rejects_scaled_atac_qpois_signal():
    config = copy.deepcopy(BASE_CONFIG)
    config["samples"][0]["peak_caller"]["spmr"] = True
    with pytest.raises(AcquisitionError, match="invalid two-ended Tn5 qpois"):
        validate_workflow_config(config)


def test_workflow_config_rejects_removed_hmmratac_caller():
    config = copy.deepcopy(BASE_CONFIG)
    config["samples"][0]["peak_caller"] = {"command": "hmmratac"}
    with pytest.raises(AcquisitionError, match="peak caller must use callpeak"):
        validate_workflow_config(config)


@pytest.mark.parametrize(
    ("field", "value"),
    [("minimum_exponent", -1), ("maximum_exponent", 1), ("minimum_length", 0)],
)
def test_workflow_config_rejects_invalid_qpois_parameters(field, value):
    config = copy.deepcopy(BASE_CONFIG)
    config["atac_qpois"] = dict(ATAC_QPOIS_DEFAULTS)
    config["atac_qpois"][field] = value
    with pytest.raises(AcquisitionError, match="ATAC qpois parameters are invalid"):
        validate_workflow_config(config)


def _add_second_replicate(config: dict) -> None:
    second = copy.deepcopy(config["samples"][0])
    second["id"] = "atac_rep2"
    second["accessions"] = ["SRR123457"]
    config["samples"].append(second)


def _enable_consensus(config: dict) -> None:
    config["atac_consensus"] = {
        "enabled": True,
        "conditions": [
            {
                "id": "example",
                "label": "example",
                "samples": ["atac_rep1", "atac_rep2"],
            }
        ],
        "minimum_replicates": 2,
        "replicate_overlap_fraction": 0.5,
    }


def test_qpois_condition_consensus_builds_master_dhs_as_final_atac_step(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    _use_reviewed_final_bams(config, tmp_path)

    output = _dry_run(tmp_path, config, "qpois-consensus")

    for rule in (
        "pool_atac_condition_insertions",
        "call_atac_condition_qpois",
        "refine_atac_condition_qpois",
        "atac_condition_pileup_bigwig",
        "atac_condition_qpois_bigwig",
        "filter_atac_qpois_replicate_support",
        "build_atac_master_dhs",
    ):
        assert rule in output
    assert "clip_bedgraph.py" in output
    assert "replicate-supported.bed" in output
    assert "master_dhs.bed" in output


@pytest.mark.parametrize(
    ("stage", "required", "forbidden"),
    [
        (
            "trimming",
            ("trim_pe", "fastqc_raw", "fastqc_trimmed"),
            ("align_lane", "filter_bam", "call_atac_replicate_qpois"),
        ),
        (
            "alignment",
            (
                "align_lane",
                "filter_bam",
                "alignment_stats",
                "export_final_bam_manifest",
            ),
            ("call_atac_replicate_qpois", "frip", "build_atac_master_dhs"),
        ),
        (
            "qc",
            (
                "refine_atac_replicate_qpois",
                "frip",
                "multiqc",
                "library_qc_review_table",
            ),
            ("pool_atac_condition_insertions", "build_atac_master_dhs"),
        ),
    ],
)
def test_output_stage_prunes_downstream_rules(tmp_path, stage, required, forbidden):
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    config["output_stage"] = stage

    output = _dry_run(tmp_path, config, f"until-{stage}")

    for rule in required:
        assert rule in output
    for rule in forbidden:
        assert rule not in output


def test_workflow_config_rejects_negative_master_summit_distance():
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    config["atac_master"] = {
        "summit_max_distance": -1,
        "minimum_summit_separation": 50,
    }

    with pytest.raises(AcquisitionError, match="master DHS parameters"):
        validate_workflow_config(config)


def test_auto_reference_preparation_dry_run(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    reference_root = tmp_path / "references" / "dm6"
    config["reference"].update(
        {
            "fasta": str(reference_root / "dm6.fa"),
            "bowtie2_index": str(reference_root / "bowtie2" / "dm6"),
            "chrom_sizes": str(reference_root / "dm6.chrom.sizes"),
            "blacklist_bed": str(reference_root / "dm6.blacklist.bed"),
            "tss_bed": str(reference_root / "dm6.tss.bed"),
            "autosomes_file": str(reference_root / "dm6.autosomes.txt"),
            "preparation": {"mode": "download", **REFERENCE_SOURCES["dm6"]},
        }
    )
    output = _dry_run(tmp_path, config, "auto-reference")
    assert "prepare_reference_fasta" in output


def test_workflow_config_rejects_unchecked_reference_source():
    config = copy.deepcopy(BASE_CONFIG)
    config["reference"]["preparation"] = {
        "mode": "download",
        **REFERENCE_SOURCES["dm6"],
    }
    config["reference"]["preparation"]["fasta"] = {
        "url": "https://example.org/dm6.fa.gz",
        "checksum": "",
    }
    with pytest.raises(AcquisitionError, match="checksum"):
        validate_workflow_config(config)


def _external_final_bam(sample: str, tmp_path: Path) -> dict:
    bam = tmp_path / "external" / f"{sample}.bam"
    bai = tmp_path / "external" / f"{sample}.bam.bai"
    bam.parent.mkdir(exist_ok=True)
    bam.write_bytes(sample.encode())
    bai.write_bytes((sample + "-index").encode())
    return {
        "library_id": sample,
        "assay": "atac",
        "context": "example",
        "role": "treatment",
        "layout": "paired",
        "bam": str(bam),
        "bai": str(bai),
        "genome": "dm6",
        "filtering_contract": "short-read-processing-final-v2",
        "bam_sha256": "a" * 64,
        "bai_sha256": "b" * 64,
        "qc_status": "accepted",
    }


def _use_reviewed_final_bams(config: dict, tmp_path: Path) -> None:
    config["input_stage"] = "final-bam"
    config["output_stage"] = "master"
    for sample in config["samples"]:
        sample.pop("r1")
        sample.pop("r2")
        sample["final_bam"] = _external_final_bam(sample["id"], tmp_path)


def test_final_bam_dry_run_prunes_all_read_processing_and_alignment(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    _use_reviewed_final_bams(config, tmp_path)

    output = _dry_run(tmp_path, config, "external-final-bams")

    assert "validate_external_final_bam" in output
    assert "export_final_bam_manifest" in output
    assert "export_master_manifest" in output
    assert "prepare_atac_tn5_insertions" in output
    assert "build_atac_master_dhs" in output
    for forbidden in (
        "fastqc_raw",
        "trim_pe",
        "trim_se",
        "align_lane",
        "merge_and_mark_duplicates",
        "filter_bam",
        "bowtie2_index",
    ):
        assert forbidden not in output


def test_qc_to_master_reuses_lenient_replicate_peaks_without_recalling(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    _use_reviewed_final_bams(config, tmp_path)
    config["start_stage"] = "qc"
    for sample in config["samples"]:
        peak = tmp_path / "qc-peaks" / f"{sample['id']}.bed"
        peak.parent.mkdir(exist_ok=True)
        peak.write_text("chr1\t10\t100\tpeak\t10\t.\n")
        sample["qc_peak"] = {
            "path": str(peak),
            "sha256": "a" * 64,
            "method": "callpeak",
        }

    output = _dry_run(tmp_path, config, "qc-to-master")

    assert "validate_external_qc_peak" in output
    assert "call_atac_condition_qpois" in output
    assert "build_atac_master_dhs" in output
    assert "call_atac_replicate_qpois" not in output
    assert "refine_atac_replicate_qpois" not in output


def test_master_config_rejects_pending_final_bam():
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    config["input_stage"] = "final-bam"
    config["output_stage"] = "master"
    for sample in config["samples"]:
        sample.pop("r1")
        sample.pop("r2")
        sample["final_bam"] = {
            "bam": f"{sample['id']}.bam",
            "bai": f"{sample['id']}.bam.bai",
            "genome": "dm6",
            "filtering_contract": "short-read-processing-final-v2",
            "bam_sha256": "a" * 64,
            "bai_sha256": "b" * 64,
            "qc_status": "pending_review",
        }

    with pytest.raises(AcquisitionError, match="invalid for output stage 'master'"):
        validate_workflow_config(config)


def test_master_reuse_dry_run_validates_but_never_reconstructs(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["input_stage"] = "master"
    config["samples"] = []
    config.pop("atac_qpois")
    artifacts = {}
    for field in (
        "master_bed",
        "summits_bed",
        "membership_tsv",
        "context_matrix_tsv",
        "stats_json",
    ):
        path = tmp_path / "master" / field
        path.parent.mkdir(exist_ok=True)
        path.write_text(field)
        artifacts[field] = str(path)
        artifacts[f"{field}_sha256"] = "a" * 64
    config["external_master"] = {
        "genome": "dm6",
        "method": "reciprocal_summit_complete_linkage_v2",
        "input_filtering_contract": "short-read-processing-final-v2",
        "source_project": "atlas",
        "source_run_id": "master-v1",
        **artifacts,
    }

    output = _dry_run(tmp_path, config, "external-master")

    assert "validate_external_master" in output
    assert "export_master_manifest" in output
    assert "resolved_config_provenance" in output
    for forbidden in (
        "validate_external_final_bam",
        "prepare_atac_tn5_insertions",
        "build_atac_master_dhs",
        "align_lane",
        "filter_bam",
    ):
        assert forbidden not in output


def test_quantification_and_catalog_dry_run_prune_read_processing(tmp_path):
    config = _activity_config(tmp_path)

    output = _dry_run(tmp_path, config, "catalog")

    for expected in (
        "validate_activity_bam",
        "validate_activity_master",
        "prepare_activity_atac_insertions",
        "prepare_activity_h3k27ac_fragments",
        "prepare_activity_h3k27ac_single_fragments",
        "count_activity_library",
        "build_activity_h3k27ac_windows",
        "count_activity_h3k27ac_windows",
        "build_activity_h3k27ac_signal",
        "build_regulatory_element_catalog",
        "build_catalog_bed_tracks",
        "build_context_mean_bigwig",
        "build_context_igv_session",
        "build_all_contexts_igv_session",
        "build_activity_background_bins",
        "count_activity_background_library",
        "build_activity_tmm_inputs",
        "calculate_activity_tmm_factors",
        "build_activity_tmm_table",
        "master_dhs_activity.tsv.gz",
        "master_elements_long.tsv.gz",
        "master_elements_wide.tsv.gz",
        "ctx.elements.tsv.gz",
        "ctx.elements.bed",
        "ctx.atac.mean.background_tmm.bw",
        "ctx.h3k27ac.mean.background_tmm.bw",
        "igv/ctx.xml",
        "all-contexts.igv.xml",
    ):
        assert expected in output
    assert "build_integrated_qc_report" not in output
    for forbidden in (
        "fastqc_raw",
        "trim_pe",
        "align_lane",
        "filter_bam",
        "build_atac_master_dhs",
        "call_atac_replicate_qpois",
    ):
        assert forbidden not in output

    invalid_method = copy.deepcopy(config)
    invalid_method["activity"]["normalization"] = "unknown"
    with pytest.raises(
        AcquisitionError,
        match="Unsupported activity normalization",
    ):
        validate_workflow_config(invalid_method)

    invalid_schema = copy.deepcopy(config)
    invalid_schema["activity"]["schema_version"] = 3
    with pytest.raises(AcquisitionError, match="activity schema version"):
        validate_workflow_config(invalid_schema)

    invalid_h3k27ac_signal = copy.deepcopy(config)
    invalid_h3k27ac_signal["activity"]["h3k27ac_signal"] = (
        "summit_max3_500bp_v1"
    )
    with pytest.raises(AcquisitionError, match="regulatory-element method"):
        validate_workflow_config(invalid_h3k27ac_signal)

    invalid_browser_extension = copy.deepcopy(config)
    invalid_browser_extension["activity"]["atac_browser_extension_bp"] = 100
    with pytest.raises(AcquisitionError, match="browser extension must be 150 bp"):
        validate_workflow_config(invalid_browser_extension)

    insufficient_libraries = copy.deepcopy(config)
    insufficient_libraries["activity"]["libraries"] = [
        library
        for library in insufficient_libraries["activity"]["libraries"]
        if library["id"] != "atlas_ctx_h3k27ac_2"
    ]
    with pytest.raises(AcquisitionError, match="background TMM requires two h3k27ac"):
        validate_workflow_config(insufficient_libraries)

    invalid_report = copy.deepcopy(config)
    invalid_report["report"]["source_files"][0]["sha256"] = "invalid"
    with pytest.raises(AcquisitionError, match="source-file SHA-256"):
        validate_workflow_config(invalid_report)

    missing_nearest = copy.deepcopy(config)
    missing_nearest.pop("nearest_tss_links")
    with pytest.raises(AcquisitionError, match="require nearest_tss_links"):
        validate_workflow_config(missing_nearest)

    invalid_nearest_method = copy.deepcopy(config)
    invalid_nearest_method["nearest_tss_links"]["candidate_assignment"] = (
        "nearest_active_tss_v0"
    )
    with pytest.raises(AcquisitionError, match="Unsupported nearest-TSS"):
        validate_workflow_config(invalid_nearest_method)

    config["output_stage"] = "report"
    report_output = _dry_run(tmp_path, config, "report")
    assert "build_nearest_tss_promoters" in report_output
    assert "build_nearest_tss_enhancer_candidates" in report_output
    assert "enhancer_nearest_tss_candidates_wide.tsv.gz" in report_output
    assert "build_integrated_qc_report" in report_output
    assert "integrated_qc_report.pdf" in report_output
    config["output_stage"] = "quantification"
    quantification_output = _dry_run(tmp_path, config, "quantification-only")
    assert "build_activity_tmm_table" in quantification_output
    assert "build_activity_h3k27ac_windows" in quantification_output
    assert "count_activity_h3k27ac_windows" in quantification_output
    assert "build_activity_h3k27ac_signal" in quantification_output
    assert "master_dhs_activity.tsv.gz" in quantification_output
    assert "h3k27ac/master_dhs_signal.tsv.gz" in quantification_output
    assert "build_regulatory_element_catalog" not in quantification_output


def test_activity_config_accepts_configured_h3k27ac_activity_threshold(tmp_path):
    config = _activity_config(tmp_path)
    config["activity"]["h3k27ac_activity_percentile_threshold"] = 0.7
    config.setdefault("provenance", {})["semantic_sha256"] = (
        workflow_semantic_sha256(config)
    )

    validate_workflow_config(config)


@pytest.mark.parametrize("threshold", [0, 1, 1.5, -0.1, True, "0.6"])
def test_activity_config_rejects_invalid_h3k27ac_activity_threshold(
    tmp_path, threshold
):
    config = _activity_config(tmp_path)
    config["activity"]["h3k27ac_activity_percentile_threshold"] = threshold

    with pytest.raises(AcquisitionError, match=r"must be in \(0, 1\)"):
        validate_workflow_config(config)


def test_catalog_import_links_dry_run_uses_only_nearest_tss_rules(tmp_path):
    config = _activity_config(tmp_path, contexts=("ctx_a", "ctx_b"))
    config["output_stage"] = "links"
    imported_paths = {}
    for field in ("catalog", "metrics", "provenance", "resolved_config", "manifest"):
        path = tmp_path / f"imported-{field}"
        path.write_bytes(field.encode())
        imported_paths[field] = path
    config.pop("activity")
    config["input_stage"] = "catalog"
    config["start_stage"] = "catalog"
    config["run_id"] = "catalog-import"
    config["catalog_import"] = {
        "schema_version": 1,
        "manifest": str(imported_paths["manifest"]),
        "manifest_sha256": sha256_file(imported_paths["manifest"]),
        "genome": "dm6",
        "method": "catalog-v1",
        "contexts": ["ctx_a", "ctx_b"],
        "source_project": "atlas",
        "source_run_id": "catalog-v1",
        "source_semantic_sha256": "a" * 64,
        "source_sample_sheet_sha256": "b" * 64,
        **{
            key: value
            for field in ("catalog", "metrics", "provenance", "resolved_config")
            for key, value in (
                (field, str(imported_paths[field])),
                (f"{field}_sha256", sha256_file(imported_paths[field])),
            )
        },
    }

    imported_output = _dry_run(tmp_path, config, "links-import")

    assert "build_nearest_tss_promoters" in imported_output
    assert "build_nearest_tss_enhancer_candidates" in imported_output
    assert str(imported_paths["catalog"]) in imported_output
    assert "build_regulatory_element_catalog" not in imported_output
    assert "build_catalog_bed_tracks" not in imported_output


def test_workflow_config_rejects_removed_contact_branch(tmp_path):
    config = _activity_config(tmp_path)
    config["contacts"] = {"schema_version": 2}

    with pytest.raises(AcquisitionError, match="contact-link branch"):
        validate_workflow_config(config)


def _footprint_config(tmp_path: Path) -> dict:
    config = copy.deepcopy(BASE_CONFIG)
    config["input_stage"] = "master"
    config["start_stage"] = "master"
    config["output_stage"] = "footprints"
    config["samples"] = []
    config.pop("atac_qpois")
    artifacts = {}
    for field in ("master_bed", "summits_bed", "membership_tsv", "context_matrix_tsv", "stats_json"):
        path = tmp_path / "master" / field
        path.parent.mkdir(exist_ok=True)
        path.write_text(field)
        artifacts[field] = str(path)
        artifacts[f"{field}_sha256"] = "a" * 64
    config["external_master"] = {
        "genome": "dm6",
        "method": "reciprocal_summit_complete_linkage_v2",
        "input_filtering_contract": "short-read-processing-final-v2",
        "source_project": "atlas",
        "source_run_id": "master-v1",
        **artifacts,
    }
    libraries = []
    for library_id, context in (("atac_rep1", "eye"), ("atac_rep2", "eye"), ("atac_rep3", "wing")):
        bam = tmp_path / "bams" / f"{library_id}.bam"
        bam.parent.mkdir(exist_ok=True)
        bam.write_text(library_id)
        Path(f"{bam}.bai").write_text(library_id)
        libraries.append({
            "id": library_id, "context": context, "layout": "paired", "genome": "dm6",
            "bam": str(bam), "bai": f"{bam}.bai", "bam_sha256": "b" * 64, "bai_sha256": "c" * 64,
            "filtering_contract": "short-read-processing-final-v2", "qc_status": "accepted",
        })
    config["footprints"] = {
        "schema_version": 1,
        "method": "tobias_0.17.3_context_merged_atacorrect_scorebigwig_v1",
        "contexts": ["eye", "wing"],
        "libraries": libraries,
        "atacorrect": {"read_shift": [4, -5], "k_flank": 12, "bg_shift": 100, "window": 100,
                       "score_mat": "DWM"},
        "scorebigwig": {"score": "footprint", "fp_min": 20, "fp_max": 50, "flank_min": 10,
                        "flank_max": 30, "extend": 100},
    }
    return config


def test_footprint_dry_run_pools_contexts_and_skips_processing(tmp_path):
    output = _dry_run(tmp_path, _footprint_config(tmp_path), "footprints")

    for rule in ("validate_external_master", "validate_footprint_bam", "merge_footprint_context_bams",
                 "tobias_atacorrect", "tobias_scorebigwig", "export_stage_checkpoint"):
        assert rule in output
    assert "TOBIAS ATACorrect" in output and "--read_shift 4 -5" in output
    assert "TOBIAS ScoreBigwig" in output
    assert output.count("rule tobias_atacorrect:") == 2
    assert re.search(r"footprints/eye\.corrected\.bw", output)
    for forbidden in ("align_lane", "build_atac_master_dhs", "export_master_manifest",
                      "prepare_activity_atac_units"):
        assert forbidden not in output


def test_footprints_require_footprint_output_stage():
    config = _footprint_config(Path("/tmp"))
    config["output_stage"] = "master"
    with pytest.raises(AcquisitionError, match="footprints mapping requires output_stage=footprints"):
        validate_workflow_config(config)


def test_footprints_reject_unaccepted_libraries():
    config = _footprint_config(Path("/tmp"))
    config["footprints"]["libraries"][0]["qc_status"] = "rejected"
    with pytest.raises(AcquisitionError, match="accepted QC is required"):
        validate_workflow_config(config)


def test_footprint_output_requires_footprints_mapping():
    config = _footprint_config(Path("/tmp"))
    config.pop("footprints")
    with pytest.raises(AcquisitionError, match="requires a footprints mapping"):
        validate_workflow_config(config)
