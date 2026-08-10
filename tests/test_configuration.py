import json
from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError, FilePlan, RunPlan
from short_read_processing.artifacts import (
    CATALOG_FILE_FIELDS,
    FINAL_BAM_FILTERING_CONTRACT,
    sha256_file,
)
from short_read_processing.configuration import (
    ATAC_QPOIS_DEFAULTS,
    generate_activity_config,
    generate_catalog_links_config,
    generate_configs,
    generate_resume_config,
)
from short_read_processing.contact_metadata import DM6_ATLAS_CONTEXT_IDS
from short_read_processing.manifest import write_manifest
from short_read_processing.workflow_config import (
    resolve_input_paths,
    validate_workflow_config,
    workflow_semantic_sha256,
)


HEADER = "accession\tlibrary_id\tassay\tcontext\trole\tcontrol_library\tpeak_caller"


def _run_plan(root: Path, requested: str, run: str, layout: str = "PAIRED") -> RunPlan:
    run_dir = root / run
    run_dir.mkdir(parents=True)
    r1 = run_dir / f"{run}_1.fastq.gz"
    r1.write_bytes(b"r1")
    files = [FilePlan("https://example/r1", "", 2, r1, "r1")]
    if layout == "PAIRED":
        r2 = run_dir / f"{run}_2.fastq.gz"
        r2.write_bytes(b"r2")
        files.append(FilePlan("https://example/r2", "", 2, r2, "r2"))
    return RunPlan(
        requested_accession=requested,
        experiment_accession=requested if requested.startswith(("SRX", "ERX")) else "SRX999999",
        run_accession=run,
        library_layout=layout,
        backend="ena",
        run_dir=run_dir,
        files=files,
        status="downloaded",
    )


def _generate(tmp_path: Path, plans: list[RunPlan], sheet_text: str, **kwargs):
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, plans)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(sheet_text)
    return generate_configs(
        manifest_path=manifest,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="baseline",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        **kwargs,
    )


def test_atac_defaults_and_contexts_group_biological_and_technical_runs(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", accession, accession)
        for accession in ("SRR123456", "SRR123457", "SRR123458")
    ]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123458\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    config = yaml.safe_load(output.read_text())
    sample = config["samples"][0]

    assert config["assay"] == "atac"
    assert config["output_stage"] == "qc"
    assert config["reference"]["name"] == "dm6"
    assert sample["accessions"] == ["SRR123456", "SRR123457"]
    assert sample["context"] == "embryo"
    assert sample["layout"] == "paired"
    assert len(sample["r1"]) == 2
    assert len(sample["r2"]) == 2
    assert sample["peak_caller"] == {
        "command": "callpeak",
        "mode": "tn5_qpois",
        "format": "BED",
        "qvalue": 0.1,
        "broad": False,
        "nomodel": True,
        "shift": -75,
        "extsize": 150,
        "write_bedgraph": True,
        "spmr": False,
    }
    assert config["atac_consensus"] == {
        "enabled": True,
        "conditions": [
            {
                "id": "embryo",
                "label": "embryo",
                "samples": ["atac_rep1", "atac_rep2"],
            }
        ],
        "minimum_replicates": 2,
        "replicate_overlap_fraction": 0.5,
    }
    assert config["atac_master"] == {
        "summit_max_distance": 150,
        "minimum_summit_separation": 50,
    }
    assert sample["parameters"]["trimming"]["adapter_preset"] == "nextera"
    preparation = config["reference"]["preparation"]
    assert preparation["mode"] == "download"
    assert preparation["fasta"]["checksum"].startswith("md5:")
    assert preparation["annotation"]["url"].endswith("dm6.ncbiRefSeq.gtf.gz")


def test_identical_config_generation_does_not_replace_file(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [
        _run_plan(tmp_path / "raw", accession, accession)
        for accession in accessions
    ]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    original = output.read_bytes()
    output.touch()
    timestamp = output.stat().st_mtime_ns

    regenerated = generate_configs(
        manifest_path=tmp_path / "manifest.tsv",
        sample_sheet_path=tmp_path / "samples.tsv",
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="baseline",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
    )[0]

    assert regenerated.read_bytes() == original
    assert regenerated.stat().st_mtime_ns == timestamp


def test_config_records_requested_output_stage(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )

    output = _generate(tmp_path, plans, sheet, output_stage="qc")[0]

    assert yaml.safe_load(output.read_text())["output_stage"] == "qc"


def test_alignment_chunking_is_recorded_but_not_scientific(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet, alignment_chunk_pairs=2)[0]
    config = yaml.safe_load(output.read_text())
    without_execution = dict(config)
    without_execution.pop("execution")

    assert config["execution"] == {"alignment_chunk_pairs": 2}
    assert workflow_semantic_sha256(config) == workflow_semantic_sha256(
        without_execution
    )
    validate_workflow_config(config)


@pytest.mark.parametrize("chunk_pairs", [0, -1])
def test_alignment_chunking_rejects_nonpositive_sizes(tmp_path, chunk_pairs):
    plan = _run_plan(tmp_path / "raw", "SRR123456", "SRR123456")
    sheet = HEADER + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t\n"

    with pytest.raises(AcquisitionError, match="chunk pairs must be positive"):
        _generate(tmp_path, [plan], sheet, alignment_chunk_pairs=chunk_pairs)


def test_accession_input_cannot_construct_master_before_manual_review(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )

    with pytest.raises(AcquisitionError, match="Cannot stop at 'master'"):
        _generate(tmp_path, plans, sheet, output_stage="master")


def test_final_bam_input_cannot_stop_at_trimming(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )

    with pytest.raises(AcquisitionError, match="Cannot stop at 'trimming'"):
        generate_configs(
            manifest_path=None,
            sample_sheet_path=sheet,
            output_dir=tmp_path / "configs",
            project="test-project",
            run_id="reuse-bams",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
            require_fastq_files=True,
            input_stage="final-bam",
            output_stage="trimming",
            final_bam_manifest_path=tmp_path / "unused.tsv",
        )


def test_atac_callpeak_override_applies_to_every_context_replicate(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    header = HEADER + "\tmacs3_qvalue\tmacs3_shift\tmacs3_extsize"
    sheet = (
        header
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\tcallpeak\t0.01\t-75\t150"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\tcallpeak\t0.01\t-75\t150\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    peaks = [sample["peak_caller"] for sample in yaml.safe_load(output.read_text())["samples"]]

    assert all(peak["format"] == "BED" for peak in peaks)
    assert all(peak["qvalue"] == 0.01 for peak in peaks)
    assert all(peak["shift"] == -75 and peak["extsize"] == 150 for peak in peaks)
    assert all(peak["spmr"] is False for peak in peaks)


def test_h3k27ac_alias_is_broad_and_resolves_matched_input(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR100001", "SRR100001"),
        _run_plan(tmp_path / "raw", "SRR100002", "SRR100002"),
    ]
    sheet = (
        HEADER
        + "\nSRR100001\th3_rep1\th3k27ac\teye\ttreatment\teye_input\t"
        + "\nSRR100002\teye_input\th3k27ac\teye\tcontrol\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet, genome="hg38")[0]
    config = yaml.safe_load(output.read_text())
    treatment = config["samples"][0]

    assert config["assay"] == "chip_histone"
    assert config["reference"]["name"] == "hg38"
    assert treatment["control"] == "eye_input"
    assert treatment["peak_caller"]["command"] == "callpeak"
    assert treatment["peak_caller"]["broad"] is True
    assert treatment["peak_caller"]["broad_cutoff"] == 0.1
    assert "peak_caller" not in config["samples"][1]


def test_matched_ip_and_input_must_have_same_resolved_layout(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR100001", "SRR100001", layout="PAIRED"),
        _run_plan(tmp_path / "raw", "SRR100002", "SRR100002", layout="SINGLE"),
    ]
    sheet = (
        HEADER
        + "\nSRR100001\th3_rep1\th3k27ac\teye\ttreatment\teye_input\t"
        + "\nSRR100002\teye_input\th3k27ac\teye\tcontrol\t\t\n"
    )

    with pytest.raises(AcquisitionError, match="same read layout"):
        _generate(tmp_path, plans, sheet)


def test_mixed_atlas_table_generates_separate_atac_and_h3k27ac_configs(tmp_path):
    accessions = ("SRR100001", "SRR100002", "SRR100003", "SRR100004")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    sheet = (
        HEADER
        + "\nSRR100001\teye_atac_rep1\tatac\teye\ttreatment\t\t"
        + "\nSRR100002\teye_atac_rep2\tatac\teye\ttreatment\t\t"
        + "\nSRR100003\teye_h3_rep1\th3k27ac\teye\ttreatment\teye_input\t"
        + "\nSRR100004\teye_input\th3k27ac\teye\tcontrol\t\t\n"
    )

    outputs = _generate(tmp_path, plans, sheet)
    configs = {yaml.safe_load(path.read_text())["assay"]: path for path in outputs}

    assert set(configs) == {"atac", "chip_histone"}
    assert configs["atac"].name == "test-project.atac.dm6.yaml"
    assert configs["chip_histone"].name == "test-project.chip_histone.dm6.yaml"
    chip = yaml.safe_load(configs["chip_histone"].read_text())
    assert chip["samples"][0]["control"] == "eye_input"


def test_h3k27ac_can_run_ip_only_from_four_columns(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001")
    sheet = (
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\th3_rep1\th3k27ac\teye\n"
    )
    output = _generate(tmp_path, [plan], sheet)[0]
    treatment = yaml.safe_load(output.read_text())["samples"][0]

    assert "control" not in treatment
    assert treatment["peak_caller"]["format"] == "BAMPE"
    assert treatment["peak_caller"]["broad"] is True
    assert treatment["parameters"]["trimming"]["adapter_preset"] == "truseq"


def test_atac_context_requires_two_biological_libraries(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001")
    sheet = (
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    with pytest.raises(AcquisitionError, match="below minimum_replicates"):
        _generate(tmp_path, [plan], sheet)


def test_removed_hmmratac_peak_caller_is_rejected(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001")
    sheet = (
        HEADER
        + "\nSRR100001\tatac_rep1\tatac\teye\ttreatment\t\thmmratac\n"
    )
    with pytest.raises(AcquisitionError, match="peak_caller must be one of callpeak"):
        _generate(tmp_path, [plan], sheet)


def _write_final_bam_manifest(
    tmp_path: Path,
    libraries: list[tuple[str, str]],
    *,
    row_overrides: dict[str, dict[str, str]] | None = None,
) -> Path:
    manifest = tmp_path / "final-bams.tsv"
    rows = [
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\tnotes"
    ]
    row_overrides = row_overrides or {}
    for library_id, context in libraries:
        bam = tmp_path / "external" / f"{library_id}.bam"
        bai = tmp_path / "external" / f"{library_id}.bam.bai"
        bam.parent.mkdir(exist_ok=True)
        bam.write_bytes(f"bam-{library_id}".encode())
        bai.write_bytes(f"bai-{library_id}".encode())
        override = row_overrides.get(library_id, {})
        rows.append(
            f"{library_id}\tatac\t{context}\ttreatment\tpaired\t{bam}\t{bai}\tdm6\t"
            f"{FINAL_BAM_FILTERING_CONTRACT}\t{sha256_file(bam)}\t"
            f"{sha256_file(bai)}\t{override.get('qc_status', 'accepted')}\t"
            f"{override.get('notes', '')}"
        )
    manifest.write_text("\n".join(rows) + "\n")
    return manifest


def _write_activity_bam_manifest(
    tmp_path: Path,
    name: str,
    libraries: list[tuple[str, str, str]],
    *,
    qc_status: str = "accepted",
    row_overrides: dict[str, dict[str, str]] | None = None,
) -> Path:
    manifest = tmp_path / f"{name}.final-bams.tsv"
    rows = [
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\t"
        "estimated_fragment_length_bp\tnotes"
    ]
    row_overrides = row_overrides or {}
    for library_id, assay, context in libraries:
        bam = tmp_path / "external" / f"{library_id}.bam"
        bai = tmp_path / "external" / f"{library_id}.bam.bai"
        bam.parent.mkdir(exist_ok=True)
        bam.write_bytes(f"bam-{library_id}".encode())
        bai.write_bytes(f"bai-{library_id}".encode())
        override = row_overrides.get(library_id, {})
        rows.append(
            f"{library_id}\t{assay}\t{context}\ttreatment\t"
            f"{override.get('layout', 'paired')}\t{bam}\t{bai}\tdm6\t"
            f"{FINAL_BAM_FILTERING_CONTRACT}\t{sha256_file(bam)}\t"
            f"{sha256_file(bai)}\t{override.get('qc_status', qc_status)}\t"
            f"{override.get('estimated_fragment_length_bp', '')}\t"
            f"{override.get('notes', '')}"
        )
    manifest.write_text("\n".join(rows) + "\n")
    return manifest


def _write_master_manifest(tmp_path: Path) -> Path:
    fields = (
        "master_bed",
        "summits_bed",
        "membership_tsv",
        "context_matrix_tsv",
        "stats_json",
    )
    artifacts = {}
    for field in fields:
        path = tmp_path / "master" / field
        path.parent.mkdir(exist_ok=True)
        path.write_text(f"{field}\n")
        artifacts[field] = path
    manifest = tmp_path / "master.tsv"
    columns = [
        "genome",
        "method",
        "input_filtering_contract",
        "source_project",
        "source_run_id",
    ]
    columns.extend(item for field in fields for item in (field, f"{field}_sha256"))
    values = [
        "dm6",
        "reciprocal_summit_complete_linkage_v2",
        FINAL_BAM_FILTERING_CONTRACT,
        "atlas",
        "master-v1",
    ]
    for field in fields:
        values.extend((str(artifacts[field]), sha256_file(artifacts[field])))
    manifest.write_text("\t".join(columns) + "\n" + "\t".join(values) + "\n")
    return manifest


def test_final_bam_mode_generates_reuse_only_samples(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path, [("atac_rep1", "eye"), ("atac_rep2", "eye")]
    )

    output = generate_configs(
        manifest_path=None,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="reuse-bams",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        input_stage="final-bam",
        final_bam_manifest_path=manifest,
    )[0]
    config = yaml.safe_load(output.read_text())

    assert config["input_stage"] == "final-bam"
    assert all("final_bam" in sample for sample in config["samples"])
    assert all("r1" not in sample and "r2" not in sample for sample in config["samples"])
    assert config["provenance"]["final_bam_manifest_sha256"] == sha256_file(manifest)
    assert len(config["provenance"]["semantic_sha256"]) == 64


def test_master_requires_every_atac_library_to_be_reviewed(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path,
        [("atac_rep1", "eye"), ("atac_rep2", "eye")],
        row_overrides={"atac_rep2": {"qc_status": "pending_review"}},
    )

    with pytest.raises(AcquisitionError, match="pending_review: atac_rep2"):
        generate_configs(
            manifest_path=None,
            sample_sheet_path=sheet,
            output_dir=tmp_path / "configs",
            project="test-project",
            run_id="master-v1",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
            require_fastq_files=True,
            input_stage="final-bam",
            output_stage="master",
            final_bam_manifest_path=manifest,
        )


def test_master_excludes_documented_rejection_and_records_provenance(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
        "SRR100003\tatac_rep3\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path,
        [("atac_rep1", "eye"), ("atac_rep2", "eye"), ("atac_rep3", "eye")],
        row_overrides={
            "atac_rep1": {
                "qc_status": "rejected",
                "notes": "failed manual ATAC QC",
            }
        },
    )

    output = generate_configs(
        manifest_path=None,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="master-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        input_stage="final-bam",
        output_stage="master",
        final_bam_manifest_path=manifest,
    )[0]
    config = yaml.safe_load(output.read_text())

    assert [sample["id"] for sample in config["samples"]] == [
        "atac_rep2",
        "atac_rep3",
    ]
    assert config["atac_consensus"]["conditions"][0]["samples"] == [
        "atac_rep2",
        "atac_rep3",
    ]
    assert config["provenance"]["excluded_master_libraries"] == [
        {
            "id": "atac_rep1",
            "assay": "atac",
            "context": "eye",
            "role": "treatment",
            "layout": "paired",
            "qc_status": "rejected",
            "reason": "failed manual ATAC QC",
        }
    ]


def test_pending_final_bams_remain_valid_for_qc_only(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path,
        [("atac_rep1", "eye"), ("atac_rep2", "eye")],
        row_overrides={
            "atac_rep1": {"qc_status": "pending_review"},
            "atac_rep2": {"qc_status": "pending_review"},
        },
    )

    output = generate_configs(
        manifest_path=None,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="qc-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        input_stage="final-bam",
        output_stage="qc",
        final_bam_manifest_path=manifest,
    )[0]

    assert yaml.safe_load(output.read_text())["output_stage"] == "qc"


def _write_qc_checkpoint(tmp_path: Path, source_config: dict) -> Path:
    artifacts = {}
    for sample in source_config["samples"]:
        peak = tmp_path / "qc-peaks" / f"{sample['id']}.bed"
        peak.parent.mkdir(exist_ok=True)
        peak.write_text("chr2L\t10\t80\tpeak\t10\t.\n")
        artifacts[f"replicate_peak.{sample['id']}"] = {
            "path": str(peak),
            "sha256": sha256_file(peak),
        }
    checkpoint = tmp_path / "qc.checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "stage": "qc",
                "source_project": "test-project",
                "source_run_id": "qc-v1",
                "semantic_sha256": source_config["provenance"]["semantic_sha256"],
                "parameters": {
                    "assay": "atac",
                    "atac_qpois": ATAC_QPOIS_DEFAULTS,
                    "peak_callers": {
                        sample["id"]: sample["peak_caller"]
                        for sample in source_config["samples"]
                    },
                },
                "artifacts": artifacts,
            },
            indent=2,
        )
        + "\n"
    )
    return checkpoint


def test_qc_checkpoint_reuses_lenient_replicate_peaks_for_master(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path, [("atac_rep1", "eye"), ("atac_rep2", "eye")]
    )
    qc_config_path = generate_configs(
        manifest_path=None,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "qc-configs",
        project="test-project",
        run_id="qc-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        input_stage="final-bam",
        start_stage="alignment",
        output_stage="qc",
        final_bam_manifest_path=manifest,
    )[0]
    qc_config = yaml.safe_load(qc_config_path.read_text())
    checkpoint = _write_qc_checkpoint(tmp_path, qc_config)

    master_path = generate_configs(
        manifest_path=None,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "master-configs",
        project="test-project",
        run_id="master-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        input_stage="final-bam",
        start_stage="qc",
        output_stage="master",
        final_bam_manifest_path=manifest,
        qc_checkpoint_manifest_path=checkpoint,
    )[0]
    master = yaml.safe_load(master_path.read_text())

    assert master["start_stage"] == "qc"
    assert master["atac_qpois"] == ATAC_QPOIS_DEFAULTS
    assert all(sample["peak_caller"]["qvalue"] == 0.1 for sample in master["samples"])
    assert all("qc_peak" in sample for sample in master["samples"])


def test_qc_checkpoint_rejects_changed_lenient_peak_parameters(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path, [("atac_rep1", "eye"), ("atac_rep2", "eye")]
    )
    source_path = generate_configs(
        manifest_path=None,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "qc-configs",
        project="test-project",
        run_id="qc-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
        input_stage="final-bam",
        start_stage="alignment",
        output_stage="qc",
        final_bam_manifest_path=manifest,
    )[0]
    source = yaml.safe_load(source_path.read_text())
    checkpoint = _write_qc_checkpoint(tmp_path, source)
    payload = json.loads(checkpoint.read_text())
    payload["parameters"]["peak_callers"]["atac_rep1"]["qvalue"] = 0.01
    checkpoint.write_text(json.dumps(payload) + "\n")

    with pytest.raises(AcquisitionError, match="peak parameters"):
        generate_configs(
            manifest_path=None,
            sample_sheet_path=sheet,
            output_dir=tmp_path / "master-configs",
            project="test-project",
            run_id="master-v1",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
            require_fastq_files=True,
            input_stage="final-bam",
            start_stage="qc",
            output_stage="master",
            final_bam_manifest_path=manifest,
            qc_checkpoint_manifest_path=checkpoint,
        )


def test_checkpoint_resume_restores_original_semantic_paths(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR100001", "SRR100001"),
        _run_plan(tmp_path / "raw", "SRR100002", "SRR100002"),
    ]
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    generated = _generate(
        tmp_path,
        plans,
        sheet.read_text(),
        output_stage="trimming",
    )[0]
    source = yaml.safe_load(generated.read_text())
    recorded_semantic = source["provenance"]["semantic_sha256"]
    resolve_input_paths(source, tmp_path)
    resolved = tmp_path / "results" / "resolved_config.json"
    resolved.parent.mkdir()
    resolved.write_text(json.dumps(source) + "\n")
    trimmed = tmp_path / "work" / "trimmed.fastq.gz"
    trimmed.parent.mkdir()
    trimmed.write_bytes(b"trimmed")
    checkpoint = tmp_path / "trimming.checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "stage": "trimming",
                "source_project": "test-project",
                "source_run_id": "baseline",
                "semantic_sha256": recorded_semantic,
                "parameters": {},
                "artifacts": {
                    "resolved_config": {
                        "path": str(resolved),
                        "sha256": sha256_file(resolved),
                    },
                    "trimmed_fastq.atac_rep1": {
                        "path": str(trimmed),
                        "sha256": sha256_file(trimmed),
                    },
                },
            }
        )
        + "\n"
    )

    resumed_path = generate_resume_config(
        checkpoint_manifest_path=checkpoint,
        start_stage="trimming",
        output_stage="alignment",
        sample_sheet_path=sheet,
        output_dir=tmp_path / "resume-configs",
        path_base=tmp_path,
    )
    resumed = yaml.safe_load(resumed_path.read_text())

    assert resumed["start_stage"] == "trimming"
    assert resumed["output_stage"] == "alignment"
    assert workflow_semantic_sha256(resumed) == recorded_semantic


def test_final_bam_mode_rejects_partial_selected_assay(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(tmp_path, [("atac_rep1", "eye")])

    with pytest.raises(AcquisitionError, match="incomplete"):
        generate_configs(
            manifest_path=None,
            sample_sheet_path=sheet,
            output_dir=tmp_path / "configs",
            project="test-project",
            run_id="reuse-bams",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
            require_fastq_files=True,
            input_stage="final-bam",
            final_bam_manifest_path=manifest,
        )


def test_identical_final_bam_config_keeps_content_and_mtime(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
        "SRR100002\tatac_rep2\tatac\teye\n"
    )
    manifest = _write_final_bam_manifest(
        tmp_path, [("atac_rep1", "eye"), ("atac_rep2", "eye")]
    )
    arguments = {
        "manifest_path": None,
        "sample_sheet_path": sheet,
        "output_dir": tmp_path / "configs",
        "project": "test-project",
        "run_id": "reuse-bams",
        "reference_root": tmp_path / "references",
        "path_base": tmp_path,
        "require_fastq_files": True,
        "input_stage": "final-bam",
        "final_bam_manifest_path": manifest,
    }
    output = generate_configs(**arguments)[0]
    original = output.read_bytes()
    output.touch()
    timestamp = output.stat().st_mtime_ns

    regenerated = generate_configs(**arguments)[0]

    assert regenerated.read_bytes() == original
    assert regenerated.stat().st_mtime_ns == timestamp


def test_activity_config_requires_complete_accepted_contexts(tmp_path):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatlas_atac\tatac\teye\n"
        "SRR100002\tatlas_h3\th3k27ac\teye\n"
        "SRR100003\tatlas_atac_2\tatac\teye\n"
        "SRR100004\tatlas_h3_2\th3k27ac\teye\n"
    )
    manifest = _write_activity_bam_manifest(
        tmp_path,
        "atlas",
        [
            ("atlas_atac", "atac", "eye"),
            ("atlas_h3", "chip_histone", "eye"),
            ("atlas_atac_2", "atac", "eye"),
            ("atlas_h3_2", "chip_histone", "eye"),
        ],
        row_overrides={
            "atlas_h3": {
                "layout": "single",
                "estimated_fragment_length_bp": "165",
            }
        },
    )
    master_manifest = _write_master_manifest(tmp_path)

    output = generate_activity_config(
        sample_sheet_path=sheet,
        final_bam_manifests=[manifest],
        master_manifest_path=master_manifest,
        output_dir=tmp_path / "configs",
        project="activity-test",
        run_id="activity-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_files=True,
        output_stage="quantification",
    )
    config = yaml.safe_load(output.read_text())

    assert config["assay"] == "activity"
    assert config["input_stage"] == "quantification"
    assert config["output_stage"] == "quantification"
    assert config["report"]["schema_version"] == 1
    assert {
        source["kind"] for source in config["report"]["source_files"]
    } >= {
        "sample_sheet",
        "sample_sheet_schema",
        "master_manifest",
        "final_bam_manifest",
    }
    assert config["samples"] == []
    assert config["activity"]["schema_version"] == 4
    assert (
        config["activity"]["h3k27ac_signal"]
        == "summit_windows_background_tmm_v2"
    )
    assert config["activity"]["contexts"] == ["eye"]
    assert config["activity"]["atac_browser_extension_bp"] == 150
    assert config["activity"]["normalization"] == "background_tmm_10kb_autosomes_v1"
    assert {
        (library["cohort"], library["assay"])
        for library in config["activity"]["libraries"]
    } == {
        ("atlas", "atac"),
        ("atlas", "h3k27ac"),
    }
    assert all(
        library["qc_status"] == "accepted"
        for library in config["activity"]["libraries"]
    )
    atlas_h3 = next(
        library
        for library in config["activity"]["libraries"]
        if library["id"] == "atlas_h3"
    )
    assert atlas_h3["layout"] == "single"
    assert atlas_h3["estimated_fragment_length_bp"] == 165
    assert len(config["provenance"]["semantic_sha256"]) == 64


def test_activity_config_requires_fragment_length_for_single_end_h3k27ac(
    tmp_path,
):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatlas_atac\tatac\teye\n"
        "SRR100002\tatlas_h3\th3k27ac\teye\n"
        "SRR100003\tatlas_atac_2\tatac\teye\n"
        "SRR100004\tatlas_h3_2\th3k27ac\teye\n"
    )
    manifest = _write_activity_bam_manifest(
        tmp_path,
        "atlas-single",
        [
            ("atlas_atac", "atac", "eye"),
            ("atlas_h3", "chip_histone", "eye"),
            ("atlas_atac_2", "atac", "eye"),
            ("atlas_h3_2", "chip_histone", "eye"),
        ],
        row_overrides={"atlas_h3": {"layout": "single"}},
    )

    with pytest.raises(AcquisitionError, match="estimated_fragment_length_bp"):
        generate_activity_config(
            sample_sheet_path=sheet,
            final_bam_manifests=[manifest],
            master_manifest_path=_write_master_manifest(tmp_path),
            output_dir=tmp_path / "configs",
            project="activity-test",
            run_id="activity-v1",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
            require_files=True,
        )


def test_activity_config_rejects_pending_review_library(tmp_path):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_1\tatac\teye\n"
        "SRR100002\tatac_2\tatac\teye\n"
        "SRR100003\th3_1\th3k27ac\teye\n"
        "SRR100004\th3_2\th3k27ac\teye\n"
    )
    manifest = _write_activity_bam_manifest(
        tmp_path,
        "pending",
        [
            ("atac_1", "atac", "eye"),
            ("atac_2", "atac", "eye"),
            ("h3_1", "chip_histone", "eye"),
            ("h3_2", "chip_histone", "eye"),
        ],
        qc_status="pending_review",
    )

    with pytest.raises(AcquisitionError, match="requires qc_status='accepted'"):
        generate_activity_config(
            sample_sheet_path=sheet,
            final_bam_manifests=[manifest],
            master_manifest_path=_write_master_manifest(tmp_path),
            output_dir=tmp_path / "configs",
            project="activity-test",
            run_id="activity-v1",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
            require_files=True,
        )


def test_activity_config_records_and_skips_documented_rejection(tmp_path):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatlas_atac\tatac\teye\n"
        "SRR100002\tatlas_h3_good\th3k27ac\teye\n"
        "SRR100003\tatlas_h3_bad\th3k27ac\teye\n"
        "SRR100004\tatlas_atac_2\tatac\teye\n"
        "SRR100005\tatlas_h3_good_2\th3k27ac\teye\n"
    )
    manifest = _write_activity_bam_manifest(
        tmp_path,
        "atlas-reviewed",
        [
            ("atlas_atac", "atac", "eye"),
            ("atlas_h3_good", "chip_histone", "eye"),
            ("atlas_h3_bad", "chip_histone", "eye"),
            ("atlas_atac_2", "atac", "eye"),
            ("atlas_h3_good_2", "chip_histone", "eye"),
        ],
        row_overrides={
            "atlas_h3_bad": {
                "qc_status": "rejected",
                "notes": "insufficient depth",
            }
        },
    )
    output = generate_activity_config(
        sample_sheet_path=sheet,
        final_bam_manifests=[manifest],
        master_manifest_path=_write_master_manifest(tmp_path),
        output_dir=tmp_path / "configs",
        project="activity-test",
        run_id="activity-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_files=True,
    )
    config = yaml.safe_load(output.read_text())

    assert "atlas_h3_bad" not in {
        library["id"] for library in config["activity"]["libraries"]
    }
    assert config["provenance"]["excluded_activity_libraries"] == [
        {
            "id": "atlas_h3_bad",
            "assay": "h3k27ac",
            "context": "eye",
            "cohort": "atlas",
            "layout": "paired",
            "qc_status": "rejected",
            "estimated_fragment_length_bp": None,
            "reason": "insufficient depth",
        }
    ]


def test_complete_dm6_atlas_config_defaults_to_nearest_tss_links(tmp_path):
    sheet_rows = ["accession\tlibrary_id\tassay\tcontext"]
    manifest_rows = []
    accession = 200000
    for context in sorted(DM6_ATLAS_CONTEXT_IDS):
        for assay, manifest_assay in (
            ("atac", "atac"),
            ("h3k27ac", "chip_histone"),
        ):
            library_id = f"{context}_{assay}"
            sheet_rows.append(
                f"SRR{accession}\t{library_id}\t{assay}\t{context}"
            )
            manifest_rows.append((library_id, manifest_assay, context))
            accession += 1
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text("\n".join(sheet_rows) + "\n", encoding="utf-8")
    manifest = _write_activity_bam_manifest(tmp_path, "atlas", manifest_rows)

    output = generate_activity_config(
        sample_sheet_path=sheet,
        final_bam_manifests=[manifest],
        master_manifest_path=_write_master_manifest(tmp_path),
        output_dir=tmp_path / "configs",
        project="activity-test",
        run_id="links-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_files=True,
        output_stage="links",
    )
    config = yaml.safe_load(output.read_text())

    assert config["output_stage"] == "links"
    assert "contacts" not in config
    assert config["nearest_tss_links"]["candidate_assignment"] == (
        "nearest_annotated_tss_all_ties_v1"
    )

    output = generate_activity_config(
        sample_sheet_path=sheet,
        final_bam_manifests=[manifest],
        master_manifest_path=_write_master_manifest(tmp_path),
        output_dir=tmp_path / "configs",
        project="activity-test",
        run_id="links-with-contacts-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_files=True,
        output_stage="links",
        include_contacts=True,
    )
    config = yaml.safe_load(output.read_text())

    assert len(config["contacts"]["contexts"]) == 9
    assert {
        row["id"]
        for row in config["contacts"]["contexts"]
        if row["strategy"] == "powerlaw"
    } == {"e13", "hid"}
    assert config["contacts"]["source_manifest"].endswith(
        "resources/atlas_contact_sources.tsv"
    )


def _write_catalog_bundle_manifest(tmp_path: Path, sheet: Path) -> Path:
    root = tmp_path / "source-catalog"
    catalog_root = root / "activity" / "catalog"
    config_root = root / "provenance" / "configs"
    catalog_root.mkdir(parents=True)
    config_root.mkdir(parents=True)
    catalog = catalog_root / "master_elements_long.tsv.gz"
    catalog.write_bytes(b"catalog")
    catalog_digest = sha256_file(catalog)
    metrics = catalog_root / "regulatory_element_metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "catalog_sha256": catalog_digest,
                "context_count": len(DM6_ATLAS_CONTEXT_IDS),
                "method": "catalog-v1",
            }
        )
    )
    provenance = catalog_root / "regulatory_element_provenance.json"
    provenance.write_text(
        json.dumps({"outputs": {"catalog": {"sha256": catalog_digest}}})
    )
    resolved = config_root / "report.resolved_config.json"
    resolved.write_text(
        json.dumps(
            {
                "project": "source-atlas",
                "run_id": "catalog-v1",
                "reference": {"name": "dm6"},
                "activity": {"contexts": sorted(DM6_ATLAS_CONTEXT_IDS)},
                "provenance": {
                    "semantic_sha256": "c" * 64,
                    "sample_sheet_sha256": sha256_file(sheet),
                },
            }
        )
    )
    files = {
        "catalog": catalog,
        "metrics": metrics,
        "provenance": provenance,
        "resolved_config": resolved,
    }
    columns = [
        "genome",
        "method",
        "contexts",
        "source_project",
        "source_run_id",
    ]
    columns.extend(
        item for field in CATALOG_FILE_FIELDS for item in (field, f"{field}_sha256")
    )
    values = [
        "dm6",
        "catalog-v1",
        ",".join(sorted(DM6_ATLAS_CONTEXT_IDS)),
        "source-atlas",
        "catalog-v1",
    ]
    values.extend(
        item
        for field, path in files.items()
        for item in (str(path), sha256_file(path))
    )
    manifest = tmp_path / "catalog-manifest.tsv"
    manifest.write_text("\t".join(columns) + "\n" + "\t".join(values) + "\n")
    return manifest


def test_catalog_import_generates_a_new_links_only_namespace(tmp_path):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\tab\n"
    )
    manifest = _write_catalog_bundle_manifest(tmp_path, sheet)

    output = generate_catalog_links_config(
        sample_sheet_path=sheet,
        catalog_manifest_path=manifest,
        output_dir=tmp_path / "configs",
        project="atlas-links",
        run_id="from-catalog-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
    )
    config = yaml.safe_load(output.read_text())

    assert config["input_stage"] == "catalog"
    assert config["output_stage"] == "links"
    assert config["samples"] == []
    assert config["catalog_import"]["source_run_id"] == "catalog-v1"
    assert config["catalog_import"]["catalog"].endswith(
        "activity/catalog/master_elements_long.tsv.gz"
    )
    assert "contacts" not in config
    assert config["nearest_tss_links"]["enhancer_classes"] == [
        "proximal_enhancer_like",
        "distal_enhancer_like",
    ]
    validate_workflow_config(config)
    resolve_input_paths(config, tmp_path)


def test_catalog_import_can_opt_into_contact_branch(tmp_path):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\tab\n"
    )
    manifest = _write_catalog_bundle_manifest(tmp_path, sheet)

    output = generate_catalog_links_config(
        sample_sheet_path=sheet,
        catalog_manifest_path=manifest,
        output_dir=tmp_path / "configs",
        project="atlas-links",
        run_id="from-catalog-contacts-v1",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        include_contacts=True,
    )
    config = yaml.safe_load(output.read_text())

    assert {row["id"] for row in config["contacts"]["contexts"]} == set(
        DM6_ATLAS_CONTEXT_IDS
    )


def test_catalog_import_rejects_a_different_sample_sheet(tmp_path):
    sheet = tmp_path / "atlas.tsv"
    sheet.write_text(
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\tab\n"
    )
    manifest = _write_catalog_bundle_manifest(tmp_path, sheet)
    sheet.write_text(sheet.read_text() + "SRR100002\tatac_rep2\tatac\tab\n")

    with pytest.raises(AcquisitionError, match="source sample sheet"):
        generate_catalog_links_config(
            sample_sheet_path=sheet,
            catalog_manifest_path=manifest,
            output_dir=tmp_path / "configs",
            project="atlas-links",
            run_id="from-catalog-v1",
            reference_root=tmp_path / "references",
            path_base=tmp_path,
        )
