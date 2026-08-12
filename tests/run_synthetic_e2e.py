#!/usr/bin/env python3
"""Run a deterministic local FASTQ-to-report integration test.

The fixture uses two contexts with two paired-end ATAC and two paired-end
H3K27ac libraries per context.  Generated data and results live below the
ignored ``results/synthetic-e2e-runtime`` directory by default.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from short_read_processing.artifacts import (  # noqa: E402
    read_final_bam_manifest,
    read_master_manifest,
)
from short_read_processing.workflow_config import (  # noqa: E402
    workflow_semantic_sha256,
)


SEED = 20260802
READ_LENGTH = 50
CHROMOSOME_LENGTH = 200_000
MITOCHONDRIAL_LENGTH = 2_000
CONTEXT_TARGETS = {
    "alpha": (20_000, 50_000, 90_000),
    "beta": (110_000, 140_000, 170_000),
}
SHARED_TARGETS = (30_000, 70_000)
SAMPLE_ROWS = (
    ("SRR99000001", "alpha_atac_rep1", "atac", "alpha"),
    ("SRR99000002", "alpha_atac_rep2", "atac", "alpha"),
    ("SRR99000003", "beta_atac_rep1", "atac", "beta"),
    ("SRR99000004", "beta_atac_rep2", "atac", "beta"),
    ("SRR99000005", "alpha_h3k27ac_rep1", "h3k27ac", "alpha"),
    ("SRR99000006", "alpha_h3k27ac_rep2", "h3k27ac", "alpha"),
    ("SRR99000007", "beta_h3k27ac_rep1", "h3k27ac", "beta"),
    ("SRR99000008", "beta_h3k27ac_rep2", "h3k27ac", "beta"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def wrapped_fasta(name: str, sequence: str) -> str:
    return ">" + name + "\n" + "\n".join(
        sequence[index : index + 80]
        for index in range(0, len(sequence), 80)
    ) + "\n"


def write_fastq_pair(
    *,
    chromosome: str,
    context: str,
    assay: str,
    replicate: int,
    r1_path: Path,
    r2_path: Path,
) -> None:
    rng = random.Random(f"{SEED}:{context}:{assay}:{replicate}")
    targets = (*SHARED_TARGETS, *CONTEXT_TARGETS[context])
    fragments: list[tuple[int, int]] = []
    if assay == "atac":
        for center in targets:
            for _ in range(600):
                length = rng.randint(105, 145)
                start = center - length // 2 + rng.randint(-110, 110)
                fragments.append((start, length))
        for _ in range(2_000):
            length = rng.randint(105, 145)
            start = rng.randint(2_500, CHROMOSOME_LENGTH - length - 2_500)
            fragments.append((start, length))
    else:
        for center in targets:
            for _ in range(900):
                length = rng.randint(180, 280)
                start = center - length // 2 + rng.randint(-280, 280)
                fragments.append((start, length))
        for _ in range(4_000):
            length = rng.randint(180, 280)
            start = rng.randint(2_500, CHROMOSOME_LENGTH - length - 2_500)
            fragments.append((start, length))
    rng.shuffle(fragments)

    quality = "I" * READ_LENGTH
    r1_path.parent.mkdir(parents=True, exist_ok=True)
    with r1_path.open("w", encoding="ascii", newline="\n") as r1_handle, r2_path.open(
        "w", encoding="ascii", newline="\n"
    ) as r2_handle:
        for index, (start, length) in enumerate(fragments, start=1):
            end = start + length
            read1 = chromosome[start : start + READ_LENGTH]
            read2 = reverse_complement(chromosome[end - READ_LENGTH : end])
            name = f"{context}_{assay}_r{replicate}_{index}"
            r1_handle.write(f"@{name}/1\n{read1}\n+\n{quality}\n")
            r2_handle.write(f"@{name}/2\n{read2}\n+\n{quality}\n")


def generate_inputs(runtime: Path) -> tuple[Path, Path, Path, dict[Path, str]]:
    print("[setup] generating deterministic synthetic genome and FASTQs", flush=True)
    rng = random.Random(SEED)
    chromosome = "".join(rng.choices("ACGT", k=CHROMOSOME_LENGTH))
    mitochondrion = "".join(rng.choices("ACGT", k=MITOCHONDRIAL_LENGTH))

    reference_root = runtime / "references"
    dm6_root = reference_root / "dm6"
    dm6_root.mkdir(parents=True, exist_ok=True)
    (dm6_root / "dm6.fa").write_text(
        wrapped_fasta("chr1", chromosome) + wrapped_fasta("chrM", mitochondrion),
        encoding="ascii",
    )
    # This interval overlaps the shared 30-kb synthetic DHS by construction.
    (dm6_root / "dm6.blacklist.bed").write_text(
        "chr1\t29750\t30250\n", encoding="ascii"
    )
    (dm6_root / "dm6.tss.bed").write_text(
        "chr1\t24999\t25000\tgene_alpha\t0\t+\n"
        "chr1\t69999\t70000\tgene_shared\t0\t+\n"
        "chr1\t149999\t150000\tgene_beta\t0\t-\n",
        encoding="ascii",
    )
    (dm6_root / "dm6.autosomes.txt").write_text("chr1\n", encoding="ascii")

    fastq_root = runtime / "raw-fastq"
    manifest_rows = []
    raw_checksums: dict[Path, str] = {}
    for accession, library_id, assay, context in SAMPLE_ROWS:
        replicate = int(library_id.rsplit("rep", 1)[1])
        r1_path = (fastq_root / accession / f"{accession}_1.fastq").resolve()
        r2_path = (fastq_root / accession / f"{accession}_2.fastq").resolve()
        write_fastq_pair(
            chromosome=chromosome,
            context=context,
            assay=assay,
            replicate=replicate,
            r1_path=r1_path,
            r2_path=r2_path,
        )
        raw_checksums[r1_path] = sha256(r1_path)
        raw_checksums[r2_path] = sha256(r2_path)
        manifest_rows.append(
            {
                "requested_accession": accession,
                "experiment_accession": accession,
                "run_accession": accession,
                "library_layout": "PAIRED",
                "backend": "synthetic-local",
                "status": "existing",
                "fastq_1": str(r1_path),
                "fastq_2": str(r2_path),
                "extra_fastqs": "",
                "md5_1": md5(r1_path),
                "md5_2": md5(r2_path),
                "extra_md5s": "",
            }
        )

    sample_sheet = runtime / "synthetic-accessions.tsv"
    with sample_sheet.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("accession", "library_id", "assay", "context"))
        writer.writerows(SAMPLE_ROWS)

    manifest = runtime / "download_manifest.tsv"
    fields = tuple(manifest_rows[0])
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    return sample_sheet, manifest, reference_root, raw_checksums


def run(command: list[str], *, env: dict[str, str], capture: bool = False) -> str:
    print("[command] " + " ".join(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    if capture and completed.stdout:
        print(completed.stdout, end="", flush=True)
    return completed.stdout or ""


def patch_config(path: Path, output_root: Path) -> dict[str, object]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    config["output_dir"] = str(output_root.resolve())
    reference = config["reference"]
    reference.pop("preparation", None)
    reference["effective_genome_size"] = CHROMOSOME_LENGTH
    reference["macs3_genome_size"] = str(CHROMOSOME_LENGTH)
    config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
    path.write_text(
        yaml.safe_dump(config, sort_keys=False, default_flow_style=False),
        encoding="utf-8",
    )
    return config


def snakemake_command(config_path: Path, cores: int) -> list[str]:
    return [
        str(REPO_ROOT / ".venv" / "bin" / "snakemake"),
        "--snakefile",
        str(REPO_ROOT / "workflow" / "Snakefile"),
        "--configfile",
        str(config_path),
        "--workflow-profile",
        str(REPO_ROOT / "profiles" / "local"),
        "--cores",
        str(cores),
        "--jobs",
        str(cores),
        "--max-threads",
        str(min(4, cores)),
    ]


def result_root(config: dict[str, object]) -> Path:
    return (
        Path(str(config["output_dir"]))
        / str(config["project"])
        / str(config["run_id"])
    )


def review_qc(
    *,
    root: Path,
    label: str,
    runtime: Path,
    python: Path,
    env: dict[str, str],
) -> Path:
    review = root / "qc" / "library-review.tsv"
    manifest = root / "provenance" / "manifests" / "final-bams.tsv"
    decision_table = runtime / "reviewed" / f"{label}.library-review.tsv"
    decision_table.parent.mkdir(parents=True, exist_ok=True)
    with review.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if not rows or "qc_decision" not in fields:
        raise AssertionError(f"Invalid QC review table: {review}")
    for row in rows:
        row["qc_decision"] = "pass"
        row["estimated_fragment_length_bp"] = ""
        row["notes"] = "accepted only for the deterministic synthetic integration test"
    with decision_table.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(
            output, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    reviewed = runtime / "reviewed" / f"{label}.final-bams.tsv"
    run(
        [
            str(python),
            "src/review_final_bam_manifest.py",
            str(manifest),
            "--review-table",
            str(decision_table),
            "--output",
            str(reviewed),
        ],
        env=env,
    )
    artifacts = read_final_bam_manifest(reviewed, require_files=True)
    if not artifacts or {item["qc_status"] for item in artifacts.values()} != {"accepted"}:
        raise AssertionError(f"Synthetic review did not accept every library: {reviewed}")
    return reviewed


def check_corrupt_fastq_rejected(runtime: Path, python: Path, env: dict[str, str]) -> None:
    corrupt = runtime / "negative-control.invalid.fastq"
    corrupt.write_text("@bad\nACGT%\n+\nIIIII\n", encoding="ascii")
    completed = subprocess.run(
        [str(python), "workflow/scripts/validate_fastq.py", str(corrupt)],
        cwd=REPO_ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode == 0:
        raise AssertionError("Strict FASTQ validation accepted an invalid base")
    print("[check] strict FASTQ validation rejected the corrupt negative control", flush=True)


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def validate_outputs(
    *,
    initial_roots: dict[str, Path],
    master_root: Path,
    activity_root: Path,
    raw_checksums: dict[Path, str],
) -> dict[str, object]:
    print("[check] validating pipeline outputs and scientific annotations", flush=True)
    for path, before in raw_checksums.items():
        if sha256(path) != before:
            raise AssertionError(f"Raw FASTQ was modified in place: {path}")

    atac_manifest = initial_roots["atac"] / "provenance/manifests/final-bams.tsv"
    chip_manifest = initial_roots["chip_histone"] / "provenance/manifests/final-bams.tsv"
    atac = read_final_bam_manifest(atac_manifest, require_files=True, allow_rejected=True)
    chip = read_final_bam_manifest(chip_manifest, require_files=True, allow_rejected=True)
    if len(atac) != 4 or len(chip) != 4:
        raise AssertionError("Expected four final ATAC and four final H3K27ac BAMs")
    if {row["qc_status"] for row in (*atac.values(), *chip.values())} != {
        "pending_review"
    }:
        raise AssertionError("Unreviewed QC BAM manifests must be pending_review")

    for root in initial_roots.values():
        for relative in (
            "qc/library-review.tsv",
            "qc/metrics.tsv",
            "qc/multiqc/multiqc_report.html",
            "provenance/manifests/qc.checkpoint.json",
        ):
            path = root / relative
            if not path.is_file() or path.stat().st_size == 0:
                raise AssertionError(f"Missing QC output: {path}")

    master_manifest = master_root / "provenance/manifests/master-dhs.tsv"
    master = read_master_manifest(master_manifest, require_files=True)
    master_bed = Path(master["master_bed"])
    master_count = sum(1 for line in master_bed.read_text(encoding="utf-8").splitlines() if line)
    if master_count < 1:
        raise AssertionError("Synthetic ATAC data produced an empty master DHS registry")

    catalog = activity_root / "activity/catalog/master_elements_long.tsv.gz"
    columns, rows = read_tsv_gz(catalog)
    required = {
        "context",
        "blacklist_overlap",
        "blacklist_overlap_bp",
        "blacklist_overlap_fraction",
        "nearest_tss_distance_bp",
        "h3k27ac_max_500_background_percentile",
        "h3k27ac_active",
    }
    if not required.issubset(columns):
        raise AssertionError("Catalog lacks required annotations: " + ", ".join(sorted(required - set(columns))))
    if {row["context"] for row in rows} != set(CONTEXT_TARGETS):
        raise AssertionError("Catalog does not contain both synthetic contexts")
    blacklist_rows = [row for row in rows if row["blacklist_overlap"] == "1"]
    if not blacklist_rows:
        raise AssertionError("Blacklist-overlapping DHS was removed instead of annotated")
    if not any(row["nearest_tss_distance_bp"] for row in rows):
        raise AssertionError("Nearest-TSS distance annotation is empty")
    percentile_values = [
        float(row["h3k27ac_max_500_background_percentile"])
        for row in rows
        if row["h3k27ac_max_500_background_percentile"]
    ]
    if any(value < 0 or value > 1 for value in percentile_values):
        raise AssertionError("H3K27ac background percentile is outside [0, 1]")

    igv_root = activity_root / "activity/catalog/igv"
    sessions = [igv_root / "alpha.xml", igv_root / "beta.xml", activity_root / "activity/catalog/all-contexts.igv.xml"]
    for session in sessions:
        tracks = [
            track
            for track in ET.parse(session).getroot().findall("./Panel/Track")
            if track.attrib.get("id") != "Reference sequence"
        ]
        if not tracks or tracks[0].attrib.get("name") != "Master DHS registry":
            raise AssertionError(f"Master DHS is not the first user track in {session}")

    for relative in (
        "activity/report/integrated_qc_report.html",
        "activity/report/integrated_qc_report.pdf",
        "activity/report/integrated_qc_report.json",
        "provenance/manifests/report.checkpoint.json",
    ):
        path = activity_root / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise AssertionError(f"Missing final report output: {path}")

    return {
        "atac_final_bams": len(atac),
        "h3k27ac_final_bams": len(chip),
        "master_dhs": master_count,
        "catalog_rows": len(rows),
        "blacklist_annotated_rows": len(blacklist_rows),
        "nonblank_h3k27ac_background_percentiles": len(percentile_values),
        "contexts": sorted(CONTEXT_TARGETS),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runtime",
        type=Path,
        default=REPO_ROOT / "results" / "synthetic-e2e-runtime",
    )
    parser.add_argument("--cores", type=int, default=6)
    args = parser.parse_args()
    if args.cores < 1:
        parser.error("--cores must be positive")
    runtime = args.runtime.resolve()
    if runtime.exists():
        raise FileExistsError(
            f"Refusing to replace an existing runtime directory: {runtime}"
        )
    runtime.mkdir(parents=True)
    (runtime / ".synthetic-e2e-runtime").write_text("owned by tests/run_synthetic_e2e.py\n")
    (runtime / "tmp").mkdir()

    python = REPO_ROOT / ".venv" / "bin" / "python"
    snakemake = REPO_ROOT / ".venv" / "bin" / "snakemake"
    if not python.is_file() or not snakemake.is_file():
        raise FileNotFoundError("Repository-local .venv is not installed")
    env = os.environ.copy()
    env.update(
        {
            "MAMBA_ROOT_PREFIX": str(REPO_ROOT / ".micromamba"),
            "CONDA_PKGS_DIRS": str(REPO_ROOT / ".micromamba" / "pkgs"),
            "CONDA_ENVS_PATH": str(REPO_ROOT / ".micromamba" / "envs"),
            "CONDARC": str(REPO_ROOT / ".condarc"),
            "XDG_CACHE_HOME": str(REPO_ROOT / ".cache"),
            "TMPDIR": str(runtime / "tmp"),
        }
    )
    sample_sheet, manifest, reference_root, raw_checksums = generate_inputs(runtime)
    check_corrupt_fastq_rejected(runtime, python, env)

    output_root = runtime / "outputs"
    initial_config_dir = runtime / "configs" / "qc"
    run(
        [
            str(python),
            "src/run_pipeline.py",
            str(sample_sheet),
            "--project",
            "synthetic-e2e",
            "--run-id",
            "qc",
            "--genome",
            "dm6",
            "--manifest",
            str(manifest),
            "--skip-download",
            "--config-dir",
            str(initial_config_dir),
            "--reference-root",
            str(reference_root),
            "--until-stage",
            "qc",
            "--config-only",
        ],
        env=env,
    )
    initial_configs: dict[str, tuple[Path, dict[str, object]]] = {}
    for path in sorted(initial_config_dir.glob("*.yaml")):
        config = patch_config(path, output_root)
        initial_configs[str(config["assay"])] = (path, config)
    if set(initial_configs) != {"atac", "chip_histone"}:
        raise AssertionError(f"Unexpected assay configs: {sorted(initial_configs)}")

    completed_configs: list[Path] = []
    initial_roots: dict[str, Path] = {}
    for assay in ("atac", "chip_histone"):
        path, config = initial_configs[assay]
        print(f"[stage] FASTQ through QC: {assay}", flush=True)
        run(snakemake_command(path, args.cores), env=env)
        completed_configs.append(path)
        initial_roots[assay] = result_root(config)

    reviewed_atac = review_qc(
        root=initial_roots["atac"],
        label="atac",
        runtime=runtime,
        python=python,
        env=env,
    )
    reviewed_h3k27ac = review_qc(
        root=initial_roots["chip_histone"],
        label="h3k27ac",
        runtime=runtime,
        python=python,
        env=env,
    )

    master_config_dir = runtime / "configs" / "master"
    atac_checkpoint = initial_roots["atac"] / "provenance/manifests/qc.checkpoint.json"
    run(
        [
            str(python),
            "src/run_pipeline.py",
            str(sample_sheet),
            "--project",
            "synthetic-e2e-master",
            "--run-id",
            "master",
            "--genome",
            "dm6",
            "--from-stage",
            "qc",
            "--until-stage",
            "master",
            "--checkpoint-manifest",
            str(atac_checkpoint),
            "--final-bam-manifest",
            str(reviewed_atac),
            "--config-dir",
            str(master_config_dir),
            "--reference-root",
            str(reference_root),
            "--config-only",
        ],
        env=env,
    )
    master_config_path = next(master_config_dir.glob("*.yaml"))
    master_config = patch_config(master_config_path, output_root)
    print("[stage] reviewed ATAC QC through master DHS", flush=True)
    run(snakemake_command(master_config_path, args.cores), env=env)
    completed_configs.append(master_config_path)
    master_root = result_root(master_config)
    master_manifest = master_root / "provenance/manifests/master-dhs.tsv"

    activity_config_dir = runtime / "configs" / "activity"
    run(
        [
            str(python),
            "src/run_pipeline.py",
            str(sample_sheet),
            "--project",
            "synthetic-e2e-activity",
            "--run-id",
            "report",
            "--genome",
            "dm6",
            "--from-stage",
            "master",
            "--until-stage",
            "report",
            "--master-manifest",
            str(master_manifest),
            "--activity-bam-manifest",
            str(reviewed_atac),
            "--activity-bam-manifest",
            str(reviewed_h3k27ac),
            "--report-source-root",
            str(initial_roots["atac"]),
            "--report-source-root",
            str(initial_roots["chip_histone"]),
            "--config-dir",
            str(activity_config_dir),
            "--reference-root",
            str(reference_root),
            "--config-only",
        ],
        env=env,
    )
    activity_config_path = next(activity_config_dir.glob("*.yaml"))
    activity_config = patch_config(activity_config_path, output_root)
    print("[stage] master DHS and reviewed BAMs through catalog and report", flush=True)
    run(snakemake_command(activity_config_path, args.cores), env=env)
    completed_configs.append(activity_config_path)
    activity_root = result_root(activity_config)

    summary = validate_outputs(
        initial_roots=initial_roots,
        master_root=master_root,
        activity_root=activity_root,
        raw_checksums=raw_checksums,
    )

    checkpoint = activity_root / "provenance/manifests/report.checkpoint.json"
    checkpoint_mtime = checkpoint.stat().st_mtime_ns
    print("[restart] rerunning every completed stage to verify idempotent restart", flush=True)
    restart_outputs = []
    for config_path in completed_configs:
        restart_outputs.append(
            run(snakemake_command(config_path, args.cores), env=env, capture=True)
        )
    if checkpoint.stat().st_mtime_ns != checkpoint_mtime:
        raise AssertionError("No-op restart changed the report checkpoint mtime")
    if not all("Nothing to be done" in output for output in restart_outputs):
        raise AssertionError("At least one completed workflow was not a no-op on restart")
    for path, before in raw_checksums.items():
        if sha256(path) != before:
            raise AssertionError(f"No-op restart modified raw FASTQ: {path}")

    summary.update(
        {
            "runtime": str(runtime),
            "activity_result_root": str(activity_root),
            "raw_fastqs_unchanged": True,
            "corrupt_fastq_rejected": True,
            "idempotent_restart": True,
            "completed_at_unix": int(time.time()),
        }
    )
    summary_path = runtime / "synthetic-e2e-summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print("[pass] synthetic end-to-end pipeline test completed", flush=True)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
