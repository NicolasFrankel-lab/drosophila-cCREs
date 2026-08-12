#!/usr/bin/env python3
"""Run a deterministic local catalog-to-contact-tables integration test."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from short_read_processing.artifacts import sha256_file  # noqa: E402
from short_read_processing.contact_metadata import (  # noqa: E402
    CONTACT_SOURCE_COLUMNS,
    DM6_ATLAS_CONTACT_CONTEXTS,
)
from short_read_processing.workflow_config import (  # noqa: E402
    workflow_semantic_sha256,
)


CONTEXTS = tuple(row["id"] for row in DM6_ATLAS_CONTACT_CONTEXTS)
OBSERVED_CONTEXTS = tuple(
    row["id"]
    for row in DM6_ATLAS_CONTACT_CONTEXTS
    if row["strategy"] == "observed"
)
POWERLAW_CONTEXTS = tuple(
    row["id"]
    for row in DM6_ATLAS_CONTACT_CONTEXTS
    if row["strategy"] == "powerlaw"
)
CHROMOSOME = "chr2L"
CHROMOSOME_LENGTH = 200_000
MAXIMUM_DISTANCE = 100_000
SOURCE_PREFIX = "synthetic_contact_e2e_v3"
CATALOG_FIELDS = (
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "context",
    "context_membership",
    "regulatory_class",
    "atac_normalized_cpm_per_kb",
    "h3k27ac_max_500_background_percentile",
    "h3k27ac_active",
    "activity_state",
    "combined_activity_max_500",
    "blacklist_overlap",
)
ELEMENTS = (
    ("DHS_SYNTH_P1", 19_750, 20_250, 20_000, "promoter_associated", 0.95, 100.0, 0),
    ("DHS_SYNTH_E1", 44_800, 45_200, 45_000, "distal_enhancer_like", 0.90, 70.0, 0),
    ("DHS_SYNTH_E2", 104_800, 105_200, 105_000, "proximal_enhancer_like", 0.80, 60.0, 1),
    ("DHS_SYNTH_P2", 119_750, 120_250, 120_000, "promoter_associated", 0.90, 80.0, 0),
    ("DHS_SYNTH_E3", 159_800, 160_200, 160_000, "distal_enhancer_like", 0.40, 50.0, 0),
    ("DHS_SYNTH_P3", 179_750, 180_250, 180_000, "promoter_associated", 0.20, 20.0, 0),
)


def sha256(path: Path) -> str:
    return sha256_file(path)


def deterministic_gzip_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(content.encode("utf-8"))


def run(
    command: list[str],
    *,
    env: dict[str, str],
    capture: bool = False,
) -> str:
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


def find_contacts_python() -> Path:
    candidates = []
    for environment in sorted((REPO_ROOT / ".snakemake" / "conda").glob("*")):
        python = environment / "bin" / "python"
        cooler = environment / "bin" / "cooler"
        pandas_records = list((environment / "conda-meta").glob("pandas-2.2*.json"))
        if python.is_file() and cooler.is_file() and pandas_records:
            candidates.append(python)
    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one installed contacts rule environment; found "
            + ", ".join(map(str, candidates))
        )
    return candidates[0]


def write_sample_sheet(runtime: Path) -> Path:
    sample_sheet = runtime / "synthetic-canonical-accessions.tsv"
    rows = []
    accession = 99200000
    for context in CONTEXTS:
        for assay in ("atac", "h3k27ac"):
            accession += 1
            rows.append(
                {
                    "accession": f"SRR{accession}",
                    "library_id": f"{context}_{assay}_synthetic",
                    "assay": assay,
                    "context": context,
                }
            )
    with sample_sheet.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("accession", "library_id", "assay", "context"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return sample_sheet


def write_catalog_bundle(runtime: Path, sample_sheet: Path) -> tuple[Path, Path]:
    source_root = runtime / "catalog-source"
    catalog = source_root / "activity" / "catalog" / "master_elements_long.tsv.gz"
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=CATALOG_FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for context in CONTEXTS:
        for master_id, start, end, summit, regulatory_class, posterior, combined, blacklist in ELEMENTS:
            writer.writerow(
                {
                    "master_dhs_id": master_id,
                    "chrom": CHROMOSOME,
                    "start": start,
                    "end": end,
                    "summit": summit,
                    "context": context,
                    "context_membership": 1,
                    "regulatory_class": regulatory_class,
                    "atac_normalized_cpm_per_kb": 100.0,
                    "h3k27ac_max_500_background_percentile": posterior,
                    "h3k27ac_active": int(posterior > 0.60 and not blacklist),
                    "activity_state": (
                        "active_h3k27ac_background_gt_0.60"
                        if posterior > 0.60 and not blacklist
                        else "below_h3k27ac_background_threshold"
                    ),
                    "combined_activity_max_500": combined,
                    "blacklist_overlap": blacklist,
                }
            )
    deterministic_gzip_text(catalog, buffer.getvalue())

    catalog_digest = sha256(catalog)
    metrics = source_root / "activity" / "catalog" / "regulatory_element_metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "method": "synthetic_contact_catalog_v1",
                "context_count": len(CONTEXTS),
                "catalog_sha256": catalog_digest,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    provenance = source_root / "activity" / "catalog" / "regulatory_element_provenance.json"
    provenance.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "method": "synthetic_contact_catalog_v1",
                "outputs": {
                    "catalog": {"path": str(catalog), "sha256": catalog_digest}
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    resolved = source_root / "provenance" / "configs" / "report.resolved_config.json"
    resolved.parent.mkdir(parents=True, exist_ok=True)
    source_config: dict[str, Any] = {
        "project": "synthetic-contact-catalog",
        "run_id": "source",
        "output_dir": str(source_root),
        "assay": "activity",
        "input_stage": "quantification",
        "start_stage": "master",
        "output_stage": "report",
        "reference": {"name": "dm6"},
        "samples": [],
        "activity": {"contexts": list(CONTEXTS)},
        "provenance": {"sample_sheet_sha256": sha256(sample_sheet)},
    }
    source_config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(
        source_config
    )
    resolved.write_text(
        json.dumps(source_config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return source_root, catalog


def write_reference(runtime: Path) -> tuple[Path, Path]:
    reference_root = runtime / "references"
    dm6 = reference_root / "dm6"
    dm6.mkdir(parents=True, exist_ok=True)
    (dm6 / "dm6.chrom.sizes").write_text(
        f"{CHROMOSOME}\t{CHROMOSOME_LENGTH}\n", encoding="ascii"
    )
    (dm6 / "dm6.fa").write_text(
        f">{CHROMOSOME}\n"
        + "\n".join(
            "N" * min(80, CHROMOSOME_LENGTH - start)
            for start in range(0, CHROMOSOME_LENGTH, 80)
        )
        + "\n",
        encoding="ascii",
    )
    (dm6 / "dm6.blacklist.bed").write_text("", encoding="ascii")
    (dm6 / "dm6.tss.bed").write_text("", encoding="ascii")
    (dm6 / "dm6.autosomes.txt").write_text(f"{CHROMOSOME}\n", encoding="ascii")
    annotation = dm6 / "synthetic.gtf.gz"
    gtf = (
        'chr2L\tsynthetic\ttranscript\t20001\t23000\t.\t+\t.\tgene_id "geneA"; gene_name "GeneA"; transcript_id "txA";\n'
        'chr2L\tsynthetic\ttranscript\t117001\t120001\t.\t-\t.\tgene_id "geneB"; gene_name "GeneB"; transcript_id "txB";\n'
        'chr2L\tsynthetic\ttranscript\t180001\t183000\t.\t+\t.\tgene_id "geneC"; gene_name "GeneC"; transcript_id "txC";\n'
    )
    deterministic_gzip_text(annotation, gtf)
    return reference_root, annotation


def source_resolution(context: str) -> int:
    return 2_000 if context == "o" else 1_000


def target_resolution(context: str) -> int:
    return 4_000 if context == "o" else 5_000


def source_count(
    first: int,
    second: int,
    *,
    source_bin_size: int,
    target_bin_size: int,
    replicate: int,
) -> int:
    offset = second - first
    count = max(1, int(round((900 + 50 * replicate) / ((offset + 1) ** 1.15))))
    factor = target_bin_size // source_bin_size
    first_target = first // factor
    second_target = second // factor
    enriched_pairs = {
        tuple(sorted((20_000 // target_bin_size, 45_000 // target_bin_size))),
        tuple(sorted((105_000 // target_bin_size, 120_000 // target_bin_size))),
    }
    if (first_target, second_target) in enriched_pairs:
        count += 5_000
    return count


def write_contact_sources(runtime: Path) -> tuple[Path, dict[str, int], dict[str, list[Path]]]:
    import cooler
    import pandas as pd

    source_root = runtime / "contact-sources"
    source_root.mkdir(parents=True, exist_ok=True)
    repository_inputs = REPO_ROOT / "data" / "raw" / "contacts"
    repository_inputs.mkdir(parents=True, exist_ok=True)
    rows = []
    expected_counts: dict[str, int] = {}
    sources_by_context: dict[str, list[Path]] = {}
    validation_pair = (0, 10)

    for context in OBSERVED_CONTEXTS:
        source_bin_size = source_resolution(context)
        target_bin_size = target_resolution(context)
        source_bin_count = CHROMOSOME_LENGTH // source_bin_size
        bins = pd.DataFrame(
            {
                "chrom": [CHROMOSOME] * source_bin_count,
                "start": [index * source_bin_size for index in range(source_bin_count)],
                "end": [(index + 1) * source_bin_size for index in range(source_bin_count)],
            }
        )
        expected = 0
        sources_by_context[context] = []
        for replicate in (1, 2):
            pixels = []
            for first in range(source_bin_count):
                for second in range(first, source_bin_count):
                    count = source_count(
                        first,
                        second,
                        source_bin_size=source_bin_size,
                        target_bin_size=target_bin_size,
                        replicate=replicate,
                    )
                    pixels.append(
                        {"bin1_id": first, "bin2_id": second, "count": count}
                    )
                    factor = target_bin_size // source_bin_size
                    if (
                        first // factor == validation_pair[0]
                        and second // factor == validation_pair[1]
                    ):
                        expected += count
            source_id = f"{SOURCE_PREFIX}_{context}_r{replicate}"
            uncompressed = source_root / f"{source_id}.cool"
            compressed = source_root / f"{source_id}.cool.gz"
            cooler.create_cooler(
                str(uncompressed),
                bins,
                pd.DataFrame(pixels),
                ordered=True,
            )
            with uncompressed.open("rb") as source, gzip.open(compressed, "wb") as target:
                shutil.copyfileobj(source, target)
            uncompressed.unlink()
            expected_path = repository_inputs / f"{source_id}.cool.gz"
            if expected_path.exists() or expected_path.is_symlink():
                raise FileExistsError(f"Refusing to replace contact input: {expected_path}")
            expected_path.symlink_to(compressed)
            sources_by_context[context].append(expected_path)
            rows.append(
                {
                    "source_id": source_id,
                    "context": context,
                    "assay": "Micro-C",
                    "replicate": f"rep{replicate}",
                    "format": "cool.gz",
                    "url": f"https://example.invalid/{source_id}.cool.gz",
                    "local_path": f"data/raw/contacts/{source_id}.cool.gz",
                    "checksum": f"sha256:{sha256(compressed)}",
                    "match_quality": "synthetic_exact",
                    "biological_context": f"synthetic {context}",
                    "caveat": "synthetic integration fixture",
                }
            )
        expected_counts[context] = expected

    manifest = runtime / "synthetic-contact-sources.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CONTACT_SOURCE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    (runtime / "expected-merged-counts.json").write_text(
        json.dumps(expected_counts, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest, expected_counts, sources_by_context


def create_catalog_manifest(
    source_root: Path,
    runtime: Path,
    python: Path,
    env: dict[str, str],
) -> Path:
    output = runtime / "synthetic-catalog-manifest.tsv"
    run(
        [
            str(python),
            "src/create_catalog_manifest.py",
            str(source_root),
            "--output",
            str(output),
        ],
        env=env,
    )
    return output


def patch_config(
    path: Path,
    *,
    output_root: Path,
    reference_root: Path,
    annotation: Path,
    contact_manifest: Path,
) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    dm6 = reference_root / "dm6"
    config["output_dir"] = str(output_root.resolve())
    config["reference"] = {
        "name": "dm6",
        "fasta": str((dm6 / "dm6.fa").resolve()),
        "bowtie2_index": str((dm6 / "bowtie2" / "dm6").resolve()),
        "chrom_sizes": str((dm6 / "dm6.chrom.sizes").resolve()),
        "blacklist_bed": str((dm6 / "dm6.blacklist.bed").resolve()),
        "tss_bed": str((dm6 / "dm6.tss.bed").resolve()),
        "autosomes_file": str((dm6 / "dm6.autosomes.txt").resolve()),
        "mitochondrial_contig": "chrM",
        "effective_genome_size": CHROMOSOME_LENGTH,
        "macs3_genome_size": str(CHROMOSOME_LENGTH),
    }
    contacts = config["contacts"]
    contacts["source_manifest"] = str(contact_manifest.resolve())
    contacts["source_manifest_sha256"] = sha256(contact_manifest)
    contacts["promoter_annotation"] = str(annotation.resolve())
    contacts["promoter_annotation_checksum"] = f"sha256:{sha256(annotation)}"
    contacts["canonical_chromosomes"] = [CHROMOSOME]
    contacts["maximum_distance_bp"] = MAXIMUM_DISTANCE
    contacts["contexts"] = [
        {
            "id": row["id"],
            "strategy": row["strategy"],
            "assay": "Micro-C" if row["strategy"] == "observed" else "distance_model",
            "match": "synthetic_exact" if row["strategy"] == "observed" else "synthetic_distance_model",
            "resolution_bp": target_resolution(row["id"]),
            "caveat": "synthetic local integration fixture",
        }
        for row in DM6_ATLAS_CONTACT_CONTEXTS
    ]
    nearest_tss_links = config["nearest_tss_links"]
    nearest_tss_links["promoter_annotation"] = str(annotation.resolve())
    nearest_tss_links["promoter_annotation_checksum"] = (
        f"sha256:{sha256(annotation)}"
    )
    config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
    path.write_text(
        yaml.safe_dump(config, sort_keys=False, default_flow_style=False),
        encoding="utf-8",
    )
    return config


def snakemake_command(config: Path, cores: int) -> list[str]:
    return [
        str(REPO_ROOT / ".venv" / "bin" / "snakemake"),
        "--snakefile",
        str(REPO_ROOT / "workflow" / "Snakefile"),
        "--configfile",
        str(config),
        "--workflow-profile",
        str(REPO_ROOT / "profiles" / "local"),
        "--cores",
        str(cores),
        "--jobs",
        str(cores),
        "--resources",
        "mem_mb=16000",
        "contact_download_slots=1",
    ]


def count_gzip_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def validate_outputs(
    *,
    result: Path,
    runtime: Path,
    expected_counts: dict[str, int],
    source_checksums: dict[Path, str],
) -> dict[str, Any]:
    import cooler

    print("[check] validating normalized contacts and candidate tables", flush=True)
    links = result / "activity" / "links"
    normalized = links / "contacts"
    contexts_root = links / "contexts"
    retry_contexts = []
    for context in OBSERVED_CONTEXTS:
        cool_path = normalized / f"{context}.balanced.cool"
        metrics_path = normalized / f"{context}.metrics.json"
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        if metrics["replicate_count"] != 2 or not metrics["balance_converged"]:
            raise AssertionError(f"Contact normalization failed for {context}")
        if metrics["source_resolutions_bp"] != [source_resolution(context)] * 2:
            raise AssertionError(f"Unexpected source resolution for {context}")
        if metrics["target_resolution_bp"] != target_resolution(context):
            raise AssertionError(f"Unexpected target resolution for {context}")
        if metrics["balance_retry_used"]:
            retry_contexts.append(context)
        matrix = cooler.Cooler(str(cool_path)).matrix(balance=False)
        observed_count = int(matrix[0, 10].item())
        if observed_count != expected_counts[context]:
            raise AssertionError(
                f"Raw-count merge mismatch for {context}: "
                f"{observed_count} != {expected_counts[context]}"
            )

    powerlaw = json.loads((normalized / "dm6_powerlaw.json").read_text(encoding="utf-8"))
    if set(powerlaw["contexts"]) != set(OBSERVED_CONTEXTS):
        raise AssertionError("Power-law model lacks observed contexts")
    if not -3 < float(powerlaw["atlas_powerlaw"]["gamma"]) < -0.05:
        raise AssertionError("Synthetic atlas power-law exponent is implausible")

    promoter_metrics = json.loads((links / "promoters.metrics.json").read_text(encoding="utf-8"))
    if promoter_metrics["promoter_count"] != 3:
        raise AssertionError("Synthetic GTF did not produce three promoters")
    wide_fields, wide_rows = count_gzip_rows(
        links / "nearest_tss" / "enhancer_nearest_tss_candidates_wide.tsv.gz"
    )
    if (
        not wide_rows
        or len({row["master_dhs_id"] for row in wide_rows}) != len(wide_rows)
        or "ab__nearest_promoter_active" not in wide_fields
        or "nearest_gene_names" not in wide_fields
    ):
        raise AssertionError("Invalid one-row-per-enhancer nearest-TSS table")

    focused_counts: dict[str, int] = {}
    nearest_counts: dict[str, int] = {}
    distance_counts: dict[str, int] = {}
    for context in CONTEXTS:
        nodes_fields, nodes = count_gzip_rows(contexts_root / f"{context}.nodes.tsv.gz")
        edge_fields, edges = count_gzip_rows(
            contexts_root / f"{context}.element_promoter_edges.tsv.gz"
        )
        _gene_fields, genes = count_gzip_rows(
            contexts_root / f"{context}.element_gene_candidates.tsv.gz"
        )
        focused_fields, focused = count_gzip_rows(
            contexts_root / f"{context}.active_contact_enhancer_gene_candidates.tsv.gz"
        )
        if not nodes or not edges or not genes or "blacklist_overlap" not in nodes_fields:
            raise AssertionError(f"Incomplete contact graph for {context}")
        if context in OBSERVED_CONTEXTS:
            if not focused:
                raise AssertionError(f"Observed focused candidate table is empty: {context}")
            if any(
                float(row["element_h3k27ac_background_percentile"]) < 0.5
                or float(row["best_observed_over_expected"]) < 1.0
                for row in focused
            ):
                raise AssertionError(f"Focused thresholds were violated in {context}")
            if any(not row["observed_balanced_contact"] for row in edges):
                raise AssertionError(f"Observed context has blank observed contacts: {context}")
        else:
            if focused:
                raise AssertionError(f"Distance context focused contact table is not header-only: {context}")
            if any(
                row["observed_balanced_contact"] or row["observed_over_expected"]
                for row in edges
            ):
                raise AssertionError(f"Distance context contains observed contacts: {context}")
        focused_counts[context] = len(focused)

        nearest_fields, nearest_rows = count_gzip_rows(
            contexts_root
            / f"{context}.nearest_active_promoter_gene_candidates.tsv.gz"
        )
        if (
            not nearest_rows
            or "active_promoter_supporting_element_ids" not in nearest_fields
            or any(
                row["evidence_type"] != "nearest_active_promoter_tss"
                or not row["active_promoter_supporting_element_ids"]
                for row in nearest_rows
            )
        ):
            raise AssertionError(
                f"Invalid nearest-active-promoter candidates: {context}"
            )
        nearest_counts[context] = len(nearest_rows)

        distance_path = contexts_root / f"{context}.active_distance_enhancer_gene_candidates.tsv.gz"
        if context in POWERLAW_CONTEXTS:
            distance_fields, distance_rows = count_gzip_rows(distance_path)
            if not distance_rows or "evidence_type" not in distance_fields:
                raise AssertionError(f"Distance candidate table is empty: {context}")
            if any(
                row["evidence_type"] != "distance_model_active_promoter"
                or row["best_observed_balanced_contact"]
                or row["best_observed_over_expected"]
                for row in distance_rows
            ):
                raise AssertionError(f"Invalid distance evidence fields: {context}")
            distance_counts[context] = len(distance_rows)
        elif distance_path.exists():
            raise AssertionError(f"Observed context received a distance-only table: {context}")

    aggregate = json.loads((links / "contact_graph_metrics.json").read_text(encoding="utf-8"))
    if (
        aggregate["context_count"] != 9
        or aggregate["observed_context_count"] != 7
        or aggregate["powerlaw_context_count"] != 2
        or aggregate["active_contact_enhancer_gene_candidate_count"] < 1
        or aggregate["nearest_active_promoter_gene_candidate_count"] < 1
        or aggregate["active_distance_enhancer_gene_candidate_count"] < 1
    ):
        raise AssertionError("Aggregate contact metrics are incomplete")
    for required in (
        links / "contact_graph_provenance.json",
        result / "provenance" / "manifests" / "links.checkpoint.json",
    ):
        if not required.is_file() or required.stat().st_size == 0:
            raise AssertionError(f"Missing contact output: {required}")
    for source, digest in source_checksums.items():
        if sha256(source) != digest:
            raise AssertionError(f"Contact source was modified: {source}")

    work_root = REPO_ROOT / "work" / "synthetic-contacts-e2e" / "links-v3" / "activity" / "contacts"
    remaining_work_files = [path for path in work_root.rglob("*") if path.is_file()] if work_root.exists() else []
    if remaining_work_files:
        raise AssertionError("Contact normalization left intermediate files in work/")
    return {
        "contexts": list(CONTEXTS),
        "observed_contexts": list(OBSERVED_CONTEXTS),
        "powerlaw_contexts": list(POWERLAW_CONTEXTS),
        "normalized_contact_maps": len(OBSERVED_CONTEXTS),
        "exact_count_merge_verified": True,
        "balance_retry_contexts": retry_contexts,
        "focused_candidate_rows": focused_counts,
        "nearest_active_promoter_candidate_rows": nearest_counts,
        "distance_candidate_rows": distance_counts,
        "aggregate_metrics": aggregate,
    }


def worker(runtime: Path, cores: int) -> int:
    import cooler  # noqa: F401
    import pandas  # noqa: F401

    if runtime.exists():
        raise FileExistsError(f"Refusing to replace existing runtime: {runtime}")
    runtime.mkdir(parents=True)
    (runtime / ".synthetic-contacts-e2e-runtime").write_text(
        "owned by tests/run_synthetic_contacts_e2e.py\n", encoding="utf-8"
    )
    (runtime / "tmp").mkdir()
    python = REPO_ROOT / ".venv" / "bin" / "python"
    snakemake = REPO_ROOT / ".venv" / "bin" / "snakemake"
    if not python.is_file() or not snakemake.is_file():
        raise FileNotFoundError("Repository-local orchestration environment is absent")
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

    print("[setup] building synthetic canonical catalog and contact maps", flush=True)
    sample_sheet = write_sample_sheet(runtime)
    catalog_source, _catalog = write_catalog_bundle(runtime, sample_sheet)
    reference_root, annotation = write_reference(runtime)
    contact_manifest, expected_counts, sources_by_context = write_contact_sources(runtime)
    source_checksums = {
        path: sha256(path)
        for paths in sources_by_context.values()
        for path in paths
    }
    catalog_manifest = create_catalog_manifest(catalog_source, runtime, python, env)

    config_dir = runtime / "configs"
    run(
        [
            str(python),
            "src/run_pipeline.py",
            str(sample_sheet),
            "--from-stage",
            "catalog",
            "--until-stage",
            "links",
            "--catalog-manifest",
            str(catalog_manifest),
            "--project",
            "synthetic-contacts-e2e",
            "--run-id",
            "links-v3",
            "--genome",
            "dm6",
            "--with-contacts",
            "--reference-root",
            str(reference_root),
            "--config-dir",
            str(config_dir),
            "--config-only",
        ],
        env=env,
    )
    config_path = next(config_dir.glob("*.yaml"))
    config = patch_config(
        config_path,
        output_root=runtime / "outputs",
        reference_root=reference_root,
        annotation=annotation,
        contact_manifest=contact_manifest,
    )
    print("[stage] catalog through normalized contacts and candidate-gene tables", flush=True)
    command = snakemake_command(config_path, cores)
    run(command, env=env)
    result = Path(config["output_dir"]) / config["project"] / config["run_id"]
    summary = validate_outputs(
        result=result,
        runtime=runtime,
        expected_counts=expected_counts,
        source_checksums=source_checksums,
    )

    checkpoint = result / "provenance" / "manifests" / "links.checkpoint.json"
    checkpoint_mtime = checkpoint.stat().st_mtime_ns
    print("[restart] checking that the completed contact branch is a no-op", flush=True)
    restart = run(command, env=env, capture=True)
    if "Nothing to be done" not in restart:
        raise AssertionError("Completed contact workflow was not a no-op")
    if checkpoint.stat().st_mtime_ns != checkpoint_mtime:
        raise AssertionError("No-op restart changed the links checkpoint")
    summary.update(
        {
            "runtime": str(runtime),
            "result_root": str(result),
            "source_files_unchanged": True,
            "idempotent_restart": True,
        }
    )
    summary_path = runtime / "synthetic-contacts-e2e-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("[pass] synthetic contact-table integration test completed", flush=True)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runtime",
        type=Path,
        default=REPO_ROOT / "results" / "synthetic-contacts-e2e-runtime-v3",
    )
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.cores < 1:
        parser.error("--cores must be positive")
    runtime = args.runtime.resolve()
    if args.worker:
        return worker(runtime, args.cores)
    contacts_python = find_contacts_python()
    command = [
        str(contacts_python),
        str(Path(__file__).resolve()),
        "--worker",
        "--runtime",
        str(runtime),
        "--cores",
        str(args.cores),
    ]
    return subprocess.run(command, cwd=REPO_ROOT, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
