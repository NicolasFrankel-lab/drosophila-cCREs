"""Generate accession-first pipeline YAML from a validated sample sheet and manifest."""

from __future__ import annotations

from collections import defaultdict
import copy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any

import yaml

from .accessions import AcquisitionError
from .artifacts import (
    CATALOG_FILE_FIELDS,
    MASTER_FILE_FIELDS,
    read_catalog_manifest,
    read_final_bam_manifest,
    read_master_manifest,
    sha256_file,
)
from .manifest import read_manifest
from .sample_sheet import DEFAULT_SCHEMA, read_sample_sheet
from .stage_checkpoints import read_stage_checkpoint
from .workflow_config import validate_stage_selection, workflow_semantic_sha256
from .integrated_report import (
    discover_report_source_files,
    report_source_record,
)


GENOME_DEFAULTS = {
    "dm6": {"effective_genome_size": 142_573_017, "macs3_genome_size": "dm"},
}
REFERENCE_SOURCES = {
    "dm6": {
        "fasta": {
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz",
            "checksum": "md5:4a5777324403eff92c3650c30af30120",
        },
        "annotation": {
            "url": (
                "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/"
                "dm6.ncbiRefSeq.gtf.gz"
            ),
            "checksum": "md5:17fa021bf9cf35a12703787cc18ea4c8",
        },
        "blacklist": {
            "url": (
                "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/v2.0/lists/"
                "dm6-blacklist.v2.bed.gz"
            ),
            "checksum": (
                "sha256:6174aa0304859a2ee836de7f1eaabc4e7ab7fc4c75ffd4b7587172b1e234026e"
            ),
        },
        "autosomes": ["chr2L", "chr2R", "chr3L", "chr3R", "chr4"],
    },
}
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

ATAC_QPOIS_DEFAULTS = {
    "fragment_maximum": 150,
    "minimum_exponent": 2,
    "maximum_exponent": 325,
    "merge_gap": 1,
    "minimum_length": 50,
    "maximum_length": 400,
}
ATAC_MASTER_DEFAULTS = {
    "summit_max_distance": 150,
    "minimum_summit_separation": 50,
}


def _safe_id(value: str, label: str) -> str:
    if not SAFE_ID_RE.fullmatch(value):
        raise AcquisitionError(
            f"Invalid {label} {value!r}; use letters, digits, periods, underscores, or hyphens"
        )
    return value


def _display_path(value: str | Path, path_base: Path) -> str:
    path = Path(value).resolve()
    try:
        return path.relative_to(path_base.resolve()).as_posix()
    except ValueError:
        return str(path)


def _split_extra_fastqs(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def _sample_fastqs(
    manifest_rows: list[dict[str, str]],
    *,
    path_base: Path,
    require_files: bool,
) -> tuple[str, list[str], list[str]]:
    layouts = {row["library_layout"].upper() for row in manifest_rows}
    if len(layouts) != 1 or layouts.pop() not in {"SINGLE", "PAIRED"}:
        raise AcquisitionError("All runs for one sample must have the same SINGLE or PAIRED layout")
    layout = manifest_rows[0]["library_layout"].upper()

    r1 = [row["fastq_1"] for row in manifest_rows if row["fastq_1"]]
    r2 = [row["fastq_2"] for row in manifest_rows if row["fastq_2"]]
    if not r1 or (layout == "PAIRED" and len(r1) != len(r2)):
        raise AcquisitionError("Manifest does not contain a complete FASTQ set for the sample")
    paths = r1 + r2 + [
        item for row in manifest_rows for item in _split_extra_fastqs(row["extra_fastqs"])
    ]
    if require_files:
        missing = [item for item in paths if not Path(item).is_file()]
        if missing:
            raise AcquisitionError("Manifest FASTQ files are missing:\n  " + "\n  ".join(missing))
    return (
        "paired" if layout == "PAIRED" else "single",
        [_display_path(item, path_base) for item in r1],
        [_display_path(item, path_base) for item in r2],
    )


def _reference_config(genome: str, reference_root: Path, path_base: Path) -> dict[str, object]:
    root = reference_root / genome
    defaults = GENOME_DEFAULTS[genome]
    return {
        "name": genome,
        "fasta": _display_path(root / f"{genome}.fa", path_base),
        "bowtie2_index": _display_path(root / "bowtie2" / genome, path_base),
        "chrom_sizes": _display_path(root / f"{genome}.chrom.sizes", path_base),
        "blacklist_bed": _display_path(root / f"{genome}.blacklist.bed", path_base),
        "tss_bed": _display_path(root / f"{genome}.tss.bed", path_base),
        "autosomes_file": _display_path(root / f"{genome}.autosomes.txt", path_base),
        "mitochondrial_contig": "chrM",
        "effective_genome_size": defaults["effective_genome_size"],
        "macs3_genome_size": defaults["macs3_genome_size"],
        "preparation": {
            "mode": "download",
            **REFERENCE_SOURCES[genome],
        },
    }


def _nearest_tss_link_config(
    *, genome: str, reference: dict[str, Any], path_base: Path
) -> dict[str, Any]:
    """Freeze the nearest-TSS link defaults."""

    annotation_url = str(reference["preparation"]["annotation"]["url"])
    annotation_name = annotation_url.rsplit("/", 1)[-1]
    reference_fasta = Path(str(reference["fasta"]))
    if not reference_fasta.is_absolute():
        reference_fasta = path_base / reference_fasta
    annotation = reference_fasta.parent / "sources" / annotation_name
    return {
        "schema_version": 2,
        "promoter_annotation": _display_path(annotation, path_base),
        "promoter_annotation_checksum": str(
            reference["preparation"]["annotation"]["checksum"]
        ),
        "chromosome_scope": "all_annotation_chromosomes_in_chrom_sizes",
        "promoter_width_bp": 1000,
        "enhancer_classes": [
            "proximal_enhancer_like",
            "distal_enhancer_like",
        ],
        "promoter_activity": "summit_within_500bp_catalog_activity_v3",
        "candidate_assignment": "nearest_annotated_tss_all_ties_v1",
    }


def _peak_config(row: dict[str, Any], layout: str) -> dict[str, object]:
    assay = str(row["assay"])

    if assay == "atac":
        if row["macs3_broad"]:
            raise AcquisitionError("ATAC two-ended qpois calling does not support broad peaks")
        return {
            "command": "callpeak",
            "mode": "tn5_qpois",
            "format": "BED",
            "qvalue": row["macs3_qvalue"],
            "broad": False,
            "nomodel": True,
            "shift": row["macs3_shift"] if row["macs3_shift"] is not None else -75,
            "extsize": row["macs3_extsize"] if row["macs3_extsize"] is not None else 150,
            "write_bedgraph": True,
            "spmr": False,
        }

    macs_format = str(row["macs3_format"])
    if macs_format == "AUTO":
        macs_format = "BAMPE" if layout == "paired" else "BAM"
    config: dict[str, object] = {
        "command": "callpeak",
        "format": macs_format,
        "qvalue": row["macs3_qvalue"],
        "broad": row["macs3_broad"],
        "nomodel": row["macs3_nomodel"],
        "shift": row["macs3_shift"],
        "extsize": row["macs3_extsize"],
        "write_bedgraph": True,
        "spmr": True,
        "bedgraph_outputs": {
            "treatment_suffix": "_treat_pileup.bdg",
            "control_suffix": "_control_lambda.bdg",
        },
    }
    if row["macs3_broad"]:
        config["broad_cutoff"] = row["macs3_broad_cutoff"]
    return config


def _processing_config(row: dict[str, Any], path_base: Path) -> dict[str, object]:
    adapter_fasta = row.get("adapter_fasta")
    return {
        "trimming": {
            "adapter_preset": row["adapter_preset"],
            "adapter_fasta": _display_path(adapter_fasta, path_base) if adapter_fasta else None,
            "quality_cutoff": row["trim_quality_cutoff"],
            "minimum_length": row["trim_minimum_length"],
            "error_rate": 0.1,
            "minimum_overlap": 3,
        },
        "alignment": {
            "preset": row["bowtie2_preset"],
            "maximum_fragment_length": row["maximum_fragment_length"],
            "mapq_minimum": row["mapq_minimum"],
        },
        "filtering": {
            "remove_duplicates": row["remove_duplicates"],
            "remove_mitochondrial": row["remove_mitochondrial"],
        },
    }


def _manifest_runs_for_accessions(
    accessions: list[str], manifest_by_request: dict[str, list[dict[str, str]]]
) -> list[dict[str, str]]:
    by_run: dict[str, dict[str, str]] = {}
    for accession in accessions:
        runs = manifest_by_request.get(accession, [])
        if not runs:
            raise AcquisitionError(f"No manifest rows found for sample-sheet accession {accession}")
        for run in runs:
            by_run[run["run_accession"]] = run
    return [by_run[key] for key in sorted(by_run)]


def _write_config_if_changed(path: Path, config: dict[str, object]) -> None:
    """Atomically write a config without changing an identical file's mtime."""

    if path.is_file():
        existing = yaml.safe_load(path.read_text(encoding="utf-8"))
        if isinstance(existing, dict):
            existing_provenance = existing.get("provenance")
            provenance = config.get("provenance")
            if isinstance(existing_provenance, dict) and isinstance(provenance, dict):
                new_generated_at = provenance.get("generated_at_utc")
                generated_at = existing_provenance.get("generated_at_utc")
                if generated_at:
                    provenance["generated_at_utc"] = generated_at
                if existing == config:
                    return
                provenance["generated_at_utc"] = new_generated_at
            elif existing == config:
                return

    content = yaml.safe_dump(config, sort_keys=False, default_flow_style=False)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{path.name}.",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def generate_configs(
    *,
    manifest_path: Path | None,
    sample_sheet_path: Path,
    output_dir: Path,
    project: str,
    run_id: str,
    reference_root: Path,
    path_base: Path,
    require_fastq_files: bool,
    schema_path: Path = DEFAULT_SCHEMA,
    genome: str = "dm6",
    atac_minimum_replicates: int = 2,
    atac_overlap_fraction: float = 0.5,
    alignment_chunk_pairs: int | None = None,
    input_stage: str = "accessions",
    start_stage: str | None = None,
    output_stage: str | None = None,
    final_bam_manifest_path: Path | None = None,
    master_manifest_path: Path | None = None,
    qc_checkpoint_manifest_path: Path | None = None,
) -> list[Path]:
    project = _safe_id(project, "project ID")
    run_id = _safe_id(run_id, "run ID")
    if genome not in GENOME_DEFAULTS:
        raise AcquisitionError(f"Unsupported genome: {genome!r}")
    if alignment_chunk_pairs is not None and alignment_chunk_pairs < 1:
        raise AcquisitionError("Alignment chunk pairs must be positive")
    if input_stage not in {"accessions", "final-bam", "master"}:
        raise AcquisitionError(f"Unsupported input stage: {input_stage!r}")
    if alignment_chunk_pairs is not None and input_stage != "accessions":
        raise AcquisitionError(
            "Alignment chunking requires accession/FASTQ input"
        )
    start_stage = start_stage or input_stage
    output_stage = validate_stage_selection(start_stage, output_stage)
    if input_stage == "accessions" and manifest_path is None:
        raise AcquisitionError("Accession mode requires a download manifest")
    if input_stage == "final-bam" and final_bam_manifest_path is None:
        raise AcquisitionError("Final-BAM mode requires a final-BAM manifest")
    if input_stage == "master" and master_manifest_path is None:
        raise AcquisitionError("Master mode requires a master manifest")
    if start_stage == "qc" and qc_checkpoint_manifest_path is None:
        raise AcquisitionError("QC reuse requires a QC checkpoint manifest")

    sheet_rows = read_sample_sheet(sample_sheet_path, schema_path=schema_path)
    sample_sheet_sha256 = sha256_file(sample_sheet_path)
    schema_sha256 = sha256_file(schema_path)
    qc_checkpoint = (
        read_stage_checkpoint(
            qc_checkpoint_manifest_path,
            expected_stage="qc",
            require_files=require_fastq_files,
        )
        if start_stage == "qc"
        else None
    )
    if qc_checkpoint is not None and qc_checkpoint["parameters"].get("assay") != "atac":
        raise AcquisitionError("Master construction requires an ATAC QC checkpoint")

    if input_stage == "master":
        if not any(row["assay"] == "atac" for row in sheet_rows):
            raise AcquisitionError("Master reuse requires ATAC metadata in the sample sheet")
        master = read_master_manifest(
            master_manifest_path,
            require_files=require_fastq_files,
        )
        if master["genome"] != genome:
            raise AcquisitionError(
                f"Master manifest genome {master['genome']!r} does not match {genome!r}"
            )
        assays = {str(row["assay"]) for row in sheet_rows}
        group_project = project if assays == {"atac"} else f"{project}.atac.{genome}"
        external_master = {
            key: (
                _display_path(value, path_base)
                if key
                in {
                    "master_bed",
                    "summits_bed",
                    "membership_tsv",
                    "context_matrix_tsv",
                    "stats_json",
                }
                else value
            )
            for key, value in master.items()
        }
        config: dict[str, Any] = {
            "project": group_project,
            "run_id": run_id,
            "output_dir": "results",
            "assay": "atac",
            "input_stage": "master",
            "start_stage": start_stage,
            "output_stage": output_stage,
            "reference": _reference_config(genome, reference_root, path_base),
            "samples": [],
            "external_master": external_master,
            "provenance": {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "input_mode": "external_master_manifest",
                "sample_sheet": _display_path(sample_sheet_path, path_base),
                "sample_sheet_sha256": sample_sheet_sha256,
                "sample_sheet_schema": _display_path(schema_path, path_base),
                "sample_sheet_schema_sha256": schema_sha256,
                "master_manifest": _display_path(master_manifest_path, path_base),
                "master_manifest_sha256": sha256_file(master_manifest_path),
            },
        }
        config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{group_project}.yaml"
        _write_config_if_changed(output_path, config)
        return [output_path.resolve()]

    manifest_rows = read_manifest(manifest_path) if input_stage == "accessions" else []
    manifest_by_request: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        manifest_by_request[row["requested_accession"]].append(row)

    reviewed_master_exclusions: list[dict[str, str]] = []
    final_bams = (
        read_final_bam_manifest(
            final_bam_manifest_path,
            require_files=require_fastq_files,
            allow_rejected=output_stage == "master",
        )
        if input_stage == "final-bam"
        else {}
    )
    if final_bams:
        sheet_library_assay = {
            str(row["library_id"]): str(row["assay"]) for row in sheet_rows
        }
        unknown = sorted(set(final_bams) - set(sheet_library_assay))
        if unknown:
            raise AcquisitionError(
                "Final-BAM manifest contains libraries absent from the sample sheet: "
                + ", ".join(unknown)
            )
        selected_assays = {
            sheet_library_assay[library_id] for library_id in final_bams
        }
        expected = {
            library_id
            for library_id, assay in sheet_library_assay.items()
            if assay in selected_assays
        }
        missing = sorted(expected - set(final_bams))
        if missing:
            raise AcquisitionError(
                "Final-BAM manifest is incomplete for its selected assay(s): "
                + ", ".join(missing)
            )
        sheet_by_library_id = {
            str(row["library_id"]): row for row in sheet_rows
        }
        for library_id, artifact in final_bams.items():
            sheet_row = sheet_by_library_id[library_id]
            for field in ("assay", "context", "role"):
                if str(artifact[field]) != str(sheet_row[field]):
                    raise AcquisitionError(
                        f"Final-BAM manifest {field} for {library_id!r} "
                        "does not match the sample sheet"
                    )
            if artifact["genome"] != genome:
                raise AcquisitionError(
                    f"Final-BAM manifest genome for {library_id!r} "
                    f"is {artifact['genome']!r}, expected {genome!r}"
                )
        if output_stage == "master":
            if selected_assays != {"atac"}:
                raise AcquisitionError(
                    "Master DHS construction requires an ATAC-only final-BAM manifest"
                )
            pending = sorted(
                library_id
                for library_id, artifact in final_bams.items()
                if artifact["qc_status"] == "pending_review"
            )
            if pending:
                raise AcquisitionError(
                    "Master DHS construction requires manual QC decisions for every "
                    "ATAC library; pending_review: " + ", ".join(pending)
                )
            reviewed_master_exclusions = [
                {
                    "id": library_id,
                    "assay": "atac",
                    "context": artifact["context"],
                    "role": artifact["role"],
                    "layout": artifact["layout"],
                    "qc_status": "rejected",
                    "reason": artifact["notes"],
                }
                for library_id, artifact in sorted(final_bams.items())
                if artifact["qc_status"] == "rejected"
            ]
            final_bams = {
                library_id: artifact
                for library_id, artifact in final_bams.items()
                if artifact["qc_status"] == "accepted"
            }
            if not final_bams:
                raise AcquisitionError(
                    "Master DHS construction has no accepted ATAC libraries"
                )
        sheet_rows = [
            row for row in sheet_rows if str(row["library_id"]) in final_bams
        ]

    sheet_by_library: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in sheet_rows:
        sheet_by_library[str(row["library_id"])].append(row)

    prepared: list[dict[str, Any]] = []
    for library_id, rows in sheet_by_library.items():
        first = rows[0]
        accessions = list(dict.fromkeys(str(row["accession"]) for row in rows))
        artifact: dict[str, str] | None = None
        if input_stage == "accessions":
            runs = _manifest_runs_for_accessions(accessions, manifest_by_request)
            layout, r1, r2 = _sample_fastqs(
                runs,
                path_base=path_base,
                require_files=require_fastq_files,
            )
        else:
            artifact = final_bams[library_id]
            for field in ("assay", "context", "role"):
                if str(artifact[field]) != str(first[field]):
                    raise AcquisitionError(
                        f"Final-BAM manifest {field} for {library_id!r} "
                        "does not match the sample sheet"
                    )
            if artifact["genome"] != genome:
                raise AcquisitionError(
                    f"Final-BAM manifest genome for {library_id!r} "
                    f"is {artifact['genome']!r}, expected {genome!r}"
                )
            layout = artifact["layout"]
            r1, r2 = [], []
        prepared.append(
            {
                "id": library_id,
                "accessions": accessions,
                "assay": first["assay"],
                "genome": genome,
                "context": first["context"],
                "role": first["role"],
                "layout": layout,
                "r1": r1,
                "r2": r2,
                "peak_caller": _peak_config(first, layout),
                "parameters": _processing_config(first, path_base),
                "artifact": artifact,
            }
        )

    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for item in prepared:
        groups[(str(item["assay"]), str(item["genome"]))].append(item)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths: list[Path] = []
    for (assay, genome), items in sorted(groups.items()):
        sample_ids = {str(item["id"]) for item in items}
        item_by_id = {str(item["id"]): item for item in items}
        samples: list[dict[str, object]] = []
        for item in items:
            sample: dict[str, object] = {
                "id": item["id"],
                "accessions": item["accessions"],
                "context": item["context"],
                "role": item["role"],
                "layout": item["layout"],
                "parameters": item["parameters"],
            }
            if input_stage == "accessions":
                sample["r1"] = item["r1"]
                if item["layout"] == "paired":
                    sample["r2"] = item["r2"]
            else:
                artifact = dict(item["artifact"])
                artifact["bam"] = _display_path(artifact["bam"], path_base)
                artifact["bai"] = _display_path(artifact["bai"], path_base)
                sample["final_bam"] = artifact
            if item["role"] == "treatment":
                sample["peak_caller"] = item["peak_caller"]
                if qc_checkpoint is not None:
                    checkpoint_callers = qc_checkpoint["parameters"].get(
                        "peak_callers", {}
                    )
                    if checkpoint_callers.get(str(item["id"])) != item["peak_caller"]:
                        raise AcquisitionError(
                            f"QC checkpoint peak parameters for {item['id']!r} differ "
                            "from the requested master configuration"
                        )
                    label = f"replicate_peak.{item['id']}"
                    if label not in qc_checkpoint["artifacts"]:
                        raise AcquisitionError(
                            f"QC checkpoint lacks replicate peak evidence for {item['id']!r}"
                        )
                    record = qc_checkpoint["artifacts"][label]
                    sample["qc_peak"] = {
                        "path": _display_path(record["path"], path_base),
                        "sha256": record["sha256"],
                        "method": item["peak_caller"]["command"],
                    }
            samples.append(sample)

        group_project = project if len(groups) == 1 else f"{project}.{assay}.{genome}"
        config = {
            "project": group_project,
            "run_id": run_id,
            "output_dir": "results",
            "assay": assay,
            "input_stage": input_stage,
            "start_stage": start_stage,
            "output_stage": output_stage,
            "reference": _reference_config(genome, reference_root, path_base),
            "samples": samples,
            "provenance": {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "input_mode": (
                    "accession_sample_sheet"
                    if input_stage == "accessions"
                    else "external_final_bam_manifest"
                ),
                "sample_sheet": _display_path(sample_sheet_path, path_base),
                "sample_sheet_sha256": sample_sheet_sha256,
                "sample_sheet_schema": _display_path(schema_path, path_base),
                "sample_sheet_schema_sha256": schema_sha256,
            },
        }
        if alignment_chunk_pairs is not None:
            config["execution"] = {
                "alignment_chunk_pairs": alignment_chunk_pairs,
            }
        if input_stage == "accessions":
            config["provenance"]["download_manifest"] = _display_path(
                manifest_path, path_base
            )
            config["provenance"]["download_manifest_sha256"] = sha256_file(
                manifest_path
            )
        else:
            config["provenance"]["final_bam_manifest"] = _display_path(
                final_bam_manifest_path, path_base
            )
            config["provenance"]["final_bam_manifest_sha256"] = sha256_file(
                final_bam_manifest_path
            )
            if output_stage == "master":
                config["provenance"]["excluded_master_libraries"] = (
                    reviewed_master_exclusions
                )
            if qc_checkpoint is not None:
                config["provenance"]["qc_checkpoint_manifest"] = _display_path(
                    qc_checkpoint_manifest_path, path_base
                )
                config["provenance"]["qc_checkpoint_manifest_sha256"] = sha256_file(
                    qc_checkpoint_manifest_path
                )
        if assay == "atac":
            from .consensus import condition_specs

            config["atac_qpois"] = dict(ATAC_QPOIS_DEFAULTS)
            if qc_checkpoint is not None and qc_checkpoint["parameters"].get(
                "atac_qpois"
            ) != config["atac_qpois"]:
                raise AcquisitionError(
                    "QC checkpoint qpois-refinement parameters differ from the "
                    "requested master configuration"
                )
            samples_by_context: dict[str, list[str]] = defaultdict(list)
            for item in items:
                samples_by_context[str(item["context"])].append(str(item["id"]))
            condition_values = [
                {
                    "id": context,
                    "label": context,
                    "samples": sorted(context_samples),
                }
                for context, context_samples in sorted(samples_by_context.items())
            ]
            conditions = condition_specs(
                condition_values,
                sample_ids=sample_ids,
                minimum_replicates=atac_minimum_replicates,
            )
            for condition in conditions:
                layouts = {
                    str(item_by_id[sample]["layout"])
                    for sample in condition.samples
                }
                if len(layouts) != 1:
                    raise AcquisitionError(
                        f"ATAC condition {condition.condition_id} mixes read layouts"
                    )
            if not 0 < atac_overlap_fraction <= 1:
                raise AcquisitionError("ATAC overlap fraction must be in (0, 1]")
            config["atac_consensus"] = {
                "enabled": True,
                "conditions": condition_values,
                "minimum_replicates": atac_minimum_replicates,
                "replicate_overlap_fraction": atac_overlap_fraction,
            }
            config["atac_master"] = dict(ATAC_MASTER_DEFAULTS)
        config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
        output_path = output_dir / f"{group_project}.yaml"
        _write_config_if_changed(output_path, config)
        output_paths.append(output_path.resolve())
    return output_paths


def _merge_final_bam_manifests(
    paths: list[Path],
    *,
    require_files: bool,
    allow_rejected: bool = False,
) -> dict[str, dict[str, str]]:
    if not paths:
        raise AcquisitionError("At least one final-BAM manifest is required")
    merged: dict[str, dict[str, str]] = {}
    for path in paths:
        rows = read_final_bam_manifest(
            path,
            require_files=require_files,
            allow_rejected=allow_rejected,
        )
        duplicate = sorted(set(merged) & set(rows))
        if duplicate:
            raise AcquisitionError(
                "Final-BAM library occurs in more than one manifest: "
                + ", ".join(duplicate)
            )
        merged.update(rows)
    return merged


def _activity_libraries(
    *,
    sheet_rows: list[dict[str, Any]],
    manifests: dict[str, dict[str, str]],
    genome: str,
    path_base: Path,
) -> tuple[list[dict[str, Any]], list[str], list[dict[str, Any]]]:
    sheet_by_library: dict[str, dict[str, Any]] = {}
    for row in sheet_rows:
        sheet_by_library.setdefault(str(row["library_id"]), row)
    unknown = sorted(set(manifests) - set(sheet_by_library))
    if unknown:
        raise AcquisitionError(
            "Activity BAM manifest contains libraries absent from the sample sheet: "
            + ", ".join(unknown)
        )

    selected = {
        library_id: row
        for library_id, row in sheet_by_library.items()
        if row["role"] == "treatment"
        and row["assay"] in {"atac", "chip_histone"}
    }
    missing = sorted(set(selected) - set(manifests))
    if missing:
        raise AcquisitionError(
            "Activity final-BAM manifests are missing treatment libraries: "
            + ", ".join(missing)
        )
    libraries: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for library_id, row in sorted(selected.items()):
        artifact = manifests[library_id]
        for field in ("context", "role"):
            if str(artifact[field]) != str(row[field]):
                raise AcquisitionError(
                    f"Activity BAM {library_id!r} {field} differs from its sample sheet"
                )
        if artifact["assay"] != row["assay"]:
            raise AcquisitionError(
                f"Activity BAM {library_id!r} assay differs from its sample sheet"
            )
        if artifact["genome"] != genome:
            raise AcquisitionError(
                f"Activity BAM {library_id!r} uses {artifact['genome']!r}, "
                f"expected {genome!r}"
            )
        qc_status = artifact["qc_status"]
        if qc_status == "pending_review":
            raise AcquisitionError(
                f"Activity BAM {library_id!r} requires qc_status='accepted'"
            )
        assay = "h3k27ac" if artifact["assay"] == "chip_histone" else "atac"
        estimated_fragment_length = artifact.get(
            "estimated_fragment_length_bp", ""
        )
        if artifact["layout"] == "single" and assay == "h3k27ac":
            if not estimated_fragment_length:
                raise AcquisitionError(
                    f"Single-end H3K27ac BAM {library_id!r} requires "
                    "estimated_fragment_length_bp"
                )
        elif estimated_fragment_length:
            raise AcquisitionError(
                f"Activity BAM {library_id!r} must leave "
                "estimated_fragment_length_bp blank unless it is single-end H3K27ac"
            )
        if qc_status == "rejected":
            exclusions.append(
                {
                    "id": library_id,
                    "assay": assay,
                    "context": str(row["context"]),
                    "cohort": "atlas",
                    "layout": artifact["layout"],
                    "qc_status": qc_status,
                    "estimated_fragment_length_bp": (
                        int(estimated_fragment_length)
                        if estimated_fragment_length
                        else None
                    ),
                    "reason": artifact["notes"],
                }
            )
            continue
        if assay == "atac" and artifact["layout"] != "paired":
            raise AcquisitionError(
                f"Activity ATAC BAM {library_id!r} must be paired-end"
            )
        library: dict[str, Any] = {
            "id": library_id,
            "assay": assay,
            "context": str(row["context"]),
            "cohort": "atlas",
            "layout": artifact["layout"],
            "genome": artifact["genome"],
            "bam": _display_path(artifact["bam"], path_base),
            "bai": _display_path(artifact["bai"], path_base),
            "bam_sha256": artifact["bam_sha256"],
            "bai_sha256": artifact["bai_sha256"],
            "filtering_contract": artifact["filtering_contract"],
            "qc_status": artifact["qc_status"],
            "source_project": artifact.get("source_project", ""),
            "source_run_id": artifact.get("source_run_id", ""),
        }
        if estimated_fragment_length:
            library["estimated_fragment_length_bp"] = int(
                estimated_fragment_length
            )
        libraries.append(library)
    contexts = sorted({library["context"] for library in libraries})
    return libraries, contexts, exclusions


def generate_activity_config(
    *,
    sample_sheet_path: Path,
    final_bam_manifests: list[Path],
    master_manifest_path: Path,
    output_dir: Path,
    project: str,
    run_id: str,
    reference_root: Path,
    path_base: Path,
    require_files: bool,
    schema_path: Path = DEFAULT_SCHEMA,
    genome: str = "dm6",
    output_stage: str = "report",
    start_stage: str = "master",
    report_source_roots: list[Path] | None = None,
    h3k27ac_activity_threshold: float = 0.60,
) -> Path:
    """Generate one resolved master-DHS quantification/catalog configuration."""

    if not 0 < float(h3k27ac_activity_threshold) < 1:
        raise AcquisitionError(
            "H3K27ac activity percentile threshold must be in (0, 1)"
        )
    h3k27ac_activity_threshold = float(h3k27ac_activity_threshold)

    project = _safe_id(project, "project ID")
    run_id = _safe_id(run_id, "run ID")
    output_stage = validate_stage_selection(start_stage, output_stage)
    if genome not in GENOME_DEFAULTS:
        raise AcquisitionError(f"Unsupported genome: {genome!r}")
    sheet_rows = read_sample_sheet(
        sample_sheet_path,
        schema_path=schema_path,
    )
    manifests = _merge_final_bam_manifests(
        final_bam_manifests,
        require_files=require_files,
        allow_rejected=True,
    )
    libraries, contexts, exclusions = _activity_libraries(
        sheet_rows=sheet_rows,
        manifests=manifests,
        genome=genome,
        path_base=path_base,
    )
    artifact_paths = [
        library[field]
        for library in libraries
        for field in ("bam", "bai")
    ]
    if len(artifact_paths) != len(set(artifact_paths)):
        raise AcquisitionError(
            "Activity libraries must use distinct BAM and BAI artifact paths"
        )
    for context in contexts:
        assays = {
            library["assay"]
            for library in libraries
            if library["context"] == context
        }
        if assays != {"atac", "h3k27ac"}:
            raise AcquisitionError(
                f"Context {context!r} does not have both ATAC and H3K27ac"
            )
    for assay in ("atac", "h3k27ac"):
        assay_count = sum(library["assay"] == assay for library in libraries)
        if assay_count < 2:
            raise AcquisitionError(
                f"Background TMM requires at least two accepted {assay} libraries"
            )

    master = read_master_manifest(
        master_manifest_path,
        require_files=require_files,
    )
    if master["genome"] != genome:
        raise AcquisitionError(
            f"Master manifest genome {master['genome']!r} does not match {genome!r}"
        )
    activity_master = {
        key: (
            _display_path(value, path_base)
            if key in MASTER_FILE_FIELDS
            else value
        )
        for key, value in master.items()
    }
    inferred_source_roots: list[Path] = []
    resolved_master_manifest = master_manifest_path.resolve()
    for parent in resolved_master_manifest.parents:
        if parent.name == "provenance":
            inferred_source_roots.append(parent.parent)
            break
    source_roots = sorted(
        {
            path.resolve()
            for path in [*(report_source_roots or []), *inferred_source_roots]
        }
    )
    report_records = {
        record["path"]: record
        for record in discover_report_source_files(source_roots)
    }
    for path, kind in (
        (sample_sheet_path, "sample_sheet"),
        (schema_path, "sample_sheet_schema"),
        (master_manifest_path, "master_manifest"),
        *((path, "final_bam_manifest") for path in final_bam_manifests),
    ):
        record = report_source_record(
            path,
            kind=kind,
            source_root=path.resolve().parent,
        )
        report_records[record["path"]] = record
    provenance: dict[str, Any] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_mode": "quantification_artifact_manifests",
        "sample_sheet": _display_path(sample_sheet_path, path_base),
        "sample_sheet_sha256": sha256_file(sample_sheet_path),
        "sample_sheet_schema": _display_path(schema_path, path_base),
        "sample_sheet_schema_sha256": sha256_file(schema_path),
        "master_manifest": _display_path(master_manifest_path, path_base),
        "master_manifest_sha256": sha256_file(master_manifest_path),
        "final_bam_manifests": [
            _display_path(path, path_base) for path in final_bam_manifests
        ],
        "final_bam_manifest_sha256": [
            sha256_file(path) for path in final_bam_manifests
        ],
        "excluded_activity_libraries": exclusions,
    }
    config: dict[str, Any] = {
        "project": project,
        "run_id": run_id,
        "output_dir": "results",
        "assay": "activity",
        "input_stage": "quantification",
        "start_stage": start_stage,
        "output_stage": output_stage,
        "reference": _reference_config(genome, reference_root, path_base),
        "samples": [],
        "activity": {
            "schema_version": 6,
            "master": activity_master,
            "contexts": contexts,
            "libraries": libraries,
            "atac_fragment_maximum": 150,
            "atac_browser_extension_bp": 150,
            "normalization": "background_tmm_10kb_autosomes_v1",
            "h3k27ac_signal": "summit_max3_background_tmm_v3",
            "h3k27ac_activity_percentile_threshold": h3k27ac_activity_threshold,
        },
        "report": {
            "schema_version": 1,
            "source_roots": [
                _display_path(path, path_base) for path in source_roots
            ],
            "source_files": [
                {
                    **record,
                    "path": _display_path(record["path"], path_base),
                    "source_root": _display_path(record["source_root"], path_base),
                }
                for _path, record in sorted(report_records.items())
            ],
        },
        "provenance": provenance,
    }
    config["nearest_tss_links"] = _nearest_tss_link_config(
        genome=genome,
        reference=config["reference"],
        path_base=path_base,
    )
    config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{project}.quantification.yaml"
    _write_config_if_changed(output_path, config)
    return output_path.resolve()


def generate_catalog_links_config(
    *,
    sample_sheet_path: Path,
    catalog_manifest_path: Path,
    output_dir: Path,
    project: str,
    run_id: str,
    reference_root: Path,
    path_base: Path,
    schema_path: Path = DEFAULT_SCHEMA,
    genome: str = "dm6",
) -> Path:
    """Generate a links-only configuration from a frozen catalog bundle."""

    project = _safe_id(project, "project ID")
    run_id = _safe_id(run_id, "run ID")
    if genome not in GENOME_DEFAULTS:
        raise AcquisitionError(f"Unsupported genome: {genome!r}")
    read_sample_sheet(sample_sheet_path, schema_path=schema_path)
    catalog = read_catalog_manifest(catalog_manifest_path, require_files=True)
    if catalog["genome"] != genome:
        raise AcquisitionError(
            f"Catalog manifest genome {catalog['genome']!r} does not match {genome!r}"
        )
    sample_sheet_digest = sha256_file(sample_sheet_path)
    if sample_sheet_digest != catalog["source_sample_sheet_sha256"]:
        raise AcquisitionError(
            "The supplied sample sheet differs from the catalog source sample sheet"
        )
    catalog_import: dict[str, Any] = {
        "schema_version": 1,
        "manifest": _display_path(catalog_manifest_path, path_base),
        "manifest_sha256": sha256_file(catalog_manifest_path),
        "genome": catalog["genome"],
        "method": catalog["method"],
        "contexts": catalog["contexts"],
        "source_project": catalog["source_project"],
        "source_run_id": catalog["source_run_id"],
        "source_semantic_sha256": catalog["source_semantic_sha256"],
        "source_sample_sheet_sha256": catalog["source_sample_sheet_sha256"],
    }
    for field in CATALOG_FILE_FIELDS:
        catalog_import[field] = _display_path(catalog[field], path_base)
        catalog_import[f"{field}_sha256"] = catalog[f"{field}_sha256"]
    config: dict[str, Any] = {
        "project": project,
        "run_id": run_id,
        "output_dir": "results",
        "assay": "activity",
        "input_stage": "catalog",
        "start_stage": "catalog",
        "output_stage": "links",
        "reference": _reference_config(genome, reference_root, path_base),
        "samples": [],
        "catalog_import": catalog_import,
        "provenance": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "input_mode": "catalog_artifact_manifest",
            "sample_sheet": _display_path(sample_sheet_path, path_base),
            "sample_sheet_sha256": sample_sheet_digest,
            "sample_sheet_schema": _display_path(schema_path, path_base),
            "sample_sheet_schema_sha256": sha256_file(schema_path),
            "catalog_manifest": _display_path(catalog_manifest_path, path_base),
            "catalog_manifest_sha256": sha256_file(catalog_manifest_path),
            "source_catalog_semantic_sha256": catalog[
                "source_semantic_sha256"
            ],
        },
    }
    config["nearest_tss_links"] = _nearest_tss_link_config(
        genome=genome,
        reference=config["reference"],
        path_base=path_base,
    )
    config["provenance"]["semantic_sha256"] = workflow_semantic_sha256(config)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{project}.catalog-to-links.yaml"
    _write_config_if_changed(output_path, config)
    return output_path.resolve()


def generate_resume_config(
    *,
    checkpoint_manifest_path: Path,
    start_stage: str,
    output_stage: str,
    sample_sheet_path: Path,
    output_dir: Path,
    path_base: Path,
    alignment_chunk_pairs: int | None = None,
) -> Path:
    """Resume a completed boundary in its original result namespace."""

    output_stage = validate_stage_selection(start_stage, output_stage)
    checkpoint = read_stage_checkpoint(
        checkpoint_manifest_path,
        expected_stage=start_stage,
        require_files=True,
    )
    record = checkpoint["artifacts"].get("resolved_config")
    if record is None:
        raise AcquisitionError("Stage checkpoint lacks the resolved_config artifact")
    try:
        source = json.loads(Path(record["path"]).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(f"Cannot read checkpoint resolved config: {error}") from error
    recorded_semantic = source.get("provenance", {}).get("semantic_sha256")
    if recorded_semantic != checkpoint["semantic_sha256"]:
        raise AcquisitionError(
            "Stage checkpoint and resolved configuration semantic digests differ"
        )
    if sha256_file(sample_sheet_path) != source.get("provenance", {}).get(
        "sample_sheet_sha256"
    ):
        raise AcquisitionError(
            "The supplied sample sheet differs from the checkpoint sample sheet"
        )
    resumed = copy.deepcopy(source)
    if alignment_chunk_pairs is not None:
        if alignment_chunk_pairs < 1:
            raise AcquisitionError("Alignment chunk pairs must be positive")
        if start_stage != "trimming" or output_stage not in {"alignment", "qc"}:
            raise AcquisitionError(
                "Alignment chunking requires a trimming checkpoint continuing "
                "through alignment or QC"
            )
        resumed["execution"] = {
            **resumed.get("execution", {}),
            "alignment_chunk_pairs": alignment_chunk_pairs,
        }

    def display(value: str) -> str:
        return _display_path(value, path_base)

    reference = resumed["reference"]
    for key in (
        "fasta",
        "bowtie2_index",
        "chrom_sizes",
        "blacklist_bed",
        "tss_bed",
        "autosomes_file",
    ):
        reference[key] = display(reference[key])
    for sample in resumed["samples"]:
        for key in ("r1", "r2"):
            if key in sample:
                sample[key] = [display(value) for value in sample[key]]
        adapter = sample["parameters"]["trimming"].get("adapter_fasta")
        if adapter:
            sample["parameters"]["trimming"]["adapter_fasta"] = display(adapter)
        if "final_bam" in sample:
            for key in ("bam", "bai"):
                sample["final_bam"][key] = display(sample["final_bam"][key])
        if "qc_peak" in sample:
            sample["qc_peak"]["path"] = display(sample["qc_peak"]["path"])
    for master_key in ("external_master",):
        master = resumed.get(master_key)
        if master:
            for key in MASTER_FILE_FIELDS:
                master[key] = display(master[key])
    activity = resumed.get("activity")
    if activity:
        for library in activity["libraries"]:
            for key in ("bam", "bai"):
                library[key] = display(library[key])
        for key in MASTER_FILE_FIELDS:
            activity["master"][key] = display(activity["master"][key])
    catalog_import = resumed.get("catalog_import")
    if catalog_import:
        for key in ("manifest", *CATALOG_FILE_FIELDS):
            catalog_import[key] = display(catalog_import[key])
    nearest_tss_links = resumed.get("nearest_tss_links")
    if nearest_tss_links:
        nearest_tss_links["promoter_annotation"] = display(
            nearest_tss_links["promoter_annotation"]
        )
    report = resumed.get("report")
    if report:
        report["source_roots"] = [display(value) for value in report["source_roots"]]
        for report_source in report["source_files"]:
            report_source["path"] = display(report_source["path"])
            report_source["source_root"] = display(report_source["source_root"])

    resumed["start_stage"] = start_stage
    resumed["output_stage"] = output_stage
    if workflow_semantic_sha256(resumed) != recorded_semantic:
        raise AcquisitionError(
            "Checkpoint resolved paths cannot be reconstructed into the original "
            "scientific configuration"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (
        f"{resumed['project']}.{start_stage}-to-{output_stage}.yaml"
    )
    _write_config_if_changed(output_path, resumed)
    return output_path.resolve()
