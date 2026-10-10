"""Semantic validation for resolved workflow YAML files."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

from .accessions import AcquisitionError
from .artifacts import (
    CATALOG_FILE_FIELDS,
    FINAL_BAM_FILTERING_CONTRACT,
    semantic_sha256,
    sha256_file,
)
from .stage_checkpoints import LOGICAL_STAGES


SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
CHECKSUM_RE = re.compile(r"^(?:md5:[0-9a-f]{32}|sha256:[0-9a-f]{64})$")
ATAC_QPOIS_FIELDS = {
    "fragment_maximum",
    "minimum_exponent",
    "maximum_exponent",
    "merge_gap",
    "minimum_length",
    "maximum_length",
}
ATAC_CONSENSUS_FIELDS = {
    "enabled",
    "conditions",
    "minimum_replicates",
    "replicate_overlap_fraction",
}
ATAC_MASTER_FIELDS = {"summit_max_distance", "minimum_summit_separation"}
FINAL_BAM_FIELDS = {
    "bam",
    "bai",
    "bam_sha256",
    "bai_sha256",
    "filtering_contract",
    "qc_status",
    "genome",
}
EXTERNAL_MASTER_FIELDS = {
    "genome",
    "method",
    "input_filtering_contract",
    "source_project",
    "source_run_id",
    "master_bed",
    "master_bed_sha256",
    "summits_bed",
    "summits_bed_sha256",
    "membership_tsv",
    "membership_tsv_sha256",
    "context_matrix_tsv",
    "context_matrix_tsv_sha256",
    "stats_json",
    "stats_json_sha256",
}
ACTIVITY_FIELDS = {
    "schema_version",
    "master",
    "contexts",
    "libraries",
    "atac_fragment_maximum",
    "atac_browser_extension_bp",
    "normalization",
    "h3k27ac_signal",
    "h3k27ac_activity_percentile_threshold",
}
ACTIVITY_LIBRARY_FIELDS = {
    "id",
    "assay",
    "context",
    "cohort",
    "layout",
    "genome",
    "bam",
    "bai",
    "bam_sha256",
    "bai_sha256",
    "filtering_contract",
    "qc_status",
}
CATALOG_IMPORT_FIELDS = {
    "schema_version",
    "manifest",
    "manifest_sha256",
    "genome",
    "method",
    "contexts",
    "source_project",
    "source_run_id",
    "source_semantic_sha256",
    "source_sample_sheet_sha256",
    *CATALOG_FILE_FIELDS,
    *(f"{field}_sha256" for field in CATALOG_FILE_FIELDS),
}
NEAREST_TSS_LINK_FIELDS = {
    "schema_version",
    "promoter_annotation",
    "promoter_annotation_checksum",
    "chromosome_scope",
    "promoter_width_bp",
    "enhancer_classes",
    "promoter_activity",
    "candidate_assignment",
}
REPORT_FIELDS = {"schema_version", "source_roots", "source_files"}
REPORT_SOURCE_FIELDS = {"path", "sha256", "kind", "source_root"}
REFERENCE_FIELDS = {
    "name",
    "fasta",
    "bowtie2_index",
    "chrom_sizes",
    "blacklist_bed",
    "tss_bed",
    "autosomes_file",
    "mitochondrial_contig",
    "effective_genome_size",
    "macs3_genome_size",
}
OUTPUT_STAGES = LOGICAL_STAGES
OUTPUT_STAGES_BY_INPUT = {
    "accessions": {"trimming", "alignment", "qc"},
    "trimming": {"trimming", "alignment", "qc"},
    "alignment": {"alignment", "qc"},
    "qc": {"qc", "master"},
    "master": {"master", "quantification", "catalog", "links", "report"},
    "quantification": {"quantification", "catalog", "links", "report"},
    "catalog": {"catalog", "links", "report"},
    "links": {"links", "report"},
    "report": {"report"},
    # Backward-compatible artifact-mode name. New commands use alignment or qc.
    "final-bam": {"alignment", "qc", "master"},
}
DEFAULT_OUTPUT_STAGE_BY_INPUT = {
    "accessions": "qc",
    "trimming": "qc",
    "alignment": "qc",
    "qc": "master",
    "master": "master",
    "quantification": "report",
    "catalog": "report",
    "links": "report",
    "report": "report",
    "final-bam": "master",
}


def wildcard_regex(values: list[str]) -> str:
    return "(?:" + "|".join(re.escape(value) for value in values) + ")" if values else r"(?!)"


def named_artifacts(paths: list[str | None], prefix: str) -> dict[str, str]:
    """Give a deterministic label to each distinct non-empty artifact path."""

    return {
        f"{prefix}.{index:04d}": path
        for index, path in enumerate(dict.fromkeys(paths), start=1)
        if path
    }


def aria2_checksum(source: dict[str, Any]) -> str:
    """Translate a catalog checksum into aria2's ALGORITHM=DIGEST form."""

    algorithm, digest = str(source["checksum"]).split(":", 1)
    return f"{'sha-256' if algorithm == 'sha256' else algorithm}={digest}"


def workflow_semantic_sha256(config: dict[str, Any]) -> str:
    """Return the timestamp-independent digest used to identify a scientific run."""

    semantic_input = {
        key: value
        for key, value in config.items()
        if key
        not in {
            "execution",
            "provenance",
            "start_stage",
            "output_stage",
            "report",
        }
    }
    provenance = config.get("provenance", {})
    semantic_input["provenance_inputs"] = {
        key: value
        for key, value in provenance.items()
        if key not in {"generated_at_utc", "semantic_sha256"}
    }
    return semantic_sha256(semantic_input)


def validate_stage_selection(input_stage: str, output_stage: str | None) -> str:
    """Resolve and validate a non-scientific workflow stopping point."""

    if input_stage not in OUTPUT_STAGES_BY_INPUT:
        raise AcquisitionError(f"Unsupported input stage: {input_stage!r}")
    resolved = output_stage or DEFAULT_OUTPUT_STAGE_BY_INPUT[input_stage]
    if resolved not in OUTPUT_STAGES:
        raise AcquisitionError(f"Unsupported output stage: {resolved!r}")
    if resolved not in OUTPUT_STAGES_BY_INPUT[input_stage]:
        allowed = ", ".join(
            stage for stage in OUTPUT_STAGES if stage in OUTPUT_STAGES_BY_INPUT[input_stage]
        )
        raise AcquisitionError(
            f"Cannot stop at {resolved!r} when starting from {input_stage!r}; "
            f"allowed output stages: {allowed}"
        )
    return resolved


def guard_result_namespace(config: dict[str, Any], resolved_config: Path) -> None:
    """Refuse to reuse a result namespace for a different scientific config."""

    if not resolved_config.is_file():
        return
    requested = config.get("provenance", {}).get("semantic_sha256")
    # Hand-authored legacy configs are outside the generated-config contract.
    # Canonical configs from run_pipeline.py always carry this digest.
    if not requested:
        return
    try:
        existing = json.loads(resolved_config.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(
            f"Cannot verify existing result namespace {resolved_config}: {error}. "
            "Use a new run_id."
        ) from error
    recorded = existing.get("provenance", {}).get("semantic_sha256")
    if not recorded:
        raise AcquisitionError(
            f"Cannot prove that existing result namespace {resolved_config.parent.parent} "
            "has the same scientific configuration. Use a new run_id."
        )
    if requested != recorded:
        raise AcquisitionError(
            f"Result namespace {resolved_config.parent.parent} belongs to semantic "
            f"configuration {recorded}, not {requested}. Use a new run_id."
        )


def _required(mapping: dict[str, Any], fields: set[str], label: str) -> None:
    missing = sorted(field for field in fields if field not in mapping or mapping[field] is None)
    if missing:
        raise AcquisitionError(f"{label} is missing: " + ", ".join(missing))


def validate_workflow_config(config: dict[str, Any]) -> None:
    """Reject inconsistent resolved configurations before Snakemake builds the DAG."""

    _required(config, {"project", "run_id", "output_dir", "assay", "reference", "samples"}, "Config")
    for field in ("project", "run_id"):
        if not SAFE_ID_RE.fullmatch(str(config[field])):
            raise AcquisitionError(f"Invalid {field}: {config[field]!r}")
    if config["assay"] not in {"atac", "chip_histone", "activity"}:
        raise AcquisitionError(f"Unsupported assay: {config['assay']!r}")
    if not isinstance(config["reference"], dict):
        raise AcquisitionError("reference must be a mapping")
    _required(config["reference"], REFERENCE_FIELDS, "Reference")
    input_stage = str(config.get("input_stage", "accessions"))
    if input_stage not in {
        "accessions",
        "final-bam",
        "master",
        "quantification",
        "catalog",
    }:
        raise AcquisitionError(f"Unsupported input_stage: {input_stage!r}")
    start_stage = str(config.get("start_stage", input_stage))
    output_stage = validate_stage_selection(start_stage, config.get("output_stage"))
    valid_input_modes = {
        "accessions": {"accessions", "trimming", "alignment"},
        "final-bam": {"alignment", "qc", "final-bam"},
        "master": {"master"},
        "quantification": {
            "master",
            "quantification",
            "catalog",
            "links",
            "report",
        },
        "catalog": {"catalog", "links"},
    }
    if start_stage not in valid_input_modes[input_stage]:
        raise AcquisitionError(
            f"start_stage {start_stage!r} is incompatible with input_stage "
            f"{input_stage!r}"
        )
    provenance = config.get("provenance")
    if provenance is not None:
        if not isinstance(provenance, dict):
            raise AcquisitionError("provenance must be a mapping")
        recorded_semantic = provenance.get("semantic_sha256")
        if recorded_semantic is not None:
            if not re.fullmatch(r"[0-9a-f]{64}", str(recorded_semantic)):
                raise AcquisitionError(
                    "provenance semantic_sha256 is not a SHA-256 digest"
                )
            calculated_semantic = workflow_semantic_sha256(config)
            if recorded_semantic != calculated_semantic:
                raise AcquisitionError(
                    "provenance semantic_sha256 does not match the resolved configuration"
                )
    execution = config.get("execution")
    if execution is not None:
        if not isinstance(execution, dict):
            raise AcquisitionError("execution must be a mapping")
        unknown_execution = sorted(set(execution) - {"alignment_chunk_pairs"})
        if unknown_execution:
            raise AcquisitionError(
                "Unsupported execution option(s): " + ", ".join(unknown_execution)
            )
        chunk_pairs = execution.get("alignment_chunk_pairs")
        if (
            chunk_pairs is not None
            and (
                isinstance(chunk_pairs, bool)
                or not isinstance(chunk_pairs, int)
                or chunk_pairs < 1
            )
        ):
            raise AcquisitionError(
                "execution alignment_chunk_pairs must be a positive integer"
            )
    external_master = config.get("external_master")
    if input_stage == "master":
        if config["assay"] != "atac":
            raise AcquisitionError("Master reuse is valid only for ATAC")
        if not isinstance(external_master, dict):
            raise AcquisitionError("Master reuse requires external_master")
        _required(external_master, EXTERNAL_MASTER_FIELDS, "External master")
        if external_master["genome"] != config["reference"].get("name"):
            raise AcquisitionError("External master and reference genomes differ")
        if external_master["input_filtering_contract"] != FINAL_BAM_FILTERING_CONTRACT:
            raise AcquisitionError("External master has an incompatible BAM filtering contract")
        for field in EXTERNAL_MASTER_FIELDS:
            if field.endswith("_sha256") and not re.fullmatch(
                r"[0-9a-f]{64}", str(external_master[field])
            ):
                raise AcquisitionError(f"External master {field} is not a SHA-256 digest")
    elif external_master is not None:
        raise AcquisitionError("external_master requires input_stage=master")
    activity = config.get("activity")
    catalog_import = config.get("catalog_import")
    if input_stage == "quantification":
        if config["assay"] != "activity" or not isinstance(activity, dict):
            raise AcquisitionError(
                "Quantification input stage requires assay=activity and an activity mapping"
            )
        _required(activity, ACTIVITY_FIELDS, "Activity")
        master = activity["master"]
        if not isinstance(master, dict):
            raise AcquisitionError("Activity master must be a mapping")
        _required(master, EXTERNAL_MASTER_FIELDS, "Activity master")
        if master["genome"] != config["reference"]["name"]:
            raise AcquisitionError("Activity master and reference genomes differ")
        if master["input_filtering_contract"] != FINAL_BAM_FILTERING_CONTRACT:
            raise AcquisitionError("Activity master has an incompatible BAM filtering contract")
        for field in EXTERNAL_MASTER_FIELDS:
            if field.endswith("_sha256") and not re.fullmatch(
                r"[0-9a-f]{64}", str(master[field])
            ):
                raise AcquisitionError(
                    f"Activity master {field} is not a SHA-256 digest"
                )
        if activity["schema_version"] != 6:
            raise AcquisitionError("Unsupported activity schema version")
        if int(activity["atac_fragment_maximum"]) < 2:
            raise AcquisitionError("Activity ATAC fragment maximum is invalid")
        if int(activity["atac_browser_extension_bp"]) != 150:
            raise AcquisitionError("Activity ATAC browser extension must be 150 bp")
        if activity["normalization"] != "background_tmm_10kb_autosomes_v1":
            raise AcquisitionError("Unsupported activity normalization method")
        if activity["h3k27ac_signal"] != "summit_max3_background_tmm_v3":
            raise AcquisitionError("Unsupported regulatory-element method")
        threshold = activity["h3k27ac_activity_percentile_threshold"]
        if (
            isinstance(threshold, bool)
            or not isinstance(threshold, (int, float))
            or not math.isfinite(float(threshold))
            or not 0 < float(threshold) < 1
        ):
            raise AcquisitionError(
                "H3K27ac activity percentile threshold must be in (0, 1)"
            )
        contexts = activity["contexts"]
        if (
            not isinstance(contexts, list)
            or not contexts
            or len(contexts) != len(set(contexts))
            or any(not SAFE_ID_RE.fullmatch(str(item)) for item in contexts)
        ):
            raise AcquisitionError("Activity contexts are invalid")
        libraries = activity["libraries"]
        if not isinstance(libraries, list) or not libraries:
            raise AcquisitionError("Activity libraries must be a non-empty list")
        library_ids = [str(library.get("id", "")) for library in libraries]
        if (
            any(not SAFE_ID_RE.fullmatch(item) for item in library_ids)
            or len(library_ids) != len(set(library_ids))
        ):
            raise AcquisitionError("Activity library IDs must be safe and unique")
        artifact_paths = [
            str(library.get(field, ""))
            for library in libraries
            for field in ("bam", "bai")
        ]
        if (
            any(not item for item in artifact_paths)
            or len(artifact_paths) != len(set(artifact_paths))
        ):
            raise AcquisitionError(
                "Activity libraries must use distinct BAM and BAI artifact paths"
            )
        for library in libraries:
            _required(
                library,
                ACTIVITY_LIBRARY_FIELDS,
                f"Activity library {library.get('id', '')}",
            )
            library_id = str(library["id"])
            if library["assay"] not in {"atac", "h3k27ac"}:
                raise AcquisitionError(
                    f"Activity library {library_id}: invalid assay"
                )
            if library["cohort"] != "atlas":
                raise AcquisitionError(
                    f"Activity library {library_id}: cohort must be atlas"
                )
            layout = library["layout"]
            if layout not in {"single", "paired"}:
                raise AcquisitionError(
                    f"Activity library {library_id}: invalid layout"
                )
            estimated_fragment_length = library.get(
                "estimated_fragment_length_bp"
            )
            if library["assay"] == "atac" and layout != "paired":
                raise AcquisitionError(
                    f"Activity library {library_id}: ATAC must be paired-end"
                )
            if library["assay"] == "h3k27ac" and layout == "single":
                if (
                    not isinstance(estimated_fragment_length, int)
                    or isinstance(estimated_fragment_length, bool)
                    or estimated_fragment_length < 1
                ):
                    raise AcquisitionError(
                        f"Activity library {library_id}: single-end H3K27ac "
                        "requires a positive estimated_fragment_length_bp"
                    )
            elif (
                estimated_fragment_length is not None
                and estimated_fragment_length != ""
            ):
                raise AcquisitionError(
                    f"Activity library {library_id}: unexpected "
                    "estimated_fragment_length_bp"
                )
            if library["genome"] != config["reference"]["name"]:
                raise AcquisitionError(
                    f"Activity library {library_id}: genome differs from reference"
                )
            if library["filtering_contract"] != FINAL_BAM_FILTERING_CONTRACT:
                raise AcquisitionError(
                    f"Activity library {library_id}: invalid filtering contract"
                )
            if library["qc_status"] != "accepted":
                raise AcquisitionError(
                    f"Activity library {library_id}: accepted QC is required"
                )
            for field in ("bam_sha256", "bai_sha256"):
                if not re.fullmatch(r"[0-9a-f]{64}", str(library[field])):
                    raise AcquisitionError(
                        f"Activity library {library_id}: invalid {field}"
                    )
            context = str(library["context"])
            if not SAFE_ID_RE.fullmatch(context):
                raise AcquisitionError(
                    f"Activity library {library_id}: invalid context"
                )
            if context not in contexts:
                raise AcquisitionError(
                    f"Activity library {library_id}: unknown context"
                )
        for context in contexts:
            assays = {
                library["assay"]
                for library in libraries
                if library["context"] == context
            }
            if assays != {"atac", "h3k27ac"}:
                raise AcquisitionError(
                    f"Activity context {context!r} lacks an assay"
                )
        for assay in ("atac", "h3k27ac"):
            count = sum(
                library["assay"] == assay
                for library in libraries
            )
            if count < 2:
                raise AcquisitionError(
                    f"Activity background TMM requires two {assay} libraries"
                )
        report = config.get("report")
        if not isinstance(report, dict):
            raise AcquisitionError("Quantification configuration requires report metadata")
        _required(report, REPORT_FIELDS, "Report")
        if report["schema_version"] != 1:
            raise AcquisitionError("Unsupported report schema version")
        roots = report["source_roots"]
        if (
            not isinstance(roots, list)
            or len(roots) != len(set(str(item) for item in roots))
        ):
            raise AcquisitionError("Report source roots must be a unique list")
        source_files = report["source_files"]
        if not isinstance(source_files, list) or not source_files:
            raise AcquisitionError("Report source files must be a non-empty list")
        source_paths = []
        for source in source_files:
            if not isinstance(source, dict):
                raise AcquisitionError("Report source-file records must be mappings")
            _required(source, REPORT_SOURCE_FIELDS, "Report source file")
            if not re.fullmatch(r"[0-9a-f]{64}", str(source["sha256"])):
                raise AcquisitionError("Report source-file SHA-256 is invalid")
            if not source["path"] or not source["kind"]:
                raise AcquisitionError("Report source-file path and kind are required")
            source_paths.append(str(source["path"]))
        if len(source_paths) != len(set(source_paths)):
            raise AcquisitionError("Report source-file paths must be unique")
    elif activity is not None:
        raise AcquisitionError("activity mapping requires input_stage=quantification")

    if input_stage == "catalog":
        if config["assay"] != "activity" or not isinstance(catalog_import, dict):
            raise AcquisitionError(
                "Catalog input stage requires assay=activity and catalog_import"
            )
        if output_stage != "links":
            raise AcquisitionError("Catalog import currently supports only links output")
        _required(catalog_import, CATALOG_IMPORT_FIELDS, "Catalog import")
        if catalog_import["schema_version"] != 1:
            raise AcquisitionError("Unsupported catalog-import schema version")
        if catalog_import["genome"] != config["reference"]["name"]:
            raise AcquisitionError("Imported catalog and reference genomes differ")
        for field in ("method", "source_project", "source_run_id"):
            if not SAFE_ID_RE.fullmatch(str(catalog_import[field])):
                raise AcquisitionError(f"Catalog import {field} is invalid")
        catalog_contexts = catalog_import["contexts"]
        if (
            not isinstance(catalog_contexts, list)
            or not catalog_contexts
            or len(catalog_contexts) != len(set(catalog_contexts))
            or any(
                not SAFE_ID_RE.fullmatch(str(context))
                for context in catalog_contexts
            )
        ):
            raise AcquisitionError("Catalog import contexts are invalid")
        digest_fields = {
            "manifest_sha256",
            "source_semantic_sha256",
            "source_sample_sheet_sha256",
            *(f"{field}_sha256" for field in CATALOG_FILE_FIELDS),
        }
        for field in digest_fields:
            if not re.fullmatch(r"[0-9a-f]{64}", str(catalog_import[field])):
                raise AcquisitionError(f"Catalog import {field} is not a SHA-256 digest")
        for field in ("manifest", *CATALOG_FILE_FIELDS):
            if not catalog_import[field]:
                raise AcquisitionError(f"Catalog import {field} path is blank")
            artifact = Path(str(catalog_import[field]))
            if artifact.is_file():
                digest_field = (
                    "manifest_sha256" if field == "manifest" else f"{field}_sha256"
                )
                if sha256_file(artifact) != catalog_import[digest_field]:
                    raise AcquisitionError(f"Catalog import {field} SHA-256 mismatch")
    elif catalog_import is not None:
        raise AcquisitionError("catalog_import requires input_stage=catalog")

    nearest_tss_links = config.get("nearest_tss_links")
    if input_stage in {"quantification", "catalog"}:
        if not isinstance(nearest_tss_links, dict):
            raise AcquisitionError(
                "Activity catalog/links configurations require nearest_tss_links"
            )
        _required(
            nearest_tss_links,
            NEAREST_TSS_LINK_FIELDS,
            "Nearest-TSS links",
        )
        if nearest_tss_links["schema_version"] != 2:
            raise AcquisitionError("Unsupported nearest-TSS links schema version")
        if not nearest_tss_links["promoter_annotation"] or not CHECKSUM_RE.fullmatch(
            str(nearest_tss_links["promoter_annotation_checksum"])
        ):
            raise AcquisitionError(
                "Nearest-TSS promoter annotation path/checksum is invalid"
            )
        if (
            nearest_tss_links["chromosome_scope"]
            != "all_annotation_chromosomes_in_chrom_sizes"
        ):
            raise AcquisitionError("Unsupported nearest-TSS chromosome scope")
        if int(nearest_tss_links["promoter_width_bp"]) != 1000:
            raise AcquisitionError(
                "Nearest-TSS promoter display width must be 1000 bp"
            )
        if nearest_tss_links["enhancer_classes"] != [
            "proximal_enhancer_like",
            "distal_enhancer_like",
        ]:
            raise AcquisitionError("Unsupported nearest-TSS enhancer classes")
        expected_nearest_methods = {
            "promoter_activity": "summit_within_500bp_catalog_activity_v3",
            "candidate_assignment": "nearest_annotated_tss_all_ties_v1",
        }
        for field, expected in expected_nearest_methods.items():
            if nearest_tss_links[field] != expected:
                raise AcquisitionError(
                    f"Unsupported nearest-TSS {field}: {nearest_tss_links[field]!r}"
                )
    elif nearest_tss_links is not None:
        raise AcquisitionError(
            "nearest_tss_links requires an activity catalog/links configuration"
        )

    if "contacts" in config:
        raise AcquisitionError(
            "The contact-link branch is not part of this workflow; remove the "
            "contacts mapping"
        )
    qpois = config.get("atac_qpois")
    if (
        config["assay"] == "atac"
        and input_stage != "master"
        and not isinstance(qpois, dict)
    ):
        raise AcquisitionError("ATAC configurations require atac_qpois parameters")
    if qpois is not None:
        if config["assay"] != "atac" or not isinstance(qpois, dict):
            raise AcquisitionError("atac_qpois is only valid for ATAC configurations")
        _required(qpois, ATAC_QPOIS_FIELDS, "ATAC qpois")
        if (
            int(qpois["fragment_maximum"]) < 2
            or int(qpois["minimum_exponent"]) < 0
            or int(qpois["maximum_exponent"]) < int(qpois["minimum_exponent"])
            or int(qpois["merge_gap"]) < 0
            or int(qpois["minimum_length"]) < 1
            or int(qpois["maximum_length"]) < int(qpois["minimum_length"])
        ):
            raise AcquisitionError("ATAC qpois parameters are invalid")
    consensus = config.get("atac_consensus")
    if consensus is not None:
        if config["assay"] != "atac" or not isinstance(consensus, dict):
            raise AcquisitionError("atac_consensus is only valid for ATAC configurations")
        _required(consensus, ATAC_CONSENSUS_FIELDS, "ATAC consensus")
        if not isinstance(consensus["enabled"], bool):
            raise AcquisitionError("ATAC consensus enabled must be true or false")
        if (
            int(consensus["minimum_replicates"]) < 2
            or not 0 < float(consensus["replicate_overlap_fraction"]) <= 1
            or not isinstance(consensus["conditions"], list)
            or not consensus["conditions"]
        ):
            raise AcquisitionError("ATAC consensus parameters are invalid")
    master = config.get("atac_master")
    if master is not None:
        if config["assay"] != "atac" or not isinstance(master, dict):
            raise AcquisitionError("atac_master is only valid for ATAC configurations")
        _required(master, ATAC_MASTER_FIELDS, "ATAC master DHS")
        if (
            int(master["summit_max_distance"]) < 0
            or int(master["minimum_summit_separation"]) < 0
            or int(master["minimum_summit_separation"])
            > int(master["summit_max_distance"])
        ):
            raise AcquisitionError("ATAC master DHS parameters are invalid")
        if not consensus or not consensus.get("enabled"):
            raise AcquisitionError("ATAC master DHS construction requires ATAC consensus")
    if int(config["reference"]["effective_genome_size"]) < 1:
        raise AcquisitionError("effective_genome_size must be positive")
    preparation = config["reference"].get("preparation")
    if preparation is not None:
        if not isinstance(preparation, dict) or preparation.get("mode") != "download":
            raise AcquisitionError("Reference preparation mode must be 'download'")
        _required(
            preparation,
            {"mode", "fasta", "annotation", "blacklist", "autosomes"},
            "Reference preparation",
        )
        for source_name in ("fasta", "annotation", "blacklist"):
            source = preparation[source_name]
            if not isinstance(source, dict):
                raise AcquisitionError(f"Reference {source_name} source must be a mapping")
            _required(source, {"url", "checksum"}, f"Reference {source_name} source")
            if not str(source["url"]).startswith("https://"):
                raise AcquisitionError(f"Reference {source_name} URL must use HTTPS")
            if not CHECKSUM_RE.fullmatch(str(source["checksum"])):
                raise AcquisitionError(
                    f"Reference {source_name} checksum must be md5:<hex> or sha256:<hex>"
                )
        autosomes = preparation["autosomes"]
        if not isinstance(autosomes, list) or not autosomes or any(not item for item in autosomes):
            raise AcquisitionError("Reference autosomes must be a non-empty list")

    samples = config["samples"]
    if not isinstance(samples, list) or (
        not samples and input_stage not in {"master", "quantification", "catalog"}
    ):
        raise AcquisitionError("samples must be a non-empty list")
    if input_stage in {"master", "quantification", "catalog"} and samples:
        raise AcquisitionError(
            f"{input_stage.capitalize()} configuration must not schedule sample processing"
        )
    sample_ids = [str(sample.get("id", "")) for sample in samples]
    if any(not SAFE_ID_RE.fullmatch(sample_id) for sample_id in sample_ids):
        raise AcquisitionError("Every sample must have a safe non-empty id")
    if len(sample_ids) != len(set(sample_ids)):
        raise AcquisitionError("Sample IDs must be unique")
    context_by_id = {str(sample["id"]): str(sample.get("context", "")) for sample in samples}

    for sample in samples:
        sample_id = str(sample["id"])
        sample_fields = {"accessions", "context", "role", "layout", "parameters"}
        if input_stage == "accessions":
            sample_fields.add("r1")
        else:
            sample_fields.add("final_bam")
        _required(sample, sample_fields, f"Sample {sample_id}")
        if not SAFE_ID_RE.fullmatch(str(sample["context"])):
            raise AcquisitionError(f"Sample {sample_id}: invalid context")
        if sample["role"] != "treatment":
            raise AcquisitionError(
                f"Sample {sample_id}: only IP/treatment libraries are supported"
            )
        if sample["layout"] not in {"single", "paired"}:
            raise AcquisitionError(f"Sample {sample_id}: invalid layout")
        if input_stage == "accessions":
            if not isinstance(sample["r1"], list) or not sample["r1"]:
                raise AcquisitionError(
                    f"Sample {sample_id}: r1 must contain at least one FASTQ"
                )
            if sample["layout"] == "paired":
                if not isinstance(sample.get("r2"), list) or len(sample["r1"]) != len(
                    sample["r2"]
                ):
                    raise AcquisitionError(
                        f"Sample {sample_id}: paired r1/r2 lane counts differ"
                    )
            elif sample.get("r2"):
                raise AcquisitionError(
                    f"Sample {sample_id}: single-end input must not contain r2"
                )
        else:
            if sample.get("r1") or sample.get("r2"):
                raise AcquisitionError(
                    f"Sample {sample_id}: final-BAM mode must not contain FASTQs"
                )
            final_bam = sample["final_bam"]
            if not isinstance(final_bam, dict):
                raise AcquisitionError(
                    f"Sample {sample_id}: final_bam must be a mapping"
                )
            _required(final_bam, FINAL_BAM_FIELDS, f"Sample {sample_id} final BAM")
            if final_bam["genome"] != config["reference"]["name"]:
                raise AcquisitionError(
                    f"Sample {sample_id}: final BAM and reference genomes differ"
                )
            allowed_qc_statuses = (
                {"accepted"}
                if output_stage == "master"
                else {"pending_review", "accepted"}
            )
            if final_bam["qc_status"] not in allowed_qc_statuses:
                raise AcquisitionError(
                    f"Sample {sample_id}: final BAM QC status "
                    f"{final_bam['qc_status']!r} is invalid for output stage "
                    f"{output_stage!r}"
                )
            if final_bam["filtering_contract"] != FINAL_BAM_FILTERING_CONTRACT:
                raise AcquisitionError(
                    f"Sample {sample_id}: unsupported final BAM filtering contract"
                )
            for field in ("bam_sha256", "bai_sha256"):
                if not re.fullmatch(r"[0-9a-f]{64}", str(final_bam[field])):
                    raise AcquisitionError(
                        f"Sample {sample_id}: {field} is not a SHA-256 digest"
                    )

        if sample["role"] == "treatment":
            peak = sample.get("peak_caller")
            if not isinstance(peak, dict) or peak.get("command") != "callpeak":
                raise AcquisitionError(
                    f"Sample {sample_id}: treatment peak caller must use callpeak"
                )
            _required(
                peak,
                {"format", "qvalue", "broad", "nomodel", "write_bedgraph", "spmr"},
                f"Sample {sample_id} callpeak",
            )
            if config["assay"] == "atac":
                if (
                    peak.get("mode") != "tn5_qpois"
                    or peak["format"] != "BED"
                    or peak["broad"]
                    or not peak["nomodel"]
                    or not peak["write_bedgraph"]
                    or peak["spmr"]
                    or peak.get("shift") is None
                    or peak.get("extsize") is None
                ):
                    raise AcquisitionError(
                        f"Sample {sample_id}: invalid two-ended Tn5 qpois configuration"
                    )
            elif not peak["broad"]:
                raise AcquisitionError(
                    f"Sample {sample_id}: H3K27ac callpeak must use broad mode"
                )
            elif not peak["write_bedgraph"] or not peak["spmr"]:
                raise AcquisitionError(
                    f"Sample {sample_id}: H3K27ac callpeak must write -B --SPMR bedGraphs"
                )
            qc_peak = sample.get("qc_peak")
            if qc_peak is not None:
                if start_stage != "qc" or config["assay"] != "atac":
                    raise AcquisitionError(
                        f"Sample {sample_id}: qc_peak is valid only for ATAC QC reuse"
                    )
                if not isinstance(qc_peak, dict):
                    raise AcquisitionError(
                        f"Sample {sample_id}: qc_peak must be a mapping"
                    )
                _required(
                    qc_peak,
                    {"path", "sha256", "method"},
                    f"Sample {sample_id} QC peak",
                )
                if qc_peak["method"] != peak["command"]:
                    raise AcquisitionError(
                        f"Sample {sample_id}: QC peak method differs from peak caller"
                    )
                if not re.fullmatch(r"[0-9a-f]{64}", str(qc_peak["sha256"])):
                    raise AcquisitionError(
                        f"Sample {sample_id}: QC peak SHA-256 is invalid"
                    )

        if "control" in sample:
            raise AcquisitionError(
                f"Sample {sample_id}: input controls are not supported"
            )

    if consensus and consensus["enabled"] and config["assay"] != "atac":
        raise AcquisitionError("ATAC consensus requires an ATAC configuration")
    if consensus and consensus["enabled"]:
        from .consensus import condition_specs

        specifications = condition_specs(
            consensus["conditions"],
            sample_ids=sample_ids,
            minimum_replicates=int(consensus["minimum_replicates"]),
        )
        for specification in specifications:
            contexts = {
                context_by_id[sample_id] for sample_id in specification.samples
            }
            if contexts != {specification.condition_id}:
                raise AcquisitionError(
                    f"ATAC condition {specification.condition_id!r} does not match "
                    "its sample contexts"
                )


def resolve_input_paths(config: dict[str, Any], base: Path) -> None:
    """Resolve relative FASTQ/reference paths against a config's launch directory in place."""

    reference = config["reference"]
    for key in ("fasta", "bowtie2_index", "chrom_sizes", "blacklist_bed", "tss_bed", "autosomes_file"):
        path = Path(reference[key])
        reference[key] = str(path if path.is_absolute() else (base / path).resolve())
    for sample in config["samples"]:
        for key in ("r1", "r2"):
            if key in sample:
                sample[key] = [
                    str(path if path.is_absolute() else (base / path).resolve())
                    for value in sample[key]
                    for path in [Path(value)]
                ]
        adapter = sample["parameters"]["trimming"].get("adapter_fasta")
        if adapter:
            path = Path(adapter)
            sample["parameters"]["trimming"]["adapter_fasta"] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
        final_bam = sample.get("final_bam")
        if final_bam:
            for key in ("bam", "bai"):
                path = Path(final_bam[key])
                final_bam[key] = str(
                    path if path.is_absolute() else (base / path).resolve()
                )
        qc_peak = sample.get("qc_peak")
        if qc_peak:
            path = Path(qc_peak["path"])
            qc_peak["path"] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
    external_master = config.get("external_master")
    if external_master:
        for key in (
            "master_bed",
            "summits_bed",
            "membership_tsv",
            "context_matrix_tsv",
            "stats_json",
        ):
            path = Path(external_master[key])
            external_master[key] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
    activity = config.get("activity")
    if activity:
        for library in activity["libraries"]:
            for key in ("bam", "bai"):
                path = Path(library[key])
                library[key] = str(
                    path if path.is_absolute() else (base / path).resolve()
                )
        for key in (
            "master_bed",
            "summits_bed",
            "membership_tsv",
            "context_matrix_tsv",
            "stats_json",
        ):
            path = Path(activity["master"][key])
            activity["master"][key] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
    catalog_import = config.get("catalog_import")
    if catalog_import:
        for key in ("manifest", *CATALOG_FILE_FIELDS):
            path = Path(catalog_import[key])
            catalog_import[key] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
    nearest_tss_links = config.get("nearest_tss_links")
    if nearest_tss_links:
        path = Path(nearest_tss_links["promoter_annotation"])
        nearest_tss_links["promoter_annotation"] = str(
            path if path.is_absolute() else (base / path).resolve()
        )
    report = config.get("report")
    if report:
        report["source_roots"] = [
            str(path if path.is_absolute() else (base / path).resolve())
            for value in report["source_roots"]
            for path in [Path(value)]
        ]
        for source in report["source_files"]:
            path = Path(source["path"])
            source["path"] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
            root = Path(source["source_root"])
            source["source_root"] = str(
                root if root.is_absolute() else (base / root).resolve()
            )
