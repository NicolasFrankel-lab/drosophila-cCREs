"""Static dm6 contact-source metadata and strict manifest validation."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import re
from typing import Any

from .accessions import AcquisitionError
from .artifacts import SAFE_ID_RE, sha256_file


CONTACT_SOURCE_COLUMNS = (
    "source_id",
    "context",
    "assay",
    "replicate",
    "format",
    "url",
    "local_path",
    "checksum",
    "match_quality",
    "biological_context",
    "caveat",
)
CONTACT_FORMAT_SUFFIXES = {
    "mcool": ".mcool",
    "cool.gz": ".cool.gz",
    "h5": ".h5",
}
CHECKSUM_RE = re.compile(r"^(?:md5:[0-9a-f]{32}|sha256:[0-9a-f]{64})$")

DM6_ATLAS_CONTACT_CONTEXTS: tuple[dict[str, Any], ...] = (
    {
        "id": "ab",
        "strategy": "observed",
        "assay": "Micro-C",
        "match": "tissue_and_stage_matched",
        "resolution_bp": 5000,
        "caveat": "Two adult-CNS biological replicates.",
    },
    {
        "id": "e11",
        "strategy": "observed",
        "assay": "Micro-C",
        "match": "stage_matched_whole_embryo",
        "resolution_bp": 5000,
        "caveat": "Whole-embryo stages 10-12.",
    },
    {
        "id": "e13",
        "strategy": "powerlaw",
        "assay": "distance_model",
        "match": "no_exact_whole_embryo_map",
        "resolution_bp": 5000,
        "caveat": "No defensible exact contact map; distance model only.",
    },
    {
        "id": "e5",
        "strategy": "observed",
        "assay": "Micro-C",
        "match": "stage_matched_whole_embryo",
        "resolution_bp": 5000,
        "caveat": "Nuclear cycle 14 is the closest stage-5 proxy.",
    },
    {
        "id": "ead",
        "strategy": "observed",
        "assay": "Micro-C",
        "match": "tissue_matched_merged",
        "resolution_bp": 5000,
        "caveat": "Eye-disc control map supplied as a four-sample merge.",
    },
    {
        "id": "hid",
        "strategy": "powerlaw",
        "assay": "distance_model",
        "match": "no_exact_haltere_map",
        "resolution_bp": 5000,
        "caveat": "No defensible exact contact map; distance model only.",
    },
    {
        "id": "lb",
        "strategy": "observed",
        "assay": "Micro-C",
        "match": "tissue_and_stage_matched",
        "resolution_bp": 5000,
        "caveat": "Two third-instar larval-CNS biological replicates.",
    },
    {
        "id": "o",
        "strategy": "observed",
        "assay": "Hi-C",
        "match": "whole_ovary_stage_subset",
        "resolution_bp": 4000,
        "caveat": "Hi-C, not Micro-C; egg-chamber stages 1-8 only.",
    },
    {
        "id": "wid",
        "strategy": "observed",
        "assay": "Micro-C",
        "match": "tissue_matched_merged",
        "resolution_bp": 5000,
        "caveat": "Wing-disc WT map supplied as a four-sample merge.",
    },
)
DM6_ATLAS_CONTEXT_IDS = frozenset(row["id"] for row in DM6_ATLAS_CONTACT_CONTEXTS)
DM6_CANONICAL_CONTACT_CHROMOSOMES = (
    "chr2L",
    "chr2R",
    "chr3L",
    "chr3R",
    "chr4",
    "chrX",
    "chrY",
)


def read_contact_source_manifest(path: Path) -> list[dict[str, str]]:
    """Return validated source rows without resolving their local paths."""

    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != CONTACT_SOURCE_COLUMNS:
                raise AcquisitionError(
                    f"Contact manifest {path} must have columns: "
                    + ", ".join(CONTACT_SOURCE_COLUMNS)
                )
            rows = list(reader)
    except OSError as error:
        raise AcquisitionError(f"Cannot read contact manifest {path}: {error}") from error
    if not rows:
        raise AcquisitionError(f"Contact manifest {path} has no data rows")

    seen_ids: set[str] = set()
    seen_paths: set[str] = set()
    replicates_by_context: dict[str, set[str]] = {}
    for line_number, row in enumerate(rows, start=2):
        source_id = row["source_id"]
        context = row["context"]
        if not SAFE_ID_RE.fullmatch(source_id) or not SAFE_ID_RE.fullmatch(context):
            raise AcquisitionError(
                f"Contact manifest {path}:{line_number} has an invalid ID"
            )
        if source_id in seen_ids:
            raise AcquisitionError(f"Duplicate contact source_id {source_id!r}")
        seen_ids.add(source_id)
        file_format = row["format"]
        if file_format not in CONTACT_FORMAT_SUFFIXES:
            raise AcquisitionError(
                f"Contact source {source_id!r} has unsupported format {file_format!r}"
            )
        expected_path = f"data/raw/contacts/{source_id}{CONTACT_FORMAT_SUFFIXES[file_format]}"
        if row["local_path"] != expected_path:
            raise AcquisitionError(
                f"Contact source {source_id!r} must use local_path={expected_path}"
            )
        if row["local_path"] in seen_paths:
            raise AcquisitionError("Contact manifest repeats a local path")
        seen_paths.add(row["local_path"])
        if not row["url"].startswith("https://"):
            raise AcquisitionError(f"Contact source {source_id!r} must use HTTPS")
        checksum = row["checksum"]
        if checksum and not CHECKSUM_RE.fullmatch(checksum):
            raise AcquisitionError(
                f"Contact source {source_id!r} has an invalid checksum"
            )
        if row["assay"] not in {"Micro-C", "Hi-C"}:
            raise AcquisitionError(f"Contact source {source_id!r} has invalid assay")
        if (
            not row["replicate"]
            or not row["match_quality"]
            or not row["biological_context"]
        ):
            raise AcquisitionError(f"Contact source {source_id!r} lacks metadata")
        context_replicates = replicates_by_context.setdefault(context, set())
        if row["replicate"] in context_replicates:
            raise AcquisitionError(
                f"Contact context {context!r} repeats replicate {row['replicate']!r}"
            )
        context_replicates.add(row["replicate"])
    for context in sorted({row["context"] for row in rows}):
        context_rows = [row for row in rows if row["context"] == context]
        for field in ("assay", "match_quality", "biological_context", "caveat"):
            if len({row[field] for row in context_rows}) != 1:
                raise AcquisitionError(
                    f"Contact context {context!r} has inconsistent {field} metadata"
                )
    return rows


def verify_reported_checksum(path: Path, checksum: str) -> None:
    """Validate one source-reported MD5/SHA-256 when a digest is available."""

    if not checksum:
        return
    if not CHECKSUM_RE.fullmatch(checksum):
        raise AcquisitionError(f"Invalid reported checksum for {path}")
    algorithm, expected = checksum.split(":", 1)
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != expected:
        raise AcquisitionError(f"Reported {algorithm} checksum mismatch: {path}")


def default_dm6_contact_config(
    *,
    contexts: list[str],
    reference: dict[str, Any],
    manifest_path: Path,
    path_base: Path,
) -> dict[str, Any] | None:
    """Return the canonical atlas contact configuration when all contexts match."""

    if set(contexts) != DM6_ATLAS_CONTEXT_IDS:
        return None
    annotation_url = str(reference["preparation"]["annotation"]["url"])
    annotation_name = annotation_url.rsplit("/", 1)[-1]
    annotation = Path(str(reference["fasta"])).parent / "sources" / annotation_name

    def display(path: Path) -> str:
        resolved = path.resolve() if path.is_absolute() else (path_base / path).resolve()
        try:
            return str(resolved.relative_to(path_base.resolve()))
        except ValueError:
            return str(resolved)

    return {
        "schema_version": 2,
        "source_manifest": display(manifest_path),
        "source_manifest_sha256": sha256_file(manifest_path),
        "promoter_annotation": display(annotation),
        "promoter_annotation_checksum": str(
            reference["preparation"]["annotation"]["checksum"]
        ),
        "canonical_chromosomes": list(DM6_CANONICAL_CONTACT_CHROMOSOMES),
        "promoter_width_bp": 1000,
        "maximum_distance_bp": 1_000_000,
        "pseudocount_fraction": 0.01,
        "candidate_observed_over_expected_threshold": 1.0,
        "normalization": "merge_counts_then_ice_retry_v2",
        "promoter_activity": "summit_within_500bp_catalog_activity_v3",
        "link_score": "contact_weight_x_promoter_activity_percentile_v3",
        "contexts": [dict(row) for row in DM6_ATLAS_CONTACT_CONTEXTS],
    }
