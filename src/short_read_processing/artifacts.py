"""Read and validate explicit reusable-artifact manifests."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

from .accessions import AcquisitionError
from .sample_sheet import ASSAY_ALIASES, read_delimited_rows


SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
FINAL_BAM_FILTERING_CONTRACT = "short-read-processing-final-v2"
FINAL_BAM_QC_STATUSES = {"pending_review", "accepted", "rejected"}
FINAL_BAM_REQUIRED_COLUMNS = {
    "library_id",
    "assay",
    "context",
    "role",
    "layout",
    "bam",
    "bai",
    "genome",
    "filtering_contract",
    "bam_sha256",
    "bai_sha256",
    "qc_status",
}
MASTER_FILE_FIELDS = (
    "master_bed",
    "summits_bed",
    "membership_tsv",
    "context_matrix_tsv",
    "stats_json",
)
MASTER_REQUIRED_COLUMNS = {
    "genome",
    "method",
    "input_filtering_contract",
    "source_project",
    "source_run_id",
    *MASTER_FILE_FIELDS,
    *(f"{field}_sha256" for field in MASTER_FILE_FIELDS),
}
CATALOG_FILE_FIELDS = (
    "catalog",
    "metrics",
    "provenance",
    "resolved_config",
)
CATALOG_REQUIRED_COLUMNS = {
    "genome",
    "method",
    "contexts",
    "source_project",
    "source_run_id",
    *CATALOG_FILE_FIELDS,
    *(f"{field}_sha256" for field in CATALOG_FILE_FIELDS),
}


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _require_columns(path: Path, columns: list[str], required: set[str]) -> None:
    missing = sorted(required - set(columns))
    if missing:
        raise AcquisitionError(f"{path} is missing required columns: " + ", ".join(missing))


def _safe_id(value: str, *, field: str, line: int) -> str:
    if not SAFE_ID_RE.fullmatch(value):
        raise AcquisitionError(f"{field} on line {line} is invalid: {value!r}")
    return value


def _sha256(value: str, *, field: str, line: int) -> str:
    normalized = value.removeprefix("sha256:").lower()
    if not SHA256_RE.fullmatch(normalized):
        raise AcquisitionError(
            f"{field} on line {line} must be a 64-character SHA-256 digest"
        )
    return normalized


def _artifact_path(value: str, *, manifest: Path, field: str, line: int) -> Path:
    if not value:
        raise AcquisitionError(f"{field} on line {line} is blank")
    path = Path(value)
    return (path if path.is_absolute() else manifest.parent / path).resolve()


def read_final_bam_manifest(
    path: Path,
    *,
    require_files: bool = True,
    allow_rejected: bool = False,
) -> dict[str, dict[str, str]]:
    """Return final-BAM rows keyed by biological library ID."""

    path = path.resolve()
    columns, rows = read_delimited_rows(path)
    _require_columns(path, columns, FINAL_BAM_REQUIRED_COLUMNS)
    if not rows:
        raise AcquisitionError(f"Final-BAM manifest {path} has no data rows")

    by_library: dict[str, dict[str, str]] = {}
    artifact_paths: set[Path] = set()
    for line, raw in enumerate(rows, start=2):
        library_id = _safe_id(raw["library_id"], field="library_id", line=line)
        if library_id in by_library:
            raise AcquisitionError(
                f"Final-BAM manifest library_id {library_id!r} occurs more than once"
            )
        assay = ASSAY_ALIASES.get(raw["assay"], raw["assay"])
        if assay not in {"atac", "chip_histone"}:
            raise AcquisitionError(f"assay on line {line} is invalid: {raw['assay']!r}")
        context = _safe_id(raw["context"], field="context", line=line)
        role = raw["role"]
        if role != "treatment":
            raise AcquisitionError(f"role on line {line} is invalid: {role!r}")
        layout = raw["layout"]
        if layout not in {"single", "paired"}:
            raise AcquisitionError(f"layout on line {line} is invalid: {layout!r}")
        if raw["qc_status"] not in FINAL_BAM_QC_STATUSES:
            raise AcquisitionError(
                f"qc_status for {library_id!r} is invalid: {raw['qc_status']!r}"
            )
        notes = raw.get("notes", "").strip()
        if raw["qc_status"] == "rejected" and not allow_rejected:
            raise AcquisitionError(
                f"Final BAM {library_id!r} has qc_status='rejected'"
            )
        if raw["qc_status"] == "rejected" and not notes:
            raise AcquisitionError(
                f"Final BAM {library_id!r} has qc_status='rejected' but no notes"
            )
        estimated_fragment_length = raw.get(
            "estimated_fragment_length_bp", ""
        ).strip()
        if estimated_fragment_length:
            try:
                estimated_fragment_length_value = int(estimated_fragment_length)
            except ValueError as error:
                raise AcquisitionError(
                    f"estimated_fragment_length_bp for {library_id!r} "
                    "must be a positive integer"
                ) from error
            if estimated_fragment_length_value <= 0:
                raise AcquisitionError(
                    f"estimated_fragment_length_bp for {library_id!r} "
                    "must be a positive integer"
                )
            estimated_fragment_length = str(estimated_fragment_length_value)
        if raw["filtering_contract"] != FINAL_BAM_FILTERING_CONTRACT:
            raise AcquisitionError(
                f"Final BAM {library_id!r} has unsupported filtering_contract "
                f"{raw['filtering_contract']!r}; expected {FINAL_BAM_FILTERING_CONTRACT!r}"
            )
        bam = _artifact_path(raw["bam"], manifest=path, field="bam", line=line)
        bai = _artifact_path(raw["bai"], manifest=path, field="bai", line=line)
        if bam == bai or bam in artifact_paths or bai in artifact_paths:
            raise AcquisitionError("Final-BAM manifest contains duplicate BAM or BAI paths")
        artifact_paths.update((bam, bai))
        if require_files:
            missing = [item for item in (bam, bai) if not item.is_file()]
            if missing:
                raise AcquisitionError(
                    f"Final BAM {library_id!r} is incomplete; missing: "
                    + ", ".join(str(item) for item in missing)
                )

        by_library[library_id] = {
            "library_id": library_id,
            "assay": assay,
            "context": context,
            "role": role,
            "layout": layout,
            "bam": str(bam),
            "bai": str(bai),
            "genome": raw["genome"],
            "filtering_contract": raw["filtering_contract"],
            "bam_sha256": _sha256(raw["bam_sha256"], field="bam_sha256", line=line),
            "bai_sha256": _sha256(raw["bai_sha256"], field="bai_sha256", line=line),
            "qc_status": raw["qc_status"],
            "estimated_fragment_length_bp": estimated_fragment_length,
            "source_project": raw.get("source_project", ""),
            "source_run_id": raw.get("source_run_id", ""),
            "notes": notes,
        }
    return by_library


def read_master_manifest(
    path: Path,
    *,
    require_files: bool = True,
) -> dict[str, str]:
    """Read the one-row immutable master-DHS bundle manifest."""

    path = path.resolve()
    columns, rows = read_delimited_rows(path)
    _require_columns(path, columns, MASTER_REQUIRED_COLUMNS)
    if len(rows) != 1:
        raise AcquisitionError(f"Master manifest {path} must contain exactly one data row")
    raw = rows[0]
    result = {
        "genome": raw["genome"],
        "method": raw["method"],
        "input_filtering_contract": raw["input_filtering_contract"],
        "source_project": _safe_id(
            raw["source_project"], field="source_project", line=2
        ),
        "source_run_id": _safe_id(
            raw["source_run_id"], field="source_run_id", line=2
        ),
    }
    if result["input_filtering_contract"] != FINAL_BAM_FILTERING_CONTRACT:
        raise AcquisitionError(
            "Master manifest has unsupported input_filtering_contract "
            f"{result['input_filtering_contract']!r}; expected "
            f"{FINAL_BAM_FILTERING_CONTRACT!r}"
        )
    for field in MASTER_FILE_FIELDS:
        artifact = _artifact_path(raw[field], manifest=path, field=field, line=2)
        if require_files and not artifact.is_file():
            raise AcquisitionError(f"Master artifact does not exist: {artifact}")
        result[field] = str(artifact)
        result[f"{field}_sha256"] = _sha256(
            raw[f"{field}_sha256"],
            field=f"{field}_sha256",
            line=2,
        )
    if len({result[field] for field in MASTER_FILE_FIELDS}) != len(MASTER_FILE_FIELDS):
        raise AcquisitionError("Master manifest artifact paths must be distinct")
    return result


def read_catalog_manifest(
    path: Path,
    *,
    require_files: bool = True,
) -> dict[str, Any]:
    """Read and verify one immutable context-resolved catalog bundle."""

    path = path.resolve()
    columns, rows = read_delimited_rows(path)
    _require_columns(path, columns, CATALOG_REQUIRED_COLUMNS)
    if len(rows) != 1:
        raise AcquisitionError(f"Catalog manifest {path} must contain exactly one data row")
    raw = rows[0]
    contexts = [item.strip() for item in raw["contexts"].split(",") if item.strip()]
    if (
        not contexts
        or len(contexts) != len(set(contexts))
        or any(not SAFE_ID_RE.fullmatch(item) for item in contexts)
    ):
        raise AcquisitionError("Catalog manifest contexts must be safe and unique")
    result: dict[str, Any] = {
        "genome": _safe_id(raw["genome"], field="genome", line=2),
        "method": _safe_id(raw["method"], field="method", line=2),
        "contexts": contexts,
        "source_project": _safe_id(
            raw["source_project"], field="source_project", line=2
        ),
        "source_run_id": _safe_id(
            raw["source_run_id"], field="source_run_id", line=2
        ),
    }
    for field in CATALOG_FILE_FIELDS:
        artifact = _artifact_path(raw[field], manifest=path, field=field, line=2)
        digest = _sha256(
            raw[f"{field}_sha256"], field=f"{field}_sha256", line=2
        )
        if require_files:
            if not artifact.is_file():
                raise AcquisitionError(f"Catalog artifact does not exist: {artifact}")
            if sha256_file(artifact) != digest:
                raise AcquisitionError(f"Catalog artifact SHA-256 mismatch: {artifact}")
        result[field] = str(artifact)
        result[f"{field}_sha256"] = digest
    if len({result[field] for field in CATALOG_FILE_FIELDS}) != len(
        CATALOG_FILE_FIELDS
    ):
        raise AcquisitionError("Catalog manifest artifact paths must be distinct")

    if require_files:
        try:
            metrics = json.loads(Path(result["metrics"]).read_text(encoding="utf-8"))
            provenance = json.loads(
                Path(result["provenance"]).read_text(encoding="utf-8")
            )
            resolved = json.loads(
                Path(result["resolved_config"]).read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as error:
            raise AcquisitionError(f"Cannot read catalog metadata: {error}") from error
        if (
            metrics.get("catalog_sha256") != result["catalog_sha256"]
            or metrics.get("method") != result["method"]
            or int(metrics.get("context_count", -1)) != len(contexts)
        ):
            raise AcquisitionError("Catalog metrics disagree with the manifest")
        if (
            provenance.get("outputs", {}).get("catalog", {}).get("sha256")
            != result["catalog_sha256"]
        ):
            raise AcquisitionError("Catalog provenance disagrees with the manifest")
        if (
            resolved.get("project") != result["source_project"]
            or resolved.get("run_id") != result["source_run_id"]
            or resolved.get("reference", {}).get("name") != result["genome"]
            or resolved.get("activity", {}).get("contexts") != contexts
        ):
            raise AcquisitionError("Catalog resolved configuration disagrees with the manifest")
        semantic = resolved.get("provenance", {}).get("semantic_sha256")
        if not isinstance(semantic, str) or not SHA256_RE.fullmatch(semantic):
            raise AcquisitionError("Catalog source semantic SHA-256 is invalid")
        result["source_semantic_sha256"] = semantic
        sample_sheet_sha256 = resolved.get("provenance", {}).get(
            "sample_sheet_sha256"
        )
        if not isinstance(sample_sheet_sha256, str) or not SHA256_RE.fullmatch(
            sample_sheet_sha256
        ):
            raise AcquisitionError("Catalog source sample-sheet SHA-256 is invalid")
        result["source_sample_sheet_sha256"] = sample_sheet_sha256
    return result
