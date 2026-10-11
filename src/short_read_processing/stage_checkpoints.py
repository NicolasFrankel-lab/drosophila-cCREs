"""Read and validate deterministic logical-stage checkpoint manifests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .accessions import AcquisitionError
from .artifacts import SAFE_ID_RE, SHA256_RE, sha256_file


STAGE_CHECKPOINT_SCHEMA_VERSION = 1
LOGICAL_STAGES = (
    "trimming",
    "alignment",
    "qc",
    "master",
    "quantification",
    "catalog",
    "links",
    "report",
    "footprints",
)


def _artifact_path(value: str, *, manifest: Path, label: str) -> Path:
    if not value:
        raise AcquisitionError(f"Checkpoint artifact {label!r} has a blank path")
    path = Path(value)
    return (path if path.is_absolute() else manifest.parent / path).resolve()


def read_stage_checkpoint(
    path: Path,
    *,
    expected_stage: str | None = None,
    require_files: bool = True,
    verify_hashes: bool = True,
) -> dict[str, Any]:
    """Return one validated logical-stage checkpoint with absolute paths."""

    path = path.resolve()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcquisitionError(f"Cannot read stage checkpoint {path}: {error}") from error
    if not isinstance(raw, dict):
        raise AcquisitionError(f"Stage checkpoint {path} must contain a JSON object")
    if raw.get("schema_version") != STAGE_CHECKPOINT_SCHEMA_VERSION:
        raise AcquisitionError(f"Unsupported stage checkpoint schema in {path}")
    stage = str(raw.get("stage", ""))
    if stage not in LOGICAL_STAGES:
        raise AcquisitionError(f"Stage checkpoint {path} has invalid stage {stage!r}")
    if expected_stage is not None and stage != expected_stage:
        raise AcquisitionError(
            f"Stage checkpoint {path} records {stage!r}, expected {expected_stage!r}"
        )
    for field in ("source_project", "source_run_id"):
        if not SAFE_ID_RE.fullmatch(str(raw.get(field, ""))):
            raise AcquisitionError(f"Stage checkpoint {path} has invalid {field}")
    semantic = str(raw.get("semantic_sha256", ""))
    if not SHA256_RE.fullmatch(semantic):
        raise AcquisitionError(f"Stage checkpoint {path} has invalid semantic_sha256")
    parameters = raw.get("parameters", {})
    if not isinstance(parameters, dict):
        raise AcquisitionError(f"Stage checkpoint {path} parameters must be a mapping")
    artifacts = raw.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise AcquisitionError(f"Stage checkpoint {path} has no artifacts")

    resolved: dict[str, dict[str, str]] = {}
    seen_paths: set[Path] = set()
    for label, record in sorted(artifacts.items()):
        if not isinstance(label, str) or not label or not isinstance(record, dict):
            raise AcquisitionError(f"Stage checkpoint {path} has an invalid artifact record")
        artifact = _artifact_path(str(record.get("path", "")), manifest=path, label=label)
        digest = str(record.get("sha256", "")).lower()
        if not SHA256_RE.fullmatch(digest):
            raise AcquisitionError(
                f"Checkpoint artifact {label!r} has an invalid SHA-256 digest"
            )
        if artifact in seen_paths:
            raise AcquisitionError(f"Stage checkpoint {path} repeats artifact path {artifact}")
        seen_paths.add(artifact)
        if require_files and not artifact.is_file():
            raise AcquisitionError(f"Checkpoint artifact does not exist: {artifact}")
        if require_files and verify_hashes and sha256_file(artifact) != digest:
            raise AcquisitionError(f"Checkpoint artifact SHA-256 mismatch: {artifact}")
        resolved[label] = {"path": str(artifact), "sha256": digest}

    return {
        "schema_version": STAGE_CHECKPOINT_SCHEMA_VERSION,
        "stage": stage,
        "source_project": str(raw["source_project"]),
        "source_run_id": str(raw["source_run_id"]),
        "semantic_sha256": semantic,
        "parameters": parameters,
        "artifacts": resolved,
        "manifest": str(path),
    }
