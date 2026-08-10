"""Build a deterministic, human-editable per-library QC review table."""

from __future__ import annotations

import csv
import io
import math
import os
from pathlib import Path
import tempfile
from typing import Any

from .accessions import AcquisitionError
from .artifacts import read_final_bam_manifest


REVIEW_FIELDS = [
    "library_id",
    "assay",
    "context",
    "role",
    "layout",
    "final_total_alignments",
    "final_mapped_alignments",
    "properly_paired_fraction",
    "frip_numerator",
    "frip_denominator",
    "frip",
    "peak_count",
    "atac_tss_enrichment",
    "atac_fragments_counted",
    "atac_median_fragment_length_bp",
    "atac_fraction_lt150bp",
    "atac_fraction_180_250bp",
    "phantompeak_read_count",
    "phantompeak_estimated_fragment_shifts_bp",
    "phantompeak_nsc",
    "phantompeak_rsc",
    "phantompeak_quality_tag",
    "atac_tss_profile",
    "atac_fragment_histogram",
    "phantompeak_cross_correlation_plot",
    "multiqc_report",
    "qc_decision",
    "estimated_fragment_length_bp",
    "notes",
]


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise AcquisitionError(f"QC table {path} has no header")
        return [dict(row) for row in reader]


def _optional_float(value: str | None) -> float | None:
    if value is None or not value.strip() or value.strip().lower() in {"na", "nan"}:
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _format_number(value: float | None) -> str:
    return "" if value is None else f"{value:.6g}"


def parse_tss_profile(path: Path) -> float | None:
    """Return max central TSS signal divided by mean signal in terminal flanks.

    The profile spans -2 kb to +2 kb in 10-bp bins. The background is the
    combined terminal 100 bp at both ends and the numerator is the maximum bin
    in the central +/-50 bp.
    """

    profile: list[float] | None = None
    with path.open(newline="", encoding="utf-8") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) < 3:
                continue
            if fields[0].strip().lower() == "bins":
                continue
            try:
                values = [float(value) for value in fields[2:] if value.strip()]
            except ValueError:
                continue
            if values and all(math.isfinite(value) for value in values):
                profile = values
                break
    if not profile:
        return None

    bin_width = 4000.0 / len(profile)
    centers = [-2000.0 + (index + 0.5) * bin_width for index in range(len(profile))]
    background_values = [
        value
        for center, value in zip(centers, profile)
        if center < -1900.0 or center >= 1900.0
    ]
    central_values = [
        value
        for center, value in zip(centers, profile)
        if -50.0 <= center < 50.0
    ]
    if not background_values or not central_values:
        return None
    background = sum(background_values) / len(background_values)
    return max(central_values) / background if background > 0 else None


def parse_fragment_histogram(path: Path) -> dict[str, str]:
    counts: list[tuple[int, int]] = []
    for line_number, row in enumerate(_read_tsv_rows(path), start=2):
        try:
            length = int(row["fragment_length"])
            count = int(row["count"])
        except (KeyError, ValueError) as error:
            raise AcquisitionError(
                f"Invalid ATAC fragment histogram row {path}:{line_number}"
            ) from error
        if length < 1 or count < 0:
            raise AcquisitionError(
                f"Invalid ATAC fragment histogram row {path}:{line_number}"
            )
        if count:
            counts.append((length, count))

    total = sum(count for _, count in counts)
    if total == 0:
        return {
            "atac_fragments_counted": "0",
            "atac_median_fragment_length_bp": "",
            "atac_fraction_lt150bp": "",
            "atac_fraction_180_250bp": "",
        }

    midpoint = (total + 1) // 2
    cumulative = 0
    median = counts[-1][0]
    for length, count in counts:
        cumulative += count
        if cumulative >= midpoint:
            median = length
            break
    below_150 = sum(count for length, count in counts if length < 150)
    between_180_250 = sum(
        count for length, count in counts if 180 <= length <= 250
    )
    return {
        "atac_fragments_counted": str(total),
        "atac_median_fragment_length_bp": str(median),
        "atac_fraction_lt150bp": _format_number(below_150 / total),
        "atac_fraction_180_250bp": _format_number(between_180_250 / total),
    }


def parse_cross_correlation(path: Path) -> dict[str, str]:
    lines = [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(lines) != 1 or len(lines[0]) < 11:
        raise AcquisitionError(f"Unexpected phantompeakqualtools output: {path}")
    fields = lines[0]
    shifts = fields[2].strip()
    return {
        "phantompeak_read_count": fields[1].strip(),
        "phantompeak_estimated_fragment_shifts_bp": shifts,
        "phantompeak_nsc": _format_number(_optional_float(fields[-3])),
        "phantompeak_rsc": _format_number(_optional_float(fields[-2])),
        "phantompeak_quality_tag": fields[-1].strip(),
        "suggested_fragment_length_bp": shifts.split(",", 1)[0] if shifts else "",
    }


def _relative(path: str | Path | None, base: Path) -> str:
    if not path:
        return ""
    return Path(os.path.relpath(Path(path).resolve(), base.resolve())).as_posix()


def _write_if_changed(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_text(encoding="utf-8") == content:
        return
    temporary_name = None
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


def build_review_table(
    *,
    final_bam_manifest: Path,
    metrics_tsv: Path,
    library_qc: dict[str, dict[str, str | None]],
    multiqc_report: Path,
    output: Path,
) -> int:
    """Write one QC review row for every final-BAM manifest library."""

    artifacts = read_final_bam_manifest(
        final_bam_manifest, require_files=True, allow_rejected=True
    )
    metric_rows: dict[str, dict[str, str]] = {}
    for row in _read_tsv_rows(metrics_tsv):
        sample = row.get("sample", "")
        if not sample or sample in metric_rows:
            raise AcquisitionError(f"Invalid or duplicate sample in {metrics_tsv}: {sample!r}")
        metric_rows[sample] = row
    missing_metrics = sorted(set(artifacts) - set(metric_rows))
    unknown_metrics = sorted(set(metric_rows) - set(artifacts))
    if missing_metrics or unknown_metrics:
        details = []
        if missing_metrics:
            details.append("missing metrics: " + ", ".join(missing_metrics))
        if unknown_metrics:
            details.append("unknown metrics: " + ", ".join(unknown_metrics))
        raise AcquisitionError("QC metrics do not match the BAM manifest; " + "; ".join(details))

    rows: list[dict[str, Any]] = []
    for library_id, artifact in sorted(artifacts.items()):
        metrics = metric_rows[library_id]
        paths = library_qc.get(library_id, {})
        row: dict[str, Any] = {field: "" for field in REVIEW_FIELDS}
        row.update(
            {
                "library_id": library_id,
                "assay": artifact["assay"],
                "context": artifact["context"],
                "role": artifact["role"],
                "layout": artifact["layout"],
                "final_total_alignments": metrics.get("total_alignments", ""),
                "final_mapped_alignments": metrics.get("mapped_alignments", ""),
                "frip_numerator": metrics.get("frip_numerator", ""),
                "frip_denominator": metrics.get("frip_denominator", ""),
                "frip": metrics.get("frip", ""),
                "peak_count": metrics.get("peak_count", ""),
                "multiqc_report": _relative(multiqc_report, output.parent),
            }
        )
        total = _optional_float(metrics.get("total_alignments"))
        paired = _optional_float(metrics.get("properly_paired"))
        if total and paired is not None:
            row["properly_paired_fraction"] = _format_number(paired / total)

        tss_profile = paths.get("tss_profile")
        fragment_histogram = paths.get("fragment_histogram")
        cross_correlation = paths.get("cross_correlation")
        cross_correlation_plot = paths.get("cross_correlation_plot")
        if tss_profile:
            row["atac_tss_enrichment"] = _format_number(
                parse_tss_profile(Path(tss_profile))
            )
            row["atac_tss_profile"] = _relative(tss_profile, output.parent)
        if fragment_histogram:
            row.update(parse_fragment_histogram(Path(fragment_histogram)))
            row["atac_fragment_histogram"] = _relative(
                fragment_histogram, output.parent
            )
        cross_metrics: dict[str, str] = {}
        if cross_correlation:
            cross_metrics = parse_cross_correlation(Path(cross_correlation))
            row.update(
                {
                    key: value
                    for key, value in cross_metrics.items()
                    if key in REVIEW_FIELDS
                }
            )
            row["phantompeak_cross_correlation_plot"] = _relative(
                cross_correlation_plot, output.parent
            )

        status = artifact["qc_status"]
        row["qc_decision"] = {
            "pending_review": "",
            "accepted": "pass",
            "rejected": "fail",
        }[status]
        row["estimated_fragment_length_bp"] = artifact.get(
            "estimated_fragment_length_bp", ""
        )
        if (
            not row["estimated_fragment_length_bp"]
            and artifact["assay"] == "chip_histone"
            and artifact["layout"] == "single"
        ):
            row["estimated_fragment_length_bp"] = cross_metrics.get(
                "suggested_fragment_length_bp", ""
            )
        row["notes"] = artifact["notes"] if status != "pending_review" else ""
        rows.append(row)

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=REVIEW_FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    _write_if_changed(output, buffer.getvalue())
    return len(rows)
