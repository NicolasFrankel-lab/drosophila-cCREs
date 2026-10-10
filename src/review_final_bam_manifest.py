#!/usr/bin/env python3
"""Apply explicit QC decisions to a complete final-BAM manifest."""

from __future__ import annotations

import argparse
import csv
import io
import os
from pathlib import Path
import tempfile

from short_read_processing.accessions import AcquisitionError
from short_read_processing.artifacts import read_final_bam_manifest
from short_read_processing.cli import cli_main
from short_read_processing.sample_sheet import ASSAY_ALIASES, read_delimited_rows


FIELDS = [
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
    "estimated_fragment_length_bp",
    "source_project",
    "source_run_id",
    "notes",
]
REVIEW_FIELDS = {
    "library_id",
    "estimated_fragment_length_bp",
    "notes",
}


def _relative(path: str, base: Path) -> str:
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


def _read_decisions(
    path: Path, *, assay: str | None = None
) -> dict[str, dict[str, str]]:
    columns, rows = read_delimited_rows(path)
    missing = sorted(REVIEW_FIELDS - set(columns))
    if missing:
        raise AcquisitionError(
            f"{path} is missing required columns: " + ", ".join(missing)
        )
    if assay is not None:
        if "assay" not in columns:
            raise AcquisitionError(f"{path} needs an assay column for --assay")
        wanted = ASSAY_ALIASES.get(assay, assay)
        rows = [
            row
            for row in rows
            if ASSAY_ALIASES.get(row["assay"], row["assay"]) == wanted
        ]
    decision_column = "qc_decision" if "qc_decision" in columns else "qc_status"
    if decision_column not in columns:
        raise AcquisitionError(
            f"{path} is missing required column: qc_decision"
        )
    decisions: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        library_id = row["library_id"]
        if not library_id:
            raise AcquisitionError(f"library_id on line {line_number} is blank")
        if library_id in decisions:
            raise AcquisitionError(f"Duplicate QC decision for {library_id!r}")
        raw_status = row[decision_column].strip().lower()
        status = {
            "pass": "accepted",
            "fail": "rejected",
            "accepted": "accepted",
            "rejected": "rejected",
        }.get(raw_status)
        if status is None:
            raise AcquisitionError(
                f"QC decision for {library_id!r} must be pass or fail"
            )
        notes = row["notes"].strip()
        if status == "rejected" and not notes:
            raise AcquisitionError(
                f"Rejected library {library_id!r} requires a reason in notes"
            )
        decisions[library_id] = {
            "qc_status": status,
            "estimated_fragment_length_bp": row[
                "estimated_fragment_length_bp"
            ].strip(),
            "notes": notes,
        }
    return decisions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_manifest", type=Path)
    review_input = parser.add_mutually_exclusive_group(required=True)
    review_input.add_argument(
        "--review-table",
        type=Path,
        help="Completed library-review.tsv with pass/fail decisions",
    )
    review_input.add_argument(
        "--decisions",
        type=Path,
        help="Legacy decision-table option",
    )
    parser.add_argument(
        "--assay",
        choices=("atac", "h3k27ac"),
        help="Use only decision rows for this assay (multi-assay decision tables)",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    source = read_final_bam_manifest(
        args.input_manifest.resolve(),
        require_files=True,
        allow_rejected=True,
    )
    review_table = args.review_table or args.decisions
    decisions = _read_decisions(review_table.resolve(), assay=args.assay)
    missing = sorted(set(source) - set(decisions))
    unknown = sorted(set(decisions) - set(source))
    if missing or unknown:
        details = []
        if missing:
            details.append("missing decisions: " + ", ".join(missing))
        if unknown:
            details.append("unknown libraries: " + ", ".join(unknown))
        raise AcquisitionError(
            "QC decisions do not match the manifest; " + "; ".join(details)
        )

    output = args.output.resolve()
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for library_id, artifact in sorted(source.items()):
        decision = decisions[library_id]
        estimated = decision["estimated_fragment_length_bp"]
        is_single_histone = (
            artifact["layout"] == "single"
            and artifact["assay"] == "chip_histone"
        )
        if is_single_histone:
            try:
                estimated_value = int(estimated)
            except ValueError as error:
                raise AcquisitionError(
                    f"Single-end histone library {library_id!r} requires a "
                    "positive estimated_fragment_length_bp"
                ) from error
            if estimated_value <= 0:
                raise AcquisitionError(
                    f"Single-end histone library {library_id!r} requires a "
                    "positive estimated_fragment_length_bp"
                )
            estimated = str(estimated_value)
        elif estimated:
            raise AcquisitionError(
                f"Library {library_id!r} must leave estimated_fragment_length_bp "
                "blank unless it is single-end histone ChIP"
            )
        writer.writerow(
            {
                **artifact,
                "bam": _relative(artifact["bam"], output.parent),
                "bai": _relative(artifact["bai"], output.parent),
                "qc_status": decision["qc_status"],
                "estimated_fragment_length_bp": estimated,
                "notes": decision["notes"],
            }
        )

    _write_if_changed(output, buffer.getvalue())
    accepted = sum(
        decision["qc_status"] == "accepted" for decision in decisions.values()
    )
    print(
        f"Wrote {len(decisions)} reviewed libraries to {output} "
        f"({accepted} accepted, {len(decisions) - accepted} rejected)"
    )
    return 0


if __name__ == "__main__":
    cli_main(main)
