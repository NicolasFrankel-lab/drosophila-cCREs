#!/usr/bin/env python3
"""Create a checksummed final-BAM manifest from exact library-ID filenames."""

from __future__ import annotations

import argparse
import csv
import io
import os
from pathlib import Path
import tempfile

from short_read_processing.artifacts import (
    FINAL_BAM_FILTERING_CONTRACT,
    FINAL_BAM_QC_STATUSES,
    sha256_file,
)
from short_read_processing.cli import cli_main
from short_read_processing.sample_sheet import ASSAY_ALIASES, DEFAULT_SCHEMA, read_sample_sheet


def _relative(path: Path, base: Path) -> str:
    return Path(os.path.relpath(path.resolve(), base.resolve())).as_posix()


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_sheet", type=Path)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--assay", required=True, choices=("atac", "h3k27ac", "chip_histone"))
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--genome", required=True, choices=("dm6",))
    parser.add_argument("--layout", required=True, choices=("single", "paired"))
    parser.add_argument(
        "--qc-status",
        required=True,
        choices=tuple(sorted(FINAL_BAM_QC_STATUSES)),
        help=(
            "QC decision for every selected library; rejected manifests cannot be reused"
        ),
    )
    parser.add_argument("--source-project", required=True)
    parser.add_argument("--source-run-id", required=True)
    args = parser.parse_args()

    selected_assay = ASSAY_ALIASES.get(args.assay, args.assay)
    rows = read_sample_sheet(args.sample_sheet.resolve(), schema_path=args.schema.resolve())
    by_library = {}
    for row in rows:
        if row["assay"] == selected_assay:
            by_library.setdefault(str(row["library_id"]), row)
    if not by_library:
        parser.error(f"sample sheet contains no {args.assay} libraries")

    output = args.output.resolve()
    bam_dir = args.bam_dir.resolve()
    buffer = io.StringIO()
    fieldnames = [
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
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for library_id, row in sorted(by_library.items()):
        bam = bam_dir / f"{library_id}.final.bam"
        bai = bam_dir / f"{library_id}.final.bam.bai"
        missing = [path for path in (bam, bai) if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                f"Library {library_id!r} is incomplete; missing: "
                + ", ".join(str(path) for path in missing)
            )
        writer.writerow(
            {
                "library_id": library_id,
                "assay": args.assay,
                "context": row["context"],
                "role": row["role"],
                "layout": args.layout,
                "bam": _relative(bam, output.parent),
                "bai": _relative(bai, output.parent),
                "genome": args.genome,
                "filtering_contract": FINAL_BAM_FILTERING_CONTRACT,
                "bam_sha256": sha256_file(bam),
                "bai_sha256": sha256_file(bai),
                "qc_status": args.qc_status,
                "estimated_fragment_length_bp": "",
                "source_project": args.source_project,
                "source_run_id": args.source_run_id,
                "notes": "",
            }
        )
    _write_if_changed(output, buffer.getvalue())
    print(f"Wrote {len(by_library)} final-BAM libraries to {output}")
    return 0


if __name__ == "__main__":
    cli_main(main)
