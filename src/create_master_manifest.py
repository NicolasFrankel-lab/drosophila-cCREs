#!/usr/bin/env python3
"""Create a checksummed manifest for one complete master-DHS bundle."""

from __future__ import annotations

import argparse
import csv
import io
import os
from pathlib import Path
import tempfile

from short_read_processing.artifacts import (
    FINAL_BAM_FILTERING_CONTRACT,
    MASTER_FILE_FIELDS,
    sha256_file,
)
from short_read_processing.cli import cli_main


FILENAMES = {
    "master_bed": "master_dhs.bed",
    "summits_bed": "master_dhs_summits.bed",
    "membership_tsv": "master_dhs_membership.tsv",
    "context_matrix_tsv": "master_dhs_context_matrix.tsv",
    "stats_json": "master_dhs.json",
}


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
    parser.add_argument("master_dir", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--genome", required=True, choices=("dm6",))
    parser.add_argument("--method", required=True)
    parser.add_argument("--source-project", required=True)
    parser.add_argument("--source-run-id", required=True)
    args = parser.parse_args()

    master_dir = args.master_dir.resolve()
    output = args.output.resolve()
    artifacts = {field: master_dir / FILENAMES[field] for field in MASTER_FILE_FIELDS}
    missing = [path for path in artifacts.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Master bundle is incomplete; missing: " + ", ".join(str(path) for path in missing)
        )

    fieldnames = [
        "genome",
        "method",
        "input_filtering_contract",
        "source_project",
        "source_run_id",
    ]
    fieldnames.extend(
        item for field in MASTER_FILE_FIELDS for item in (field, f"{field}_sha256")
    )
    row = {
        "genome": args.genome,
        "method": args.method,
        "input_filtering_contract": FINAL_BAM_FILTERING_CONTRACT,
        "source_project": args.source_project,
        "source_run_id": args.source_run_id,
    }
    for field, path in artifacts.items():
        row[field] = Path(os.path.relpath(path, output.parent)).as_posix()
        row[f"{field}_sha256"] = sha256_file(path)
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerow(row)
    _write_if_changed(output, buffer.getvalue())
    print(f"Wrote master-DHS bundle manifest to {output}")
    return 0


if __name__ == "__main__":
    cli_main(main)
