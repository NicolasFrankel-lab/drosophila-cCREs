#!/usr/bin/env python3
"""Generate processing YAML from a validated accession sample sheet and download manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.cli import cli_main
from short_read_processing.configuration import generate_configs
from short_read_processing.sample_sheet import DEFAULT_SCHEMA


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_sheet", type=Path, help="Canonical CSV or TSV sample sheet")
    parser.add_argument("--manifest", type=Path, required=True, help="Completed download manifest")
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--output-dir", type=Path, default=Path("configs"))
    parser.add_argument("--project", default="short-read")
    parser.add_argument("--run-id", default="baseline")
    parser.add_argument("--genome", choices=("dm6",), default="dm6")
    parser.add_argument("--reference-root", type=Path, default=Path("references"))
    parser.add_argument("--atac-minimum-replicates", type=int, default=2)
    parser.add_argument("--atac-overlap-fraction", type=float, default=0.5)
    parser.add_argument(
        "--until-stage",
        choices=("trimming", "alignment", "qc"),
        default="qc",
        help=(
            "Logical workflow stopping point; master construction requires a "
            "separate run from a manually reviewed final-BAM manifest"
        ),
    )
    parser.add_argument(
        "--path-base",
        type=Path,
        default=Path.cwd(),
        help="Base for portable FASTQ/reference paths (default: current directory)",
    )
    parser.add_argument(
        "--allow-missing-fastqs",
        action="store_true",
        help="Allow configuration from a dry-run manifest",
    )
    args = parser.parse_args()

    outputs = generate_configs(
        manifest_path=args.manifest,
        sample_sheet_path=args.sample_sheet,
        output_dir=args.output_dir,
        project=args.project,
        run_id=args.run_id,
        reference_root=args.reference_root,
        path_base=args.path_base,
        require_fastq_files=not args.allow_missing_fastqs,
        schema_path=args.schema,
        genome=args.genome,
        atac_minimum_replicates=args.atac_minimum_replicates,
        atac_overlap_fraction=args.atac_overlap_fraction,
        output_stage=args.until_stage,
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    cli_main(main)
