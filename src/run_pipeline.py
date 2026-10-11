#!/usr/bin/env python3
"""Run accession sample sheet -> FASTQs -> resolved YAML -> Snakemake outputs."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from short_read_processing.accessions import AcquisitionError
from short_read_processing.cli import (
    add_download_arguments,
    cli_main,
    execute_download,
)
from short_read_processing.configuration import (
    generate_activity_config,
    generate_catalog_links_config,
    generate_configs,
    generate_footprint_config,
    generate_resume_config,
)
from short_read_processing.sample_sheet import DEFAULT_SCHEMA, sample_sheet_accessions
from short_read_processing.stage_checkpoints import read_stage_checkpoint
from short_read_processing.workflow_config import (
    OUTPUT_STAGES,
    validate_stage_selection,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def rule_threads(value: str) -> str:
    """Validate a Snakemake RULE=THREADS override."""
    rule, separator, threads = value.rpartition("=")
    if not separator or not rule or not threads.isdigit() or int(threads) < 1:
        raise argparse.ArgumentTypeError("use RULE=THREADS with a positive thread count")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_sheet", type=Path, help="Canonical accession CSV/TSV")
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--project", default="short-read")
    parser.add_argument("--run-id", default="baseline")
    parser.add_argument("--genome", choices=("dm6",), default="dm6")
    parser.add_argument("--config-dir", type=Path, default=Path("configs"))
    parser.add_argument("--reference-root", type=Path, default=Path("references"))
    parser.add_argument("--snakefile", type=Path, default=REPO_ROOT / "workflow" / "Snakefile")
    parser.add_argument(
        "--workflow-profile",
        type=Path,
        default=REPO_ROOT / "profiles" / "local",
    )
    parser.add_argument(
        "--cores",
        type=int,
        help="Maximum aggregate cores (local or across submitted cluster jobs)",
    )
    parser.add_argument("--jobs", type=int, help="Maximum concurrent cluster jobs")
    parser.add_argument(
        "--max-threads",
        type=int,
        help="Maximum threads/CPUs requested by any individual rule",
    )
    parser.add_argument(
        "--set-threads",
        action="append",
        default=[],
        type=rule_threads,
        metavar="RULE=THREADS",
        help="Override one rule's thread count; repeat for additional rules",
    )
    parser.add_argument(
        "--alignment-chunk-pairs",
        type=int,
        help=(
            "Align each lane in restartable chunks of this many read pairs "
            "(or reads for single-end data)"
        ),
    )
    parser.add_argument("--skip-download", action="store_true", help="Reuse --manifest")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--config-only", action="store_true")
    parser.add_argument("--snakemake-dry-run", action="store_true")
    parser.add_argument(
        "--from-stage",
        choices=(
            "accessions",
            *(stage for stage in OUTPUT_STAGES if stage != "footprints"),
            "final-bam",
        ),
        default="accessions",
        help=(
            "Completed logical boundary to reuse. Reuse is strict and has no "
            "upstream fallback; final-bam remains a compatibility alias."
        ),
    )
    parser.add_argument(
        "--until-stage",
        choices=OUTPUT_STAGES,
        help=(
            "Stop at this boundary. Defaults to QC from accessions/alignment, "
            "master from QC/master, and report from quantification/catalog/links. "
            "footprints (from master, with --final-bam-manifest) runs TOBIAS "
            "ATAC footprinting."
        ),
    )
    parser.add_argument(
        "--checkpoint-manifest",
        type=Path,
        help=(
            "Checksummed JSON checkpoint exported at --from-stage; required for "
            "trimming, QC-to-master, quantification, catalog, links, and report resumes"
        ),
    )
    parser.add_argument(
        "--final-bam-manifest",
        type=Path,
        help="Complete final-BAM manifest for an alignment or reviewed QC boundary",
    )
    parser.add_argument(
        "--master-manifest",
        type=Path,
        help="Immutable master-DHS bundle manifest for the master boundary",
    )
    parser.add_argument(
        "--catalog-manifest",
        type=Path,
        help=(
            "Immutable catalog bundle manifest for a downstream-only "
            "catalog-to-links run"
        ),
    )
    parser.add_argument(
        "--activity-bam-manifest",
        type=Path,
        action="append",
        default=[],
        help=(
            "Accepted ATAC/H3K27ac final-BAM manifest for quantification; "
            "repeat for assay-specific manifests"
        ),
    )
    parser.add_argument(
        "--report-source-root",
        type=Path,
        action="append",
        default=[],
        help=(
            "Completed upstream result root whose QC manifests, JSONs, and plots "
            "should be included in the final report; repeat as needed"
        ),
    )
    parser.add_argument(
        "--h3k27ac-activity-threshold",
        type=float,
        default=0.60,
        help=(
            "H3K27ac activity cutoff: a context-member DHS is active when its "
            "max3 signal background percentile is strictly greater than this "
            "value (default 0.60). Changing it requires a new --run-id."
        ),
    )
    parser.add_argument("--atac-minimum-replicates", type=int, default=2)
    parser.add_argument("--atac-overlap-fraction", type=float, default=0.5)
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Additional single Snakemake argument; repeat as needed",
    )
    add_download_arguments(parser)
    args = parser.parse_args()
    if args.cores is not None and args.cores < 1:
        parser.error("--cores must be positive")
    if args.jobs is not None and args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.max_threads is not None and args.max_threads < 1:
        parser.error("--max-threads must be positive")
    if args.alignment_chunk_pairs is not None and args.alignment_chunk_pairs < 1:
        parser.error("--alignment-chunk-pairs must be positive")
    if not 0 < args.h3k27ac_activity_threshold < 1:
        parser.error("--h3k27ac-activity-threshold must be in (0, 1)")
    if args.download_only and args.config_only:
        parser.error("--download-only and --config-only are mutually exclusive")
    try:
        output_stage = validate_stage_selection(args.from_stage, args.until_stage)
    except AcquisitionError as error:
        parser.error(str(error))
    if args.alignment_chunk_pairs is not None and not (
        args.from_stage in {"accessions", "trimming"}
        and output_stage in {"alignment", "qc"}
    ):
        parser.error(
            "--alignment-chunk-pairs requires an accessions/trimming run "
            "continuing through alignment or QC"
        )
    if args.from_stage == "accessions":
        if (
            args.checkpoint_manifest
            or args.final_bam_manifest
            or args.master_manifest
            or args.catalog_manifest
            or args.activity_bam_manifest
            or args.report_source_root
        ):
            parser.error(
                "artifact manifests require the matching reuse --from-stage"
            )
    else:
        if args.catalog_manifest and args.from_stage != "catalog":
            parser.error("--catalog-manifest requires --from-stage catalog")
        incompatible = (
            args.skip_download
            or args.download_only
            or args.dry_run
            or args.manifest is not None
        )
        if incompatible:
            parser.error(
                "--manifest, --skip-download, --download-only, and acquisition --dry-run "
                "are valid only with --from-stage accessions"
            )
        if args.from_stage == "trimming" and not args.checkpoint_manifest:
            parser.error("--from-stage trimming requires --checkpoint-manifest")
        if args.from_stage == "alignment":
            if bool(args.checkpoint_manifest) == bool(args.final_bam_manifest):
                parser.error(
                    "--from-stage alignment requires exactly one of "
                    "--checkpoint-manifest or --final-bam-manifest"
                )
        if args.from_stage == "qc":
            if not args.checkpoint_manifest:
                parser.error("--from-stage qc requires --checkpoint-manifest")
            if output_stage == "master" and not args.final_bam_manifest:
                parser.error(
                    "QC-to-master reuse also requires the reviewed "
                    "--final-bam-manifest"
                )
        if args.from_stage == "final-bam":
            if not args.final_bam_manifest:
                parser.error("--from-stage final-bam requires --final-bam-manifest")
            if args.master_manifest:
                parser.error("--master-manifest requires --from-stage master")
            if args.activity_bam_manifest:
                parser.error("activity options require --from-stage quantification")
            if args.report_source_root:
                parser.error("report sources require --from-stage quantification")
        if args.from_stage == "master" and output_stage == "footprints":
            if not args.master_manifest or not args.final_bam_manifest:
                parser.error(
                    "--until-stage footprints requires --master-manifest and the "
                    "reviewed ATAC --final-bam-manifest"
                )
            if args.activity_bam_manifest or args.report_source_root:
                parser.error(
                    "activity BAM manifests and report sources do not apply to footprints"
                )
        elif args.from_stage == "master":
            if not args.master_manifest:
                parser.error("--from-stage master requires --master-manifest")
            if args.final_bam_manifest:
                parser.error(
                    "--from-stage master consumes the frozen master bundle, not BAMs"
                )
            if output_stage != "master" and not args.activity_bam_manifest:
                parser.error(
                    "continuing after master requires --activity-bam-manifest"
                )
            if output_stage == "master" and args.activity_bam_manifest:
                parser.error(
                    "--activity-bam-manifest is used only when continuing after master"
                )
        if args.from_stage == "quantification":
            legacy_inputs = bool(args.master_manifest or args.activity_bam_manifest)
            if args.checkpoint_manifest and legacy_inputs:
                parser.error(
                    "quantification resume uses either --checkpoint-manifest or the "
                    "legacy master/activity manifests, not both"
                )
            if not args.checkpoint_manifest and (
                not args.master_manifest or not args.activity_bam_manifest
            ):
                parser.error(
                    "--from-stage quantification requires --checkpoint-manifest; "
                    "legacy commands require both master and activity BAM manifests"
                )
            if args.final_bam_manifest:
                parser.error(
                    "quantification mode uses --activity-bam-manifest"
                )
        if args.from_stage == "catalog":
            if bool(args.checkpoint_manifest) == bool(args.catalog_manifest):
                parser.error(
                    "--from-stage catalog requires exactly one of "
                    "--checkpoint-manifest or --catalog-manifest"
                )
            if args.catalog_manifest:
                if output_stage != "links":
                    parser.error(
                        "--catalog-manifest currently supports only --until-stage links"
                    )
                if (
                    args.final_bam_manifest
                    or args.master_manifest
                    or args.activity_bam_manifest
                    or args.report_source_root
                ):
                    parser.error(
                        "--catalog-manifest cannot be combined with upstream artifact "
                        "manifests or report sources"
                    )
        if args.from_stage in {"links", "report"} and not args.checkpoint_manifest:
            parser.error(f"--from-stage {args.from_stage} requires --checkpoint-manifest")

    generates_activity_config = (
        args.from_stage == "master" and output_stage not in {"master", "footprints"}
    ) or (args.from_stage == "quantification" and not args.checkpoint_manifest)
    if args.h3k27ac_activity_threshold != 0.60 and not generates_activity_config:
        parser.error(
            "--h3k27ac-activity-threshold applies only when a new master or "
            "quantification run generates the activity configuration"
        )

    sample_sheet = args.sample_sheet.resolve()
    if args.from_stage == "alignment" and args.checkpoint_manifest:
        alignment_checkpoint = read_stage_checkpoint(
            args.checkpoint_manifest.resolve(), expected_stage="alignment"
        )
        try:
            args.final_bam_manifest = Path(
                alignment_checkpoint["artifacts"]["final_bam_manifest"]["path"]
            )
        except KeyError as error:
            parser.error(
                "alignment checkpoint does not contain a final_bam_manifest artifact"
            )
    manifest: Path | None = None
    if args.from_stage == "accessions":
        accessions = sample_sheet_accessions(
            sample_sheet, schema_path=args.schema.resolve()
        )
        manifest = (args.manifest or args.output_dir / "download_manifest.tsv").resolve()
        if args.skip_download:
            if not manifest.is_file():
                raise FileNotFoundError(f"Manifest does not exist: {manifest}")
            print(f"Reusing download manifest: {manifest}")
        else:
            manifest = execute_download(accessions, args)
            if args.dry_run:
                print("Download dry-run complete; processing was not started")
                return 0
        if args.download_only:
            return 0
    else:
        manifests = [
            artifact
            for artifact in (
                args.checkpoint_manifest,
                args.final_bam_manifest,
                args.master_manifest,
                args.catalog_manifest,
                *args.activity_bam_manifest,
            )
            if artifact is not None
        ]
        print(f"Strict reuse-only mode: starting from {args.from_stage}")
        for artifact in manifests:
            print(f"  manifest={artifact.resolve()}")
        print("Fallback before the selected checkpoint is disabled")

    resume_in_place = bool(args.checkpoint_manifest) and (
        args.from_stage in {
            "trimming",
            "quantification",
            "catalog",
            "links",
            "report",
        }
        or (args.from_stage == "qc" and output_stage == "qc")
    )
    if args.catalog_manifest:
        configs = [
            generate_catalog_links_config(
                sample_sheet_path=sample_sheet,
                catalog_manifest_path=args.catalog_manifest.resolve(),
                output_dir=args.config_dir.resolve(),
                project=args.project,
                run_id=args.run_id,
                reference_root=args.reference_root,
                path_base=REPO_ROOT,
                schema_path=args.schema.resolve(),
                genome=args.genome,
            )
        ]
    elif output_stage == "footprints":
        configs = [
            generate_footprint_config(
                sample_sheet_path=sample_sheet,
                final_bam_manifests=[args.final_bam_manifest.resolve()],
                master_manifest_path=args.master_manifest.resolve(),
                output_dir=args.config_dir.resolve(),
                project=args.project,
                run_id=args.run_id,
                reference_root=args.reference_root,
                path_base=REPO_ROOT,
                schema_path=args.schema.resolve(),
                genome=args.genome,
            )
        ]
    elif resume_in_place:
        configs = [
            generate_resume_config(
                checkpoint_manifest_path=args.checkpoint_manifest.resolve(),
                start_stage=args.from_stage,
                output_stage=output_stage,
                sample_sheet_path=sample_sheet,
                output_dir=args.config_dir.resolve(),
                path_base=REPO_ROOT,
                alignment_chunk_pairs=args.alignment_chunk_pairs,
            )
        ]
    elif generates_activity_config:
        configs = [
            generate_activity_config(
                sample_sheet_path=sample_sheet,
                final_bam_manifests=[
                    path.resolve() for path in args.activity_bam_manifest
                ],
                master_manifest_path=args.master_manifest.resolve(),
                output_dir=args.config_dir.resolve(),
                project=args.project,
                run_id=args.run_id,
                reference_root=args.reference_root,
                path_base=REPO_ROOT,
                require_files=True,
                schema_path=args.schema.resolve(),
                genome=args.genome,
                output_stage=output_stage,
                start_stage=args.from_stage,
                report_source_roots=[
                    path.resolve() for path in args.report_source_root
                ],
                h3k27ac_activity_threshold=args.h3k27ac_activity_threshold,
            )
        ]
    else:
        configs = generate_configs(
            manifest_path=manifest,
            sample_sheet_path=sample_sheet,
            output_dir=args.config_dir.resolve(),
            project=args.project,
            run_id=args.run_id,
            reference_root=args.reference_root,
            path_base=REPO_ROOT,
            require_fastq_files=True,
            schema_path=args.schema.resolve(),
            genome=args.genome,
            atac_minimum_replicates=args.atac_minimum_replicates,
            atac_overlap_fraction=args.atac_overlap_fraction,
            alignment_chunk_pairs=args.alignment_chunk_pairs,
            input_stage=(
                "final-bam"
                if args.from_stage in {"alignment", "qc", "final-bam"}
                else args.from_stage
            ),
            start_stage=args.from_stage,
            output_stage=output_stage,
            final_bam_manifest_path=(
                args.final_bam_manifest.resolve()
                if args.final_bam_manifest
                else None
            ),
            master_manifest_path=(
                args.master_manifest.resolve() if args.master_manifest else None
            ),
            qc_checkpoint_manifest_path=(
                args.checkpoint_manifest.resolve()
                if args.from_stage == "qc"
                else None
            ),
        )
    for config_path in configs:
        print(f"Resolved workflow config: {config_path}")
    if args.config_only:
        return 0

    snakemake = shutil.which("snakemake") or str(Path(sys.executable).with_name("snakemake"))
    if not Path(snakemake).is_file() and not shutil.which(snakemake):
        raise FileNotFoundError("snakemake is not available in PATH")
    for config_path in configs:
        command = [
            snakemake,
            "--snakefile",
            str(args.snakefile.resolve()),
            "--configfile",
            str(config_path),
            "--workflow-profile",
            str(args.workflow_profile.resolve()),
        ]
        if args.cores is not None:
            command.extend(["--cores", str(args.cores)])
        if args.jobs is not None:
            command.extend(["--jobs", str(args.jobs)])
        if args.max_threads is not None:
            command.extend(["--max-threads", str(args.max_threads)])
        if args.set_threads:
            command.extend(["--set-threads", *args.set_threads])
        if args.snakemake_dry_run:
            command.append("--dry-run")
        command.extend(args.snakemake_arg)
        print("Running: " + " ".join(command))
        subprocess.run(command, cwd=REPO_ROOT, check=True)
    return 0


if __name__ == "__main__":
    cli_main(main)
