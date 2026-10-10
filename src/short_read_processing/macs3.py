"""Render validated MACS3 callpeak arguments and declared bedGraph outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .accessions import AcquisitionError


def callpeak_bedgraph_paths(output_dir: Path, name: str) -> tuple[Path, Path]:
    return (
        output_dir / f"{name}_treat_pileup.bdg",
        output_dir / f"{name}_control_lambda.bdg",
    )


def callpeak_arguments(
    peak_config: dict[str, Any],
    *,
    treatment_bam: Path,
    name: str,
    genome_size: str | int,
    output_dir: Path,
) -> list[str]:
    """Return a shell-safe argv list for IP-only q-value MACS3 callpeak."""

    if peak_config.get("command") != "callpeak":
        raise AcquisitionError("MACS3 callpeak arguments require command=callpeak")
    if peak_config.get("qvalue") is None:
        raise AcquisitionError("MACS3 callpeak requires a qvalue")
    if not peak_config.get("write_bedgraph") or not peak_config.get("spmr"):
        raise AcquisitionError("MACS3 callpeak requires write_bedgraph=true and spmr=true")

    arguments = [
        "macs3",
        "callpeak",
        "-t",
        str(treatment_bam),
        "-f",
        str(peak_config["format"]),
        "-g",
        str(genome_size),
        "-n",
        name,
        "--outdir",
        str(output_dir),
        "-q",
        str(peak_config["qvalue"]),
        "-B",
        "--SPMR",
    ]
    if peak_config.get("broad"):
        arguments.extend(["--broad", "--broad-cutoff", str(peak_config["broad_cutoff"])])
    if peak_config.get("nomodel"):
        arguments.append("--nomodel")
    if peak_config.get("shift") is not None:
        arguments.extend(["--shift", str(peak_config["shift"])])
    if peak_config.get("extsize") is not None:
        arguments.extend(["--extsize", str(peak_config["extsize"])])
    return arguments


def atac_qpois_callpeak_arguments(
    peak_config: dict[str, Any],
    *,
    insertion_bed: Path,
    name: str,
    genome_size: str | int,
    output_dir: Path,
) -> list[str]:
    """Return MACS3 arguments for the two-ended Tn5 qpois ATAC branch."""

    if peak_config.get("command") != "callpeak":
        raise AcquisitionError("ATAC qpois arguments require command=callpeak")
    if peak_config.get("format") != "BED":
        raise AcquisitionError("ATAC qpois callpeak requires two-ended insertion BED input")
    if peak_config.get("broad") or not peak_config.get("nomodel"):
        raise AcquisitionError("ATAC qpois callpeak requires narrow, nomodel peaks")
    if peak_config.get("spmr"):
        raise AcquisitionError("ATAC qpois must use unscaled pileup and local lambda")
    for field in ("qvalue", "shift", "extsize"):
        if peak_config.get(field) is None:
            raise AcquisitionError(f"ATAC qpois callpeak requires {field}")

    return [
        "macs3",
        "callpeak",
        "-t",
        str(insertion_bed),
        "-f",
        "BED",
        "-g",
        str(genome_size),
        "-n",
        name,
        "--outdir",
        str(output_dir),
        "-q",
        str(peak_config["qvalue"]),
        "--nomodel",
        "--shift",
        str(peak_config["shift"]),
        "--extsize",
        str(peak_config["extsize"]),
        "--keep-dup",
        "all",
        "-B",
    ]
