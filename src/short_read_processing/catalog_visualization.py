"""Build BED and normalized mean-signal tracks for catalog visualization."""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Iterable, Iterator

from .activity import read_master_elements, sha256_file, write_json_if_changed
from .activity_tmm import (
    TMM_BACKGROUND_METHOD,
    TMM_FACTOR_FIELDS,
    _atomic_text_if_changed,
)


REGULATORY_CLASS_COLORS = {
    "promoter_associated": "215,48,39",
    "proximal_enhancer_like": "253,174,97",
    "distal_enhancer_like": "69,117,180",
    "unclassified_no_tss_on_contig": "117,112,179",
}
ATAC_BROWSER_EXTENSION_BP = 150


def _read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected chromosome and size")
            chrom, size_text = fields
            size = int(size_text)
            if not chrom or chrom in seen or size < 1:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            rows.append((chrom, size))
            seen.add(chrom)
    if not rows:
        raise ValueError(f"Chromosome-size file is empty: {path}")
    return rows


def _read_context_matrix(
    path: Path,
    contexts: list[str],
    elements,
) -> dict[tuple[str, str], bool]:
    membership: dict[tuple[str, str], bool] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "master_dhs_id",
            "chrom",
            "start",
            "end",
            "summit",
            "context_n",
            *contexts,
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: context matrix lacks required columns")
        rows = list(reader)
    if len(rows) != len(elements):
        raise ValueError("Context matrix and master registry have different row counts")
    for line_number, (row, element) in enumerate(zip(rows, elements), start=2):
        if (
            row["master_dhs_id"] != element.master_id
            or row["chrom"] != element.chrom
            or int(row["start"]) != element.start
            or int(row["end"]) != element.end
            or int(row["summit"]) != element.summit
        ):
            raise ValueError(f"{path}:{line_number}: master coordinates differ")
        values = []
        for context in contexts:
            if row[context] not in {"0", "1"}:
                raise ValueError(f"{path}:{line_number}: invalid {context} membership")
            value = row[context] == "1"
            values.append(value)
            membership[(element.master_id, context)] = value
        matrix_context_n = sum(
            int(row[field])
            for field in reader.fieldnames
            if field not in {"master_dhs_id", "chrom", "start", "end", "summit", "context_n"}
        )
        if int(row["context_n"]) != matrix_context_n:
            raise ValueError(f"{path}:{line_number}: context_n is inconsistent")
    return membership


def _element_bed_content(
    *,
    path: Path,
    context: str,
    element_by_id: dict[str, Any],
    membership: dict[tuple[str, str], bool],
) -> tuple[str, dict[str, int]]:
    lines = []
    class_counts: dict[str, int] = {}
    seen = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "master_dhs_id",
            "chrom",
            "start",
            "end",
            "summit",
            "context",
            "context_membership",
            "regulatory_class",
            "nearest_tss_distance_bp",
            "blacklist_overlap",
            "h3k27ac_max_500_background_percentile",
            "h3k27ac_active",
            "activity_state",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: context-element table lacks required columns")
        for line_number, row in enumerate(reader, start=2):
            master_id = row["master_dhs_id"]
            element = element_by_id.get(master_id)
            if element is None or master_id in seen:
                raise ValueError(f"{path}:{line_number}: unknown or duplicate master DHS")
            if (
                row["context"] != context
                or row["context_membership"] != "1"
                or not membership[(master_id, context)]
                or row["chrom"] != element.chrom
                or int(row["start"]) != element.start
                or int(row["end"]) != element.end
                or int(row["summit"]) != element.summit
            ):
                raise ValueError(f"{path}:{line_number}: inconsistent context element")
            element_class = row["regulatory_class"]
            blacklisted = row["blacklist_overlap"] == "1"
            if row["blacklist_overlap"] not in {"0", "1"}:
                raise ValueError(f"{path}:{line_number}: invalid blacklist annotation")
            color = (
                "0,0,0"
                if blacklisted
                else REGULATORY_CLASS_COLORS.get(element_class, "105,105,105")
            )
            percentile_text = row["h3k27ac_max_500_background_percentile"]
            percentile = float(percentile_text) if percentile_text else None
            if percentile is not None and (
                not math.isfinite(percentile) or not 0 <= percentile <= 1
            ):
                raise ValueError(f"{path}:{line_number}: invalid background percentile")
            if row["h3k27ac_active"] not in {"0", "1"}:
                raise ValueError(f"{path}:{line_number}: invalid H3K27ac activity call")
            active = row["h3k27ac_active"] == "1"
            distance = row["nearest_tss_distance_bp"] or "NA"
            percentile_label = "NA" if percentile is None else f"{percentile:.6g}"
            state = row["activity_state"]
            activity_label = "active" if active else "inactive"
            blacklist_label = "blacklist_overlap" if blacklisted else "blacklist_clear"
            name = (
                f"{master_id}|{element_class}|tss_distance={distance}|"
                f"max3_background_percentile={percentile_label}|{activity_label}|{state}|"
                f"{blacklist_label}"
            )
            score = 0 if percentile is None else round(1000 * percentile)
            lines.append(
                f"{element.chrom}\t{element.start}\t{element.end}\t{name}\t"
                f"{score}\t.\t{element.summit}\t{element.summit + 1}\t{color}\n"
            )
            class_counts[element_class] = class_counts.get(element_class, 0) + 1
            seen.add(master_id)
    return "".join(lines), class_counts


def build_catalog_beds(
    *,
    master_bed: Path,
    summit_bed: Path,
    context_matrix: Path,
    element_paths: dict[str, Path],
    output_master_bed: Path,
    output_context_dhs: dict[str, Path],
    output_element_beds: dict[str, Path],
    output_manifest: Path,
) -> dict[str, Any]:
    """Write portable master, context-DHS, and activity-annotated BED tracks."""

    contexts = list(element_paths)
    if (
        not contexts
        or set(output_context_dhs) != set(contexts)
        or set(output_element_beds) != set(contexts)
    ):
        raise ValueError("Catalog BED context mappings must be complete and identical")
    elements = read_master_elements(master_bed, summit_bed)
    element_by_id = {element.master_id: element for element in elements}
    membership = _read_context_matrix(context_matrix, contexts, elements)

    _atomic_text_if_changed(
        output_master_bed,
        master_bed.read_text(encoding="utf-8"),
    )
    context_metrics = {}
    for context in contexts:
        dhs_lines = []
        for element in elements:
            if membership[(element.master_id, context)]:
                dhs_lines.append(
                    f"{element.chrom}\t{element.start}\t{element.end}\t"
                    f"{element.master_id}\t0\t.\t{element.summit}\t"
                    f"{element.summit + 1}\t0,145,130\n"
                )
        _atomic_text_if_changed(output_context_dhs[context], "".join(dhs_lines))
        element_content, class_counts = _element_bed_content(
            path=element_paths[context],
            context=context,
            element_by_id=element_by_id,
            membership=membership,
        )
        _atomic_text_if_changed(output_element_beds[context], element_content)
        element_n = sum(class_counts.values())
        if element_n != len(dhs_lines):
            raise ValueError(f"Context-element count differs from context DHSs for {context}")
        context_metrics[context] = {
            "context_dhs_count": len(dhs_lines),
            "context_element_count": element_n,
            "elements_by_regulatory_class": dict(sorted(class_counts.items())),
        }

    outputs = {
        "master_dhs_bed": {
            "path": str(output_master_bed.resolve()),
            "sha256": sha256_file(output_master_bed),
        },
        "contexts": {
            context: {
                "context_dhs_bed": {
                    "path": str(output_context_dhs[context].resolve()),
                    "sha256": sha256_file(output_context_dhs[context]),
                },
                "context_elements_bed": {
                    "path": str(output_element_beds[context].resolve()),
                    "sha256": sha256_file(output_element_beds[context]),
                },
            }
            for context in contexts
        },
    }
    manifest = {
        "status": "ok",
        "schema_version": 2,
        "method": "catalog_igv_bed_tracks_v2",
        "context_dhs_definition": "master DHS coordinates with context_matrix membership=1",
        "context_element_definition": "all context-member master DHSs; background-percentile activity is annotation, not an export filter",
        "bed_format": "BED9; score is max3 background percentile times 1000 and thick interval marks the representative summit",
        "blacklist_item_rgb": "0,0,0",
        "regulatory_class_item_rgb": REGULATORY_CLASS_COLORS,
        "master_dhs_count": len(elements),
        "context_metrics": context_metrics,
        "inputs": {
            "master_bed": sha256_file(master_bed),
            "summit_bed": sha256_file(summit_bed),
            "context_matrix": sha256_file(context_matrix),
            "context_element_tables": {
                context: sha256_file(element_paths[context]) for context in contexts
            },
        },
        "outputs": outputs,
    }
    write_json_if_changed(output_manifest, manifest)
    return manifest


def read_track_factors(
    *,
    factor_path: Path,
    library_ids: list[str],
    assay: str,
    context: str,
) -> dict[str, dict[str, float | int | str]]:
    """Return validated TMM effective sizes for one context and assay."""

    selected: dict[str, dict[str, float | int | str]] = {}
    wanted = set(library_ids)
    if not library_ids or len(wanted) != len(library_ids):
        raise ValueError("Mean BigWig libraries must be non-empty and unique")
    with factor_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TMM_FACTOR_FIELDS:
            raise ValueError(f"{factor_path}: unexpected TMM-factor columns")
        for row in reader:
            library_id = row["library_id"]
            if library_id not in wanted:
                continue
            effective_size = float(row["effective_library_size"])
            total_units = int(row["total_units"])
            normalization_factor = float(row["tmm_normalization_factor"])
            if (
                library_id in selected
                or row["assay"] != assay
                or row["context"] != context
                or row["normalization_method"] != TMM_BACKGROUND_METHOD
                or not math.isfinite(effective_size)
                or effective_size <= 0
                or not math.isfinite(normalization_factor)
                or normalization_factor <= 0
                or total_units < 1
            ):
                raise ValueError(f"Invalid TMM factor for {library_id}")
            selected[library_id] = {
                "total_units": total_units,
                "tmm_normalization_factor": normalization_factor,
                "effective_library_size": effective_size,
                "scale": 1_000_000.0 / effective_size,
            }
    if set(selected) != wanted:
        missing = sorted(wanted - set(selected))
        raise ValueError("TMM factors are missing libraries: " + ", ".join(missing))
    return {library_id: selected[library_id] for library_id in library_ids}


def mean_unionbedg_rows(
    lines: Iterable[str],
    *,
    chromosome_sizes: list[tuple[str, int]],
    library_n: int,
) -> Iterator[tuple[str, int, int, float]]:
    """Validate unionbedg rows and yield nonzero arithmetic mean segments."""

    if library_n < 1:
        raise ValueError("library_n must be positive")
    chrom_order = {chrom: index for index, (chrom, _) in enumerate(chromosome_sizes)}
    sizes = dict(chromosome_sizes)
    previous: tuple[int, int, int] | None = None
    for line_number, line in enumerate(lines, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 3 + library_n:
            raise ValueError(f"unionbedg line {line_number}: unexpected column count")
        chrom = fields[0]
        if chrom not in chrom_order:
            raise ValueError(f"unionbedg line {line_number}: unknown chromosome {chrom}")
        start, end = int(fields[1]), int(fields[2])
        values = [float(value) for value in fields[3:]]
        if (
            start < 0
            or end <= start
            or end > sizes[chrom]
            or any(not math.isfinite(value) or value < 0 for value in values)
        ):
            raise ValueError(f"unionbedg line {line_number}: invalid interval or value")
        current = (chrom_order[chrom], start, end)
        if previous is not None and (
            current[0] < previous[0]
            or (current[0] == previous[0] and start < previous[2])
        ):
            raise ValueError("unionbedg rows are unsorted or overlapping")
        previous = current
        mean_value = sum(values) / library_n
        if mean_value > 0:
            yield chrom, start, end, mean_value


def centered_extension_flanks(target_width: int) -> tuple[int, int]:
    """Return left/right additions that expand a one-base event to target width."""

    if target_width < 1:
        raise ValueError("Target width must be positive")
    left = target_width // 2
    right = target_width - left - 1
    return left, right


def _coverage_bedgraph(
    *,
    unit_path: Path,
    chrom_sizes_path: Path,
    scale: float,
    output_path: Path,
    threads: int,
    centered_extension_bp: int | None = None,
) -> None:
    if centered_extension_bp is not None and centered_extension_bp < 1:
        raise ValueError("Centered extension must be positive")
    with output_path.open("wb") as output_handle:
        decompressor = subprocess.Popen(
            ["pigz", "-p", str(max(1, threads)), "-dc", str(unit_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert decompressor.stdout is not None
        expander = None
        coverage_input = decompressor.stdout
        if centered_extension_bp is not None:
            left, right = centered_extension_flanks(centered_extension_bp)
            expander = subprocess.Popen(
                [
                    "bedtools",
                    "slop",
                    "-l",
                    str(left),
                    "-r",
                    str(right),
                    "-g",
                    str(chrom_sizes_path),
                    "-i",
                    "stdin",
                ],
                stdin=coverage_input,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            coverage_input.close()
            assert expander.stdout is not None
            coverage_input = expander.stdout
        coverage = subprocess.Popen(
            [
                "bedtools",
                "genomecov",
                "-bg",
                "-scale",
                f"{scale:.17g}",
                "-i",
                "stdin",
                "-g",
                str(chrom_sizes_path),
            ],
            stdin=coverage_input,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        coverage_input.close()
        assert coverage.stdout is not None
        sorter = subprocess.Popen(
            [
                "sort",
                "-k1,1",
                "-k2,2n",
                "--temporary-directory",
                str(output_path.parent),
            ],
            stdin=coverage.stdout,
            stdout=output_handle,
            stderr=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C"},
        )
        coverage.stdout.close()
        _, coverage_stderr = coverage.communicate()
        expander_stderr = b""
        if expander is not None:
            _, expander_stderr = expander.communicate()
        _, decompressor_stderr = decompressor.communicate()
        _, sorter_stderr = sorter.communicate()
    if (
        decompressor.returncode
        or (expander is not None and expander.returncode)
        or coverage.returncode
        or sorter.returncode
    ):
        raise RuntimeError(
            f"Could not make coverage for {unit_path}: "
            + decompressor_stderr.decode(errors="replace")
            + expander_stderr.decode(errors="replace")
            + coverage_stderr.decode(errors="replace")
            + sorter_stderr.decode(errors="replace")
        )


def _replace_if_changed(temporary: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if (
        output.is_file()
        and output.stat().st_size == temporary.stat().st_size
        and sha256_file(output) == sha256_file(temporary)
    ):
        temporary.unlink()
        return
    os.replace(temporary, output)


def build_context_mean_bigwig(
    *,
    unit_paths: dict[str, Path],
    factor_path: Path,
    chrom_sizes_path: Path,
    assay: str,
    context: str,
    output_bigwig: Path,
    output_metrics: Path,
    atac_extension_bp: int = ATAC_BROWSER_EXTENSION_BP,
    threads: int = 2,
) -> dict[str, Any]:
    """Write the mean of per-library background-TMM normalized coverage tracks."""

    if assay not in {"atac", "h3k27ac"}:
        raise ValueError(f"Unsupported track assay: {assay}")
    if atac_extension_bp < 1:
        raise ValueError("ATAC browser extension must be positive")
    library_ids = list(unit_paths)
    factors = read_track_factors(
        factor_path=factor_path,
        library_ids=library_ids,
        assay=assay,
        context=context,
    )
    chromosome_sizes = sorted(_read_chrom_sizes(chrom_sizes_path))
    output_bigwig.parent.mkdir(parents=True, exist_ok=True)
    temporary_bigwig: Path | None = None
    try:
        with tempfile.TemporaryDirectory(
            prefix=f".{context}.{assay}.", dir=output_bigwig.parent
        ) as temporary_dir_text:
            temporary_dir = Path(temporary_dir_text)
            bedgraphs = []
            for index, library_id in enumerate(library_ids):
                bedgraph = temporary_dir / f"{index:04d}.bedgraph"
                _coverage_bedgraph(
                    unit_path=unit_paths[library_id],
                    chrom_sizes_path=chrom_sizes_path,
                    scale=float(factors[library_id]["scale"]),
                    output_path=bedgraph,
                    threads=threads,
                    centered_extension_bp=(
                        atac_extension_bp if assay == "atac" else None
                    ),
                )
                bedgraphs.append(bedgraph)

            try:
                import pyBigWig
            except ImportError as error:  # pragma: no cover - rule environment owns this
                raise RuntimeError("pyBigWig is required to write catalog tracks") from error

            union = None
            if len(bedgraphs) == 1:
                mean_stream = bedgraphs[0].open(encoding="utf-8")
            else:
                union = subprocess.Popen(
                    [
                        "bedtools",
                        "unionbedg",
                        "-filler",
                        "0",
                        "-i",
                        *map(str, bedgraphs),
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                assert union.stdout is not None
                mean_stream = union.stdout

            with tempfile.NamedTemporaryFile(
                prefix=f".{output_bigwig.name}.",
                dir=output_bigwig.parent,
                delete=False,
            ) as handle:
                temporary_bigwig = Path(handle.name)
            bigwig = pyBigWig.open(str(temporary_bigwig), "w")
            if bigwig is None:
                raise RuntimeError(f"Could not create BigWig: {temporary_bigwig}")
            try:
                bigwig.addHeader(chromosome_sizes)
                chrom_buffer: list[str] = []
                start_buffer: list[int] = []
                end_buffer: list[int] = []
                value_buffer: list[float] = []
                for chrom, start, end, value in mean_unionbedg_rows(
                    mean_stream,
                    chromosome_sizes=chromosome_sizes,
                    library_n=len(library_ids),
                ):
                    chrom_buffer.append(chrom)
                    start_buffer.append(start)
                    end_buffer.append(end)
                    value_buffer.append(value)
                    if len(chrom_buffer) == 100_000:
                        bigwig.addEntries(
                            chrom_buffer,
                            start_buffer,
                            ends=end_buffer,
                            values=value_buffer,
                        )
                        chrom_buffer.clear()
                        start_buffer.clear()
                        end_buffer.clear()
                        value_buffer.clear()
                if chrom_buffer:
                    bigwig.addEntries(
                        chrom_buffer,
                        start_buffer,
                        ends=end_buffer,
                        values=value_buffer,
                    )
            finally:
                bigwig.close()
                mean_stream.close()
            if union is not None:
                union_stderr = union.stderr.read() if union.stderr is not None else ""
                if union.wait() != 0:
                    raise RuntimeError("bedtools unionbedg failed: " + union_stderr)
            _replace_if_changed(temporary_bigwig, output_bigwig)
            temporary_bigwig = None
    finally:
        if temporary_bigwig is not None and temporary_bigwig.exists():
            temporary_bigwig.unlink()

    unit_semantics = (
        f"Tn5 insertion records expanded to centered "
        f"{atac_extension_bp}-bp intervals"
        if assay == "atac"
        else "ChIP fragments"
    )
    metrics = {
        "status": "ok",
        "schema_version": 2,
        "method": "mean_background_tmm_normalized_browser_coverage_v2",
        "assay": assay,
        "context": context,
        "unit_semantics": unit_semantics,
        "atac_insertion_extension_bp": (
            atac_extension_bp if assay == "atac" else None
        ),
        "value_definition": (
            "arithmetic mean across accepted libraries of basewise coverage "
            "multiplied by 1e6/effective_library_size; ATAC insertions are "
            f"expanded to centered, chromosome-clipped "
            f"{atac_extension_bp}-bp intervals before coverage"
            if assay == "atac"
            else "arithmetic mean across accepted libraries of basewise "
            "fragment coverage multiplied by 1e6/effective_library_size"
        ),
        "normalization_method": TMM_BACKGROUND_METHOD,
        "library_n": len(library_ids),
        "libraries": {
            library_id: {
                **factors[library_id],
                "unit_bed": str(unit_paths[library_id].resolve()),
                "unit_bed_sha256": sha256_file(unit_paths[library_id]),
            }
            for library_id in library_ids
        },
        "inputs": {
            "normalization_factors_sha256": sha256_file(factor_path),
            "chrom_sizes_sha256": sha256_file(chrom_sizes_path),
        },
        "output": {
            "path": str(output_bigwig.resolve()),
            "sha256": sha256_file(output_bigwig),
        },
    }
    write_json_if_changed(output_metrics, metrics)
    return metrics
