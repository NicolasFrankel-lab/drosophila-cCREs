"""Build a context-resolved regulatory-element catalog from master DHSs."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
import csv
import gzip
import math
from pathlib import Path
import statistics
from typing import Any

from .activity import (
    _tsv_content,
    read_master_elements,
    sha256_file,
    write_deterministic_gzip,
    write_json_if_changed,
)
from .activity_tmm import (
    TMM_ACTIVITY_FIELDS,
    TMM_BACKGROUND_METHOD,
    TMM_FACTOR_FIELDS,
    _atomic_text_if_changed,
    _write_gzip_dict_rows,
)


MAX_WINDOW_ORDER = ("center_500", "left_500", "right_500")
WINDOW_ORDER = MAX_WINDOW_ORDER
CATALOG_METHOD = "background_tmm_summit_max3_nonmember_percentile_v6"
WINDOW_OFFSETS = {
    "left_500": (-750, -250),
    "center_500": (-250, 250),
    "right_500": (250, 750),
}
DEFAULT_H3K27AC_ACTIVITY_PERCENTILE_THRESHOLD = 0.60
WINDOW_FIELDS = [
    "window_id",
    "master_dhs_id",
    "chrom",
    "summit",
    "window",
    "start",
    "end",
    "width_bp",
]
WINDOW_COUNT_FIELDS = [
    *WINDOW_FIELDS,
    "library_id",
    "context",
    "raw_count",
    "total_units",
]
WINDOW_SIGNAL_VALUE_FIELDS = [
    "width_bp",
    "raw_count_sum",
    "raw_count_mean",
    "cpm",
    "cpm_per_kb",
    "normalized_cpm",
    "normalized_cpm_per_kb",
]
H3K27AC_QUANTIFICATION_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "blacklist_overlap",
    "blacklist_overlap_bp",
    "blacklist_overlap_fraction",
    "context",
    "context_membership",
    "normalization_method",
    "h3k27ac_library_n",
    "h3k27ac_library_ids",
    *(
        f"h3k27ac_{window}_{field}"
        for window in WINDOW_ORDER
        for field in WINDOW_SIGNAL_VALUE_FIELDS
    ),
    "h3k27ac_max_500_window",
    "h3k27ac_max_500_raw_count_sum",
    "h3k27ac_max_500_raw_count_mean",
    "h3k27ac_max_500_cpm",
    "h3k27ac_max_500_cpm_per_kb",
    "h3k27ac_max_500_normalized_cpm",
    "h3k27ac_max_500_normalized_cpm_per_kb",
    "h3k27ac_max_500_background_n",
    "h3k27ac_max_500_background_percentile",
    "h3k27ac_max_500_background_upper_tail_probability",
]
CATALOG_FIELDS = [
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "width_bp",
    "blacklist_overlap",
    "blacklist_overlap_bp",
    "blacklist_overlap_fraction",
    "context",
    "context_membership",
    "nearest_tss_distance_bp",
    "nearest_tss_ids",
    "regulatory_class",
    "atac_normalized_cpm_per_kb",
    "h3k27ac_library_n",
    "h3k27ac_library_ids",
    "h3k27ac_left_500_width_bp",
    "h3k27ac_left_500_raw_count_sum",
    "h3k27ac_left_500_cpm_per_kb",
    "h3k27ac_left_500_normalized_cpm_per_kb",
    "h3k27ac_center_500_width_bp",
    "h3k27ac_center_500_raw_count_sum",
    "h3k27ac_center_500_cpm_per_kb",
    "h3k27ac_center_500_normalized_cpm_per_kb",
    "h3k27ac_right_500_width_bp",
    "h3k27ac_right_500_raw_count_sum",
    "h3k27ac_right_500_cpm_per_kb",
    "h3k27ac_right_500_normalized_cpm_per_kb",
    "h3k27ac_max_500_raw_count_sum",
    "h3k27ac_max_500_raw_count_mean",
    "h3k27ac_max_500_cpm_per_kb",
    "h3k27ac_max_500_normalized_cpm_per_kb",
    "h3k27ac_max_500_window",
    "h3k27ac_max_500_background_n",
    "h3k27ac_max_500_background_percentile",
    "h3k27ac_max_500_background_upper_tail_probability",
    "h3k27ac_active",
    "activity_state",
    "combined_activity_max_500",
]
SUMMARY_FIELDS = [
    "context",
    "regulatory_class",
    "blacklist_overlap",
    "h3k27ac_active",
    "activity_state",
    "element_n",
]
WIDE_CONTEXT_FIELDS = [
    "context_membership",
    "atac_normalized_cpm_per_kb",
    "h3k27ac_max_500_raw_count_sum",
    "h3k27ac_max_500_normalized_cpm_per_kb",
    "h3k27ac_max_500_window",
    "h3k27ac_max_500_background_percentile",
    "h3k27ac_max_500_background_upper_tail_probability",
    "h3k27ac_active",
    "activity_state",
    "combined_activity_max_500",
]


def _read_chrom_sizes(path: Path) -> tuple[list[str], dict[str, int]]:
    order: list[str] = []
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"{path}:{line_number}: invalid chromosome sizes")
            chrom, size_text = fields
            size = int(size_text)
            if size < 1 or chrom in sizes:
                raise ValueError(f"{path}:{line_number}: invalid chromosome sizes")
            order.append(chrom)
            sizes[chrom] = size
    if not order:
        raise ValueError(f"{path}: chromosome sizes are empty")
    return order, sizes


def window_rows(
    *,
    master_bed: Path,
    summit_bed: Path,
    chrom_sizes: Path,
) -> list[dict[str, Any]]:
    """Return clipped left/center/right 500-bp windows around each summit."""

    order, sizes = _read_chrom_sizes(chrom_sizes)
    rank = {chrom: index for index, chrom in enumerate(order)}
    elements = read_master_elements(master_bed, summit_bed)
    if any(element.chrom not in sizes for element in elements):
        raise ValueError("Master DHS contains a chromosome absent from chromosome sizes")
    rows = []
    for element in elements:
        for window in WINDOW_ORDER:
            offset_start, offset_end = WINDOW_OFFSETS[window]
            start = min(max(element.summit + offset_start, 0), sizes[element.chrom])
            end = min(max(element.summit + offset_end, 0), sizes[element.chrom])
            rows.append(
                {
                    "window_id": f"{element.master_id}|{window}",
                    "master_dhs_id": element.master_id,
                    "chrom": element.chrom,
                    "summit": element.summit,
                    "window": window,
                    "start": start,
                    "end": end,
                    "width_bp": end - start,
                    "_rank": rank[element.chrom],
                }
            )
    rows.sort(
        key=lambda row: (
            row["_rank"],
            row["start"],
            row["end"],
            row["master_dhs_id"],
            WINDOW_ORDER.index(row["window"]),
        )
    )
    return rows


def write_window_definitions(
    *,
    master_bed: Path,
    summit_bed: Path,
    chrom_sizes: Path,
    output_table: Path,
    output_bed: Path,
) -> dict[str, int]:
    """Write all window metadata and the positive-width BED used for counting."""

    rows = window_rows(
        master_bed=master_bed,
        summit_bed=summit_bed,
        chrom_sizes=chrom_sizes,
    )
    table_rows = [{field: row[field] for field in WINDOW_FIELDS} for row in rows]
    write_deterministic_gzip(output_table, _tsv_content(WINDOW_FIELDS, table_rows))
    bed_lines = [
        f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['window_id']}\t0\t.\n"
        for row in rows
        if row["width_bp"] > 0
    ]
    _atomic_text_if_changed(output_bed, "".join(bed_lines))
    return {
        "window_count": len(rows),
        "countable_window_count": sum(row["width_bp"] > 0 for row in rows),
        "zero_width_window_count": sum(row["width_bp"] == 0 for row in rows),
    }


def _read_window_table(path: Path) -> list[dict[str, Any]]:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != WINDOW_FIELDS:
            raise ValueError(f"{path}: unexpected window columns")
        for row in reader:
            parsed = {
                **row,
                "summit": int(row["summit"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "width_bp": int(row["width_bp"]),
            }
            if (
                parsed["window"] not in WINDOW_ORDER
                or parsed["window_id"]
                != f"{parsed['master_dhs_id']}|{parsed['window']}"
                or parsed["end"] - parsed["start"] != parsed["width_bp"]
                or parsed["width_bp"] < 0
            ):
                raise ValueError(f"{path}: invalid window row")
            rows.append(parsed)
    if not rows:
        raise ValueError(f"{path}: window table is empty")
    return rows


def write_window_counts(
    *,
    window_table: Path,
    coverage_path: Path,
    total_units_path: Path,
    library_id: str,
    context: str,
    output_path: Path,
) -> dict[str, int]:
    """Join bedtools coverage counts to all windows, including clipped empties."""

    rows = _read_window_table(window_table)
    observed: dict[str, int] = {}
    with coverage_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 7:
                raise ValueError(f"{coverage_path}:{line_number}: expected BED6 + count")
            window_id = fields[3]
            count = int(fields[6])
            if window_id in observed or count < 0:
                raise ValueError(f"{coverage_path}:{line_number}: invalid coverage row")
            observed[window_id] = count
    expected = {row["window_id"] for row in rows if row["width_bp"] > 0}
    if set(observed) != expected:
        raise ValueError("Coverage output and positive-width windows differ")
    total_units = int(total_units_path.read_text(encoding="utf-8").strip())
    if total_units < 1:
        raise ValueError("H3K27ac total unit count must be positive")
    output_rows = [
        {
            **{field: row[field] for field in WINDOW_FIELDS},
            "library_id": library_id,
            "context": context,
            "raw_count": observed.get(row["window_id"], 0),
            "total_units": total_units,
        }
        for row in rows
    ]
    write_deterministic_gzip(
        output_path,
        _tsv_content(WINDOW_COUNT_FIELDS, output_rows),
    )
    return {
        "window_count": len(output_rows),
        "overlap_count_sum": sum(row["raw_count"] for row in output_rows),
        "total_units": total_units,
    }


def _read_factor_table(path: Path) -> dict[str, dict[str, Any]]:
    factors: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TMM_FACTOR_FIELDS:
            raise ValueError(f"{path}: unexpected TMM-factor columns")
        for row in reader:
            library_id = row["library_id"]
            total_units = int(row["total_units"])
            effective_size = float(row["effective_library_size"])
            if (
                library_id in factors
                or row["normalization_method"] != TMM_BACKGROUND_METHOD
                or total_units < 1
                or effective_size <= 0
                or not math.isfinite(effective_size)
            ):
                raise ValueError(f"{path}: invalid background-TMM factor")
            factors[library_id] = {
                **row,
                "total_units": total_units,
                "effective_library_size": effective_size,
            }
    return factors


def _read_window_counts(path: Path) -> list[dict[str, Any]]:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != WINDOW_COUNT_FIELDS:
            raise ValueError(f"{path}: unexpected window-count columns")
        for row in reader:
            rows.append(
                {
                    **row,
                    "summit": int(row["summit"]),
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "width_bp": int(row["width_bp"]),
                    "raw_count": int(row["raw_count"]),
                    "total_units": int(row["total_units"]),
                }
            )
    if not rows:
        raise ValueError(f"{path}: window-count table is empty")
    return rows


def _read_activity_table(
    path: Path,
) -> tuple[list[tuple[str, str]], dict[tuple[str, str], dict[str, Any]]]:
    keys = []
    rows = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TMM_ACTIVITY_FIELDS:
            raise ValueError(f"{path}: unexpected background-TMM activity columns")
        for row in reader:
            if row["normalization_method"] != TMM_BACKGROUND_METHOD:
                raise ValueError(f"{path}: expected background-TMM activity")
            key = (row["master_dhs_id"], row["context"])
            if key in rows:
                raise ValueError(f"{path}: duplicate master/context row")
            keys.append(key)
            rows[key] = row
    if not rows:
        raise ValueError(f"{path}: activity table is empty")
    return keys, rows


def _read_context_matrix(
    path: Path,
    *,
    contexts: list[str],
) -> dict[tuple[str, str], int]:
    membership = {}
    expected_prefix = ["master_dhs_id", "chrom", "start", "end", "summit", "context_n"]
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != [*expected_prefix, *contexts]:
            raise ValueError(f"{path}: context matrix columns or order differ")
        for row in reader:
            context_n = int(row["context_n"])
            values = []
            for context in contexts:
                value = int(row[context])
                if value not in {0, 1}:
                    raise ValueError(f"{path}: context membership must be binary")
                membership[(row["master_dhs_id"], context)] = value
                values.append(value)
            if sum(values) != context_n:
                raise ValueError(f"{path}: context_n differs from memberships")
    return membership


def _read_tss(path: Path) -> dict[str, tuple[list[int], dict[int, list[str]]]]:
    by_chrom: dict[str, dict[int, set[str]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected BED4+")
            chrom, start_text, end_text, identifier = fields[:4]
            start, end = int(start_text), int(end_text)
            if start < 0 or end != start + 1 or not identifier:
                raise ValueError(f"{path}:{line_number}: invalid one-base TSS")
            by_chrom.setdefault(chrom, {}).setdefault(start, set()).add(identifier)
    return {
        chrom: (
            sorted(values),
            {coordinate: sorted(identifiers) for coordinate, identifiers in values.items()},
        )
        for chrom, values in by_chrom.items()
    }


def _read_blacklist(path: Path) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track ", "browser ")):
                continue
            fields = stripped.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected BED3+")
            chrom, start_text, end_text = fields[:3]
            start, end = int(start_text), int(end_text)
            if not chrom or start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid blacklist interval")
            by_chrom.setdefault(chrom, []).append((start, end))

    merged_by_chrom = {}
    for chrom, intervals in by_chrom.items():
        merged: list[tuple[int, int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        merged_by_chrom[chrom] = ([start for start, _end in merged], merged)
    return merged_by_chrom


def blacklist_overlap(
    chrom: str,
    start: int,
    end: int,
    blacklist: dict[str, tuple[list[int], list[tuple[int, int]]]],
) -> tuple[int, float]:
    """Return unioned blacklist overlap in bases and as an interval fraction."""

    if start < 0 or end <= start:
        raise ValueError("Element coordinates must define a positive-width interval")
    if chrom not in blacklist:
        return 0, 0.0
    starts, intervals = blacklist[chrom]
    stop = bisect_left(starts, end)
    begin = max(0, bisect_left(starts, start) - 1)
    overlap_bp = sum(
        max(0, min(end, interval_end) - max(start, interval_start))
        for interval_start, interval_end in intervals[begin:stop]
    )
    return overlap_bp, overlap_bp / (end - start)


def nearest_tss(
    chrom: str,
    summit: int,
    tss: dict[str, tuple[list[int], dict[int, list[str]]]],
) -> tuple[int | None, str]:
    if chrom not in tss:
        return None, ""
    coordinates, identifiers = tss[chrom]
    index = bisect_left(coordinates, summit)
    candidates = []
    if index < len(coordinates):
        candidates.append(coordinates[index])
    if index:
        candidates.append(coordinates[index - 1])
    distance = min(abs(coordinate - summit) for coordinate in candidates)
    tied = sorted(coordinate for coordinate in candidates if abs(coordinate - summit) == distance)
    names = sorted({name for coordinate in tied for name in identifiers[coordinate]})
    return distance, ";".join(names)


def regulatory_class(distance: int | None) -> str:
    if distance is None:
        return "unclassified_no_tss_on_contig"
    if distance <= 500:
        return "promoter_associated"
    if distance <= 1000:
        return "proximal_enhancer_like"
    return "distal_enhancer_like"


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _empirical_background_scores(
    signal: float,
    ordered_background: list[float],
) -> tuple[float | None, float | None]:
    """Return a midrank percentile and conservative upper-tail probability."""

    if not ordered_background:
        return None, None
    lower = bisect_left(ordered_background, signal)
    upper = bisect_right(ordered_background, signal)
    size = len(ordered_background)
    percentile = (lower + 0.5 * (upper - lower)) / size
    upper_tail = (size - lower + 1) / (size + 1)
    return percentile, upper_tail


def validate_activity_threshold(threshold: float) -> float:
    """Return a background-percentile threshold strictly inside (0, 1)."""

    value = float(threshold)
    if not math.isfinite(value) or not 0 < value < 1:
        raise ValueError("H3K27ac activity percentile threshold must be in (0, 1)")
    return value


def activity_state_label(threshold: float) -> str:
    """Name the active state with the configured threshold, e.g. ``gt_0.60``."""

    value = validate_activity_threshold(threshold)
    text = f"{value:.2f}" if math.isclose(value * 100, round(value * 100)) else f"{value:g}"
    return f"active_h3k27ac_background_gt_{text}"


def _is_h3k27ac_active(
    *,
    member: bool,
    background_percentile: float | None,
    threshold: float = DEFAULT_H3K27AC_ACTIVITY_PERCENTILE_THRESHOLD,
) -> bool:
    return bool(
        member
        and background_percentile is not None
        and background_percentile > threshold
    )


def build_h3k27ac_quantification(
    *,
    master_bed: Path,
    summit_bed: Path,
    context_matrix: Path,
    blacklist_bed: Path,
    window_table: Path,
    window_count_paths: dict[str, Path],
    factor_table: Path,
    contexts: list[str],
    output_signal: Path,
    output_metrics: Path,
    output_provenance: Path,
) -> dict[str, Any]:
    """Write raw and background-TMM summit-window signals for all DHS contexts."""

    if not contexts or len(contexts) != len(set(contexts)):
        raise ValueError("Contexts must be non-empty and unique")
    elements = read_master_elements(master_bed, summit_bed)
    definitions = _read_window_table(window_table)
    definition_keys = [(row["window_id"], row["width_bp"]) for row in definitions]
    factors = _read_factor_table(factor_table)
    h3_libraries = sorted(
        library_id
        for library_id, factor in factors.items()
        if factor["assay"] == "h3k27ac" and factor["context"] in contexts
    )
    if set(window_count_paths) != set(h3_libraries):
        raise ValueError("Window count paths and atlas H3K27ac factors differ")

    values_by_library: dict[
        str, dict[tuple[str, str], dict[str, float | int | None]]
    ] = {}
    libraries_by_context: dict[str, list[str]] = {context: [] for context in contexts}
    for library_id in h3_libraries:
        factor = factors[library_id]
        rows = _read_window_counts(window_count_paths[library_id])
        if [(row["window_id"], row["width_bp"]) for row in rows] != definition_keys:
            raise ValueError(f"Window definitions differ for library {library_id}")
        if any(
            row["library_id"] != library_id
            or row["context"] != factor["context"]
            or row["total_units"] != factor["total_units"]
            for row in rows
        ):
            raise ValueError(f"Window-count metadata differ for library {library_id}")
        library_values = {}
        for row in rows:
            width = row["width_bp"]
            raw = row["raw_count"]
            cpm = raw * 1_000_000.0 / factor["total_units"] if width else None
            normalized_cpm = (
                raw * 1_000_000.0 / factor["effective_library_size"]
                if width
                else None
            )
            library_values[(row["master_dhs_id"], row["window"])] = {
                "raw_count": raw,
                "cpm": cpm,
                "cpm_per_kb": cpm * 1000.0 / width if cpm is not None else None,
                "normalized_cpm": normalized_cpm,
                "normalized_cpm_per_kb": (
                    normalized_cpm * 1000.0 / width
                    if normalized_cpm is not None
                    else None
                ),
                "width_bp": width,
            }
        values_by_library[library_id] = library_values
        libraries_by_context[factor["context"]].append(library_id)
    if any(not libraries_by_context[context] for context in contexts):
        raise ValueError("Every context requires at least one H3K27ac library")

    expected_keys = [
        (element.master_id, context) for context in contexts for element in elements
    ]
    membership = _read_context_matrix(context_matrix, contexts=contexts)
    if set(membership) != set(expected_keys):
        raise ValueError("Context membership rows differ from the master registry")
    blacklist = _read_blacklist(blacklist_bed)
    blacklist_by_element = {}
    for element in elements:
        overlap_bp, overlap_fraction = blacklist_overlap(
            element.chrom,
            element.start,
            element.end,
            blacklist,
        )
        blacklist_by_element[element.master_id] = {
            "blacklist_overlap": int(overlap_bp > 0),
            "blacklist_overlap_bp": overlap_bp,
            "blacklist_overlap_fraction": overlap_fraction,
        }

    quantified: dict[tuple[str, str], dict[str, Any]] = {}
    for context in contexts:
        library_ids = sorted(libraries_by_context[context])
        for element in elements:
            windows = {}
            for window in WINDOW_ORDER:
                replicate_values = [
                    values_by_library[library_id][(element.master_id, window)]
                    for library_id in library_ids
                ]
                widths = {int(value["width_bp"]) for value in replicate_values}
                if len(widths) != 1:
                    raise ValueError("Window widths differ between libraries")
                width = widths.pop()
                raw_values = [int(value["raw_count"]) for value in replicate_values]
                window_values: dict[str, Any] = {
                    "width_bp": width,
                    "raw_count_sum": sum(raw_values),
                    "raw_count_mean": statistics.fmean(raw_values),
                }
                for field in (
                    "cpm",
                    "cpm_per_kb",
                    "normalized_cpm",
                    "normalized_cpm_per_kb",
                ):
                    numeric = [
                        float(value[field])
                        for value in replicate_values
                        if value[field] is not None
                    ]
                    window_values[field] = (
                        statistics.fmean(numeric) if numeric else None
                    )
                windows[window] = window_values

            eligible = [
                window
                for window in MAX_WINDOW_ORDER
                if windows[window]["normalized_cpm_per_kb"] is not None
            ]
            if not eligible:
                raise ValueError("Every master DHS requires a countable max3 window")
            winner = max(
                eligible,
                key=lambda window: (
                    windows[window]["normalized_cpm_per_kb"],
                    -MAX_WINDOW_ORDER.index(window),
                ),
            )
            quantified[(element.master_id, context)] = {
                "windows": windows,
                "winner": winner,
                "library_ids": library_ids,
            }

    background_by_context: dict[str, list[float]] = {}
    for context in contexts:
        background_by_context[context] = sorted(
            float(
                quantified[(element.master_id, context)]["windows"][
                    quantified[(element.master_id, context)]["winner"]
                ]["normalized_cpm_per_kb"]
            )
            for element in elements
            if not membership[(element.master_id, context)]
            and not blacklist_by_element[element.master_id]["blacklist_overlap"]
        )

    rows = []
    for context in contexts:
        background = background_by_context[context]
        for element in elements:
            key = (element.master_id, context)
            aggregate = quantified[key]
            windows = aggregate["windows"]
            winner = aggregate["winner"]
            max3_signal = float(windows[winner]["normalized_cpm_per_kb"])
            percentile, upper_tail = _empirical_background_scores(
                max3_signal,
                background,
            )
            row = {
                "master_dhs_id": element.master_id,
                "chrom": element.chrom,
                "start": element.start,
                "end": element.end,
                "summit": element.summit,
                "width_bp": element.width,
                **blacklist_by_element[element.master_id],
                "context": context,
                "context_membership": membership[key],
                "normalization_method": TMM_BACKGROUND_METHOD,
                "h3k27ac_library_n": len(aggregate["library_ids"]),
                "h3k27ac_library_ids": ";".join(aggregate["library_ids"]),
                "h3k27ac_max_500_window": winner,
                "h3k27ac_max_500_background_n": len(background),
                "h3k27ac_max_500_background_percentile": percentile,
                "h3k27ac_max_500_background_upper_tail_probability": upper_tail,
            }
            for window in WINDOW_ORDER:
                for field in WINDOW_SIGNAL_VALUE_FIELDS:
                    row[f"h3k27ac_{window}_{field}"] = windows[window][field]
            for field in WINDOW_SIGNAL_VALUE_FIELDS[1:]:
                row[f"h3k27ac_max_500_{field}"] = windows[winner][field]
            rows.append(row)

    _write_gzip_dict_rows(output_signal, H3K27AC_QUANTIFICATION_FIELDS, rows)
    background_metrics = {
        context: {
            "eligible_nonmember_n": len(values),
            "zero_fraction": (
                sum(value == 0 for value in values) / len(values) if values else None
            ),
            "quantiles": {
                f"q{int(probability * 100):02d}": (
                    _quantile(values, probability) if values else None
                )
                for probability in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
            },
        }
        for context, values in background_by_context.items()
    }
    metrics = {
        "method": "background_tmm_h3k27ac_summit_max3_background_v3",
        "master_dhs_count": len(elements),
        "context_count": len(contexts),
        "row_count": len(rows),
        "window_definition": {
            "coordinate_system": "zero_based_half_open",
            "left_500": "[summit-750,summit-250)",
            "center_500": "[summit-250,summit+250)",
            "right_500": "[summit+250,summit+750)",
            "chromosome_boundary_policy": "clip; zero-width windows retain raw zero and blank normalized values",
            "replicate_aggregation": "equal_weight_mean_of_per-library_values",
            "maximum_tie_order": list(MAX_WINDOW_ORDER),
        },
        "max3_definition": (
            "maximum equal-weight context-mean background-TMM normalized CPM/kb "
            "across the left, center, and right 500-bp summit windows"
        ),
        "background_definition": {
            "population": "same-context master DHSs with context_membership=0",
            "blacklist_policy": (
                "exclude master DHSs overlapping the reference blacklist from "
                "the background; percentiles are reported for every row"
            ),
            "signal": "h3k27ac_max_500_normalized_cpm_per_kb",
            "percentile": "midrank empirical CDF",
            "upper_tail_probability": "(1 + count(background >= signal)) / (1 + background_n)",
            "by_context": background_metrics,
        },
        "output_signal_sha256": sha256_file(output_signal),
    }
    provenance = {
        "method": metrics["method"],
        "inputs": {
            "master_bed": {"path": str(master_bed.resolve()), "sha256": sha256_file(master_bed)},
            "summit_bed": {"path": str(summit_bed.resolve()), "sha256": sha256_file(summit_bed)},
            "context_matrix": {"path": str(context_matrix.resolve()), "sha256": sha256_file(context_matrix)},
            "blacklist_bed": {"path": str(blacklist_bed.resolve()), "sha256": sha256_file(blacklist_bed)},
            "window_table": {"path": str(window_table.resolve()), "sha256": sha256_file(window_table)},
            "factor_table": {"path": str(factor_table.resolve()), "sha256": sha256_file(factor_table)},
            "window_counts": {
                library_id: {"path": str(path.resolve()), "sha256": sha256_file(path)}
                for library_id, path in sorted(window_count_paths.items())
            },
        },
        "parameters": {
            key: metrics[key]
            for key in (
                "window_definition",
                "max3_definition",
                "background_definition",
            )
        },
        "outputs": {
            "signal": {
                "path": str(output_signal.resolve()),
                "sha256": sha256_file(output_signal),
            }
        },
    }
    write_json_if_changed(output_metrics, metrics)
    write_json_if_changed(output_provenance, provenance)
    return metrics


def _read_h3k27ac_quantification(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    rows: dict[tuple[str, str], dict[str, Any]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != H3K27AC_QUANTIFICATION_FIELDS:
            raise ValueError(f"{path}: unexpected H3K27ac quantification columns")
        for source in reader:
            key = (source["master_dhs_id"], source["context"])
            if key in rows or source["normalization_method"] != TMM_BACKGROUND_METHOD:
                raise ValueError(f"{path}: invalid H3K27ac quantification row")
            row: dict[str, Any] = dict(source)
            for field in (
                "start",
                "end",
                "summit",
                "width_bp",
                "blacklist_overlap",
                "blacklist_overlap_bp",
                "context_membership",
                "h3k27ac_library_n",
                "h3k27ac_max_500_background_n",
            ):
                row[field] = int(source[field])
            row["blacklist_overlap_fraction"] = float(
                source["blacklist_overlap_fraction"]
            )
            for window in WINDOW_ORDER:
                for field in WINDOW_SIGNAL_VALUE_FIELDS:
                    name = f"h3k27ac_{window}_{field}"
                    if field in {"width_bp", "raw_count_sum"}:
                        row[name] = int(source[name])
                    else:
                        row[name] = float(source[name]) if source[name] else None
            for field in (
                "h3k27ac_max_500_raw_count_sum",
            ):
                row[field] = int(source[field])
            for field in (
                "h3k27ac_max_500_raw_count_mean",
                "h3k27ac_max_500_cpm",
                "h3k27ac_max_500_cpm_per_kb",
                "h3k27ac_max_500_normalized_cpm",
                "h3k27ac_max_500_normalized_cpm_per_kb",
                "h3k27ac_max_500_background_percentile",
                "h3k27ac_max_500_background_upper_tail_probability",
            ):
                row[field] = float(source[field]) if source[field] else None
            if source["h3k27ac_max_500_window"] not in MAX_WINDOW_ORDER:
                raise ValueError(f"{path}: invalid maximum H3K27ac window")
            rows[key] = row
    if not rows:
        raise ValueError(f"{path}: H3K27ac quantification table is empty")
    return rows


def build_regulatory_catalog(
    *,
    master_bed: Path,
    summit_bed: Path,
    context_matrix: Path,
    tss_bed: Path,
    blacklist_bed: Path,
    h3k27ac_quantification: Path,
    activity_table: Path,
    contexts: list[str],
    output_catalog: Path,
    output_wide: Path,
    output_element_paths: dict[str, Path],
    output_summary: Path,
    output_metrics: Path,
    output_provenance: Path,
    activity_threshold: float = DEFAULT_H3K27AC_ACTIVITY_PERCENTILE_THRESHOLD,
) -> dict[str, Any]:
    """Call max3 background-percentile activity and write the catalog.

    Every context-member DHS is tested against the blacklist-clear,
    context-nonmember background. Blacklist-overlapping members receive the
    same call and are flagged by ``blacklist_overlap``.
    """

    activity_threshold = validate_activity_threshold(activity_threshold)
    active_state = activity_state_label(activity_threshold)

    if not contexts or len(contexts) != len(set(contexts)):
        raise ValueError("Contexts must be non-empty and unique")
    if set(output_element_paths) != set(contexts):
        raise ValueError("Context-element output paths must match contexts")
    elements = read_master_elements(master_bed, summit_bed)
    activity_keys, activity_rows = _read_activity_table(activity_table)
    expected_keys = [
        (element.master_id, context) for context in contexts for element in elements
    ]
    if set(activity_keys) != set(expected_keys):
        raise ValueError("Background-TMM activity master/context records differ")
    membership = _read_context_matrix(context_matrix, contexts=contexts)
    if set(membership) != set(expected_keys):
        raise ValueError("Context membership rows differ from activity rows")
    tss = _read_tss(tss_bed)
    blacklist = _read_blacklist(blacklist_bed)
    blacklist_by_element = {}
    for element in elements:
        overlap_bp, overlap_fraction = blacklist_overlap(
            element.chrom,
            element.start,
            element.end,
            blacklist,
        )
        blacklist_by_element[element.master_id] = {
            "blacklist_overlap": int(overlap_bp > 0),
            "blacklist_overlap_bp": overlap_bp,
            "blacklist_overlap_fraction": overlap_fraction,
        }

    quantified = _read_h3k27ac_quantification(h3k27ac_quantification)
    if set(quantified) != set(expected_keys):
        raise ValueError("H3K27ac quantification master/context records differ")
    aggregated: dict[tuple[str, str], dict[str, Any]] = {}
    for context in contexts:
        for element in elements:
            key = (element.master_id, context)
            source = quantified[key]
            annotation = blacklist_by_element[element.master_id]
            if (
                source["chrom"] != element.chrom
                or source["start"] != element.start
                or source["end"] != element.end
                or source["summit"] != element.summit
                or source["context_membership"] != membership[key]
                or any(source[field] != annotation[field] for field in annotation)
            ):
                raise ValueError(
                    "H3K27ac quantification annotations differ from catalog inputs"
                )
            windows = {
                window: {
                    field: source[f"h3k27ac_{window}_{field}"]
                    for field in WINDOW_SIGNAL_VALUE_FIELDS
                }
                for window in WINDOW_ORDER
            }
            aggregated[key] = {
                "windows": windows,
                "max_value": source[
                    "h3k27ac_max_500_normalized_cpm_per_kb"
                ],
                "winner": source["h3k27ac_max_500_window"],
                "library_n": source["h3k27ac_library_n"],
                "library_ids": source["h3k27ac_library_ids"],
                "background_n": source[
                    "h3k27ac_max_500_background_n"
                ],
                "background_percentile": source[
                    "h3k27ac_max_500_background_percentile"
                ],
                "background_upper_tail_probability": source[
                    "h3k27ac_max_500_background_upper_tail_probability"
                ],
            }

    summary_counts: dict[tuple[str, str, int, int, str], int] = {}

    def catalog_row(context: str, element) -> dict[str, Any]:
        key = (element.master_id, context)
        aggregate = aggregated[key]
        activity = activity_rows[key]
        if (
            activity["chrom"] != element.chrom
            or int(activity["start"]) != element.start
            or int(activity["end"]) != element.end
            or int(activity["summit"]) != element.summit
        ):
            raise ValueError("Activity coordinates differ from master registry")
        member = membership[key]
        distance, tss_ids = nearest_tss(element.chrom, element.summit, tss)
        element_class = regulatory_class(distance)
        active = _is_h3k27ac_active(
            member=bool(member),
            background_percentile=aggregate["background_percentile"],
            threshold=activity_threshold,
        )
        if not member:
            state = "inactive_in_context"
        elif active:
            state = active_state
        elif aggregate["max_value"] <= 0:
            state = "accessible_no_h3k27ac_signal"
        else:
            state = "below_h3k27ac_background_threshold"
        windows = aggregate["windows"]
        atac = float(activity["atac_normalized_cpm_per_kb"])
        combined = math.sqrt(atac * aggregate["max_value"])
        row = {
            "master_dhs_id": element.master_id,
            "chrom": element.chrom,
            "start": element.start,
            "end": element.end,
            "summit": element.summit,
            "width_bp": element.width,
            **blacklist_by_element[element.master_id],
            "context": context,
            "context_membership": member,
            "nearest_tss_distance_bp": distance,
            "nearest_tss_ids": tss_ids,
            "regulatory_class": element_class,
            "atac_normalized_cpm_per_kb": atac,
            "h3k27ac_library_n": aggregate["library_n"],
            "h3k27ac_library_ids": aggregate["library_ids"],
            "h3k27ac_max_500_raw_count_sum": windows[aggregate["winner"]][
                "raw_count_sum"
            ],
            "h3k27ac_max_500_raw_count_mean": windows[aggregate["winner"]][
                "raw_count_mean"
            ],
            "h3k27ac_max_500_cpm_per_kb": windows[aggregate["winner"]][
                "cpm_per_kb"
            ],
            "h3k27ac_max_500_normalized_cpm_per_kb": aggregate["max_value"],
            "h3k27ac_max_500_window": aggregate["winner"],
            "h3k27ac_max_500_background_n": aggregate["background_n"],
            "h3k27ac_max_500_background_percentile": aggregate[
                "background_percentile"
            ],
            "h3k27ac_max_500_background_upper_tail_probability": aggregate[
                "background_upper_tail_probability"
            ],
            "h3k27ac_active": int(active),
            "activity_state": state,
            "combined_activity_max_500": combined,
        }
        for window in ("left_500", "center_500", "right_500"):
            prefix = f"h3k27ac_{window}"
            row[f"{prefix}_width_bp"] = windows[window]["width_bp"]
            row[f"{prefix}_raw_count_sum"] = windows[window]["raw_count_sum"]
            row[f"{prefix}_cpm_per_kb"] = windows[window]["cpm_per_kb"]
            row[f"{prefix}_normalized_cpm_per_kb"] = windows[window][
                "normalized_cpm_per_kb"
            ]
        return row

    def catalog_rows():
        for context in contexts:
            for element in elements:
                row = catalog_row(context, element)
                summary_key = (
                    context,
                    row["regulatory_class"],
                    row["blacklist_overlap"],
                    row["h3k27ac_active"],
                    row["activity_state"],
                )
                summary_counts[summary_key] = summary_counts.get(summary_key, 0) + 1
                yield row

    _write_gzip_dict_rows(output_catalog, CATALOG_FIELDS, catalog_rows())

    wide_base_fields = [
        "master_dhs_id",
        "chrom",
        "start",
        "end",
        "summit",
        "width_bp",
        "blacklist_overlap",
        "blacklist_overlap_bp",
        "blacklist_overlap_fraction",
        "nearest_tss_distance_bp",
        "nearest_tss_ids",
        "regulatory_class",
    ]
    wide_fields = [
        *wide_base_fields,
        *(
            f"{context}__{field}"
            for context in contexts
            for field in WIDE_CONTEXT_FIELDS
        ),
    ]

    def wide_rows():
        for element in elements:
            first = catalog_row(contexts[0], element)
            row = {field: first[field] for field in wide_base_fields}
            for context in contexts:
                context_row = first if context == contexts[0] else catalog_row(context, element)
                for field in WIDE_CONTEXT_FIELDS:
                    row[f"{context}__{field}"] = context_row[field]
            yield row

    _write_gzip_dict_rows(output_wide, wide_fields, wide_rows())
    for context in contexts:
        _write_gzip_dict_rows(
            output_element_paths[context],
            CATALOG_FIELDS,
            (
                row
                for element in elements
                if (row := catalog_row(context, element))["context_membership"]
            ),
        )

    summary_rows = [
        {
            "context": context,
            "regulatory_class": element_class,
            "blacklist_overlap": blacklisted,
            "h3k27ac_active": active,
            "activity_state": state,
            "element_n": count,
        }
        for (
            context,
            element_class,
            blacklisted,
            active,
            state,
        ), count in sorted(summary_counts.items())
    ]
    summary_content = _tsv_content(SUMMARY_FIELDS, summary_rows)
    _atomic_text_if_changed(output_summary, summary_content)

    metrics = {
        "method": CATALOG_METHOD,
        "master_dhs_count": len(elements),
        "context_count": len(contexts),
        "catalog_row_count": len(elements) * len(contexts),
        "window_definition": {
            "coordinate_system": "zero_based_half_open",
            "left_500": "[summit-750,summit-250)",
            "center_500": "[summit-250,summit+250)",
            "right_500": "[summit+250,summit+750)",
            "chromosome_boundary_policy": "clip; zero-width windows retained with blank normalized values",
            "replicate_aggregation": "read from quantification-stage equal-weight context means",
            "maximum_tie_order": list(MAX_WINDOW_ORDER),
        },
        "activity_definition": {
            "signal": "h3k27ac_max_500_normalized_cpm_per_kb",
            "background_population": (
                "same-context context-nonmember master DHSs excluding "
                "blacklist-overlapping elements"
            ),
            "percentile": "midrank empirical CDF",
            "threshold": activity_threshold,
            "comparison": "strictly_greater_than",
            "active_state": active_state,
            "call_population": (
                "all context-member master DHSs; blacklist-overlapping members "
                "are called and flagged with blacklist_overlap=1"
            ),
        },
        "regulatory_class_definition": {
            "coordinate": "absolute master_summit_to_nearest_tss_distance_bp",
            "promoter_associated": "distance <= 500",
            "proximal_enhancer_like": "500 < distance <= 1000",
            "distal_enhancer_like": "distance > 1000",
            "unclassified_no_tss_on_contig": "no TSS on contig",
        },
        "blacklist_annotation_definition": {
            "coordinate": "master_dhs_interval",
            "blacklist_overlap": "1 when at least one master-DHS base overlaps the reference blacklist",
            "blacklist_overlap_bp": "unioned number of overlapping bases",
            "blacklist_overlap_fraction": "blacklist_overlap_bp divided by master-DHS width",
            "filtering_policy": (
                "retain, test and annotate all rows; blacklist-overlapping DHSs "
                "are excluded only from the activity background"
            ),
        },
        "catalog_sha256": sha256_file(output_catalog),
        "wide_catalog_sha256": sha256_file(output_wide),
        "blacklist_overlapping_master_dhs_count": sum(
            annotation["blacklist_overlap"]
            for annotation in blacklist_by_element.values()
        ),
        "context_element_count_by_context": {
            context: sum(
                membership[(element.master_id, context)] for element in elements
            )
            for context in contexts
        },
        "active_element_count_by_context": {
            context: sum(
                count
                for (
                    summary_context,
                    _element_class,
                    _blacklisted,
                    active,
                    _state,
                ), count in summary_counts.items()
                if summary_context == context and active
            )
            for context in contexts
        },
        "active_blacklist_overlap_element_count_by_context": {
            context: sum(
                count
                for (
                    summary_context,
                    _element_class,
                    blacklisted,
                    active,
                    _state,
                ), count in summary_counts.items()
                if summary_context == context and active and blacklisted
            )
            for context in contexts
        },
    }
    provenance = {
        "method": metrics["method"],
        "inputs": {
            "master_bed": {"path": str(master_bed.resolve()), "sha256": sha256_file(master_bed)},
            "summit_bed": {"path": str(summit_bed.resolve()), "sha256": sha256_file(summit_bed)},
            "context_matrix": {"path": str(context_matrix.resolve()), "sha256": sha256_file(context_matrix)},
            "tss_bed": {"path": str(tss_bed.resolve()), "sha256": sha256_file(tss_bed)},
            "blacklist_bed": {"path": str(blacklist_bed.resolve()), "sha256": sha256_file(blacklist_bed)},
            "h3k27ac_quantification": {
                "path": str(h3k27ac_quantification.resolve()),
                "sha256": sha256_file(h3k27ac_quantification),
            },
            "activity_table": {"path": str(activity_table.resolve()), "sha256": sha256_file(activity_table)},
        },
        "parameters": {
            key: metrics[key]
            for key in (
                "window_definition",
                "activity_definition",
                "regulatory_class_definition",
                "blacklist_annotation_definition",
            )
        },
        "outputs": {
            "catalog": {"path": str(output_catalog.resolve()), "sha256": sha256_file(output_catalog)},
            "wide_catalog": {"path": str(output_wide.resolve()), "sha256": sha256_file(output_wide)},
            "context_elements": {
                context: {
                    "path": str(path.resolve()),
                    "sha256": sha256_file(path),
                }
                for context, path in sorted(output_element_paths.items())
            },
            "summary": {"path": str(output_summary.resolve()), "sha256": sha256_file(output_summary)},
        },
    }
    write_json_if_changed(output_metrics, metrics)
    write_json_if_changed(output_provenance, provenance)
    return metrics
