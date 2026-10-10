"""Annotate GTF promoters and assign enhancers to their nearest annotated TSS."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from contextlib import contextmanager
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Iterable

from .accessions import AcquisitionError
from .artifacts import sha256_file


ATTRIBUTE_RE = re.compile(
    r'(?:^|;\s*)(gene_id|gene_name|transcript_id)\s+"([^"]+)"'
)
PROMOTER_SUPPORT_DISTANCE_BP = 500


PROMOTER_FIELDS = (
    "promoter_id",
    "gene_id",
    "gene_name",
    "chrom",
    "tss",
    "strand",
    "start",
    "end",
    "transcript_ids",
    "transcript_n",
    "selection_rule",
)



@dataclass(frozen=True)
class Promoter:
    promoter_id: str
    gene_id: str
    gene_name: str
    chrom: str
    tss: int
    strand: str
    start: int
    end: int
    transcript_ids: str


@dataclass(frozen=True)
class Element:
    master_dhs_id: str
    chrom: str
    start: int
    end: int
    summit: int
    regulatory_class: str
    activity_state: str
    atac_signal: float
    h3k27ac_background_percentile: float | None
    h3k27ac_active: int
    combined_activity: float
    blacklist_overlap: int


@dataclass(frozen=True)
class PromoterActivity:
    dhs_ids: str
    accessible: int
    atac_signal: float
    h3k27ac_background_percentile: float | None
    combined_activity: float
    activity_score: float
    active: int


def _open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8", newline="")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8", newline="")
    )


def _attributes(value: str) -> dict[str, str]:
    return {key: item for key, item in ATTRIBUTE_RE.findall(value)}


def _chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or int(fields[1]) < 1:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            sizes[canonical_chromosome(fields[0])] = int(fields[1])
    return sizes


def _centered_window(center: int, width: int, chromosome_size: int) -> tuple[int, int]:
    start = max(0, center - width // 2)
    end = min(chromosome_size, center + (width - width // 2))
    return start, end


def _atomic_tsv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with _output_text(temporary, compressed=path.suffix == ".gz") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fields),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _output_text(path: Path, *, compressed: bool):
    """Open a deterministic UTF-8 output stream for one temporary path."""

    if not compressed:
        with path.open("w", encoding="utf-8", newline="") as handle:
            yield handle
        return
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as compressed_handle:
            with io.TextIOWrapper(
                compressed_handle, encoding="utf-8", newline=""
            ) as text_handle:
                yield text_handle



CHECKSUM_RE = re.compile(r"^(?:md5:[0-9a-f]{32}|sha256:[0-9a-f]{64})$")


def canonical_chromosome(value: str) -> str:
    return value if value.startswith("chr") else f"chr{value}"


def verify_reported_checksum(path: Path, checksum: str) -> None:
    """Validate one source-reported MD5/SHA-256 when a digest is available."""

    if not checksum:
        return
    if not CHECKSUM_RE.fullmatch(checksum):
        raise AcquisitionError(f"Invalid reported checksum for {path}")
    algorithm, expected = checksum.split(":", 1)
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != expected:
        raise AcquisitionError(f"Reported {algorithm} checksum mismatch: {path}")


def write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".json", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def build_promoter_table(
    *,
    annotation: Path,
    chrom_sizes: Path,
    canonical_chromosomes: list[str] | None,
    promoter_width: int,
    promoter_id_prefix: str = "DM6PROM",
    annotation_checksum: str | None = None,
    output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Create one promoter node per distinct gene/TSS from the current GTF."""

    if promoter_width < 1:
        raise ValueError("Promoter width must be positive")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", promoter_id_prefix):
        raise ValueError("Promoter ID prefix is invalid")
    if annotation_checksum:
        verify_reported_checksum(annotation, annotation_checksum)
    sizes = _chrom_sizes(chrom_sizes)
    selected_chromosomes = (
        [canonical_chromosome(value) for value in canonical_chromosomes]
        if canonical_chromosomes is not None
        else list(sizes)
    )
    canonical = set(selected_chromosomes)
    missing = canonical - set(sizes)
    if missing:
        raise ValueError(f"Chromosome sizes lack canonical sequences: {sorted(missing)}")

    grouped: dict[tuple[str, str, int, str], dict[str, Any]] = {}
    with _open_text(annotation) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            chromosome = canonical_chromosome(fields[0])
            if chromosome not in canonical:
                continue
            attributes = _attributes(fields[8])
            gene_id = attributes.get("gene_id", "")
            transcript_id = attributes.get("transcript_id", "")
            strand = fields[6]
            if not gene_id or not transcript_id or strand not in {"+", "-"}:
                continue
            tss = int(fields[3]) - 1 if strand == "+" else int(fields[4]) - 1
            key = (gene_id, chromosome, tss, strand)
            record = grouped.setdefault(
                key,
                {
                    "gene_name": attributes.get("gene_name", gene_id),
                    "transcripts": set(),
                },
            )
            record["transcripts"].add(transcript_id)
    if not grouped:
        raise ValueError(f"No canonical transcript promoters found in {annotation}")

    chromosome_order = {
        canonical_chromosome(chromosome): index
        for index, chromosome in enumerate(selected_chromosomes)
    }
    ordered = sorted(
        grouped,
        key=lambda key: (chromosome_order[key[1]], key[2], key[0], key[3]),
    )
    rows = []
    for index, (gene_id, chromosome, tss, strand) in enumerate(ordered, start=1):
        start, end = _centered_window(tss, promoter_width, sizes[chromosome])
        record = grouped[(gene_id, chromosome, tss, strand)]
        transcripts = sorted(record["transcripts"])
        rows.append(
            {
                "promoter_id": f"{promoter_id_prefix}{index:08d}",
                "gene_id": gene_id,
                "gene_name": record["gene_name"],
                "chrom": chromosome,
                "tss": tss,
                "strand": strand,
                "start": start,
                "end": end,
                "transcript_ids": ";".join(transcripts),
                "transcript_n": len(transcripts),
                "selection_rule": "all_distinct_gene_tss",
            }
        )
    _atomic_tsv(output, PROMOTER_FIELDS, rows)
    metrics = {
        "schema_version": 1,
        "annotation": str(annotation),
        "annotation_sha256": sha256_file(annotation),
        "annotation_reported_checksum": annotation_checksum,
        "chrom_sizes": str(chrom_sizes),
        "chrom_sizes_sha256": sha256_file(chrom_sizes),
        "canonical_chromosomes": canonical_chromosomes,
        "chromosomes": selected_chromosomes,
        "promoter_width_bp": promoter_width,
        "promoter_id_prefix": promoter_id_prefix,
        "promoter_count": len(rows),
        "gene_count": len({row["gene_id"] for row in rows}),
        "output": str(output),
        "output_sha256": sha256_file(output),
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def read_promoters(path: Path) -> list[Promoter]:
    promoters = []
    seen: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(PROMOTER_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} lacks promoter fields: {sorted(missing)}")
        for row in reader:
            promoter_id = row["promoter_id"]
            if promoter_id in seen:
                raise ValueError(f"Duplicate promoter ID {promoter_id}")
            seen.add(promoter_id)
            promoters.append(
                Promoter(
                    promoter_id=promoter_id,
                    gene_id=row["gene_id"],
                    gene_name=row["gene_name"],
                    chrom=row["chrom"],
                    tss=int(row["tss"]),
                    strand=row["strand"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    transcript_ids=row["transcript_ids"],
                )
            )
    if not promoters:
        raise ValueError(f"No promoters in {path}")
    return promoters


def read_context_elements(path: Path, context: str) -> list[Element]:
    required = {
        "blacklist_overlap",
        "master_dhs_id",
        "chrom",
        "start",
        "end",
        "summit",
        "context",
        "context_membership",
        "regulatory_class",
        "atac_normalized_cpm_per_kb",
        "h3k27ac_max_500_background_percentile",
        "h3k27ac_active",
        "activity_state",
        "combined_activity_max_500",
    }
    elements = []
    seen: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        missing = required - fieldnames
        if missing:
            raise ValueError(f"{path} lacks catalog fields: {sorted(missing)}")
        for row in reader:
            if row["context"] != context or row["context_membership"] != "1":
                continue
            master_id = row["master_dhs_id"]
            if master_id in seen:
                raise ValueError(f"{path} repeats {master_id}")
            seen.add(master_id)
            percentile_text = row["h3k27ac_max_500_background_percentile"]
            percentile = float(percentile_text) if percentile_text else None
            elements.append(
                Element(
                    master_dhs_id=master_id,
                    chrom=row["chrom"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    summit=int(row["summit"]),
                    regulatory_class=row["regulatory_class"],
                    activity_state=row["activity_state"],
                    atac_signal=float(row["atac_normalized_cpm_per_kb"]),
                    h3k27ac_background_percentile=percentile,
                    h3k27ac_active=int(row["h3k27ac_active"]),
                    combined_activity=float(row["combined_activity_max_500"]),
                    blacklist_overlap=int(row["blacklist_overlap"]),
                )
            )
    if not elements:
        raise ValueError(f"No context-member elements for {context} in {path}")
    return elements


def promoter_activities(
    promoters: list[Promoter],
    elements: list[Element],
) -> dict[str, PromoterActivity]:
    """Summarize DHSs whose summits lie within 500 bp of each promoter TSS.

    Every supporting DHS is listed, but promoter activity, signal and score use
    only blacklist-clear supporting DHSs.
    """

    by_chromosome: dict[str, list[Element]] = defaultdict(list)
    for element in elements:
        by_chromosome[element.chrom].append(element)
    summits: dict[str, list[int]] = {}
    for chromosome, values in by_chromosome.items():
        values.sort(key=lambda value: (value.summit, value.master_dhs_id))
        summits[chromosome] = [value.summit for value in values]

    result = {}
    for promoter in promoters:
        values = by_chromosome.get(promoter.chrom, [])
        chromosome_summits = summits.get(promoter.chrom, [])
        start = bisect_left(
            chromosome_summits, promoter.tss - PROMOTER_SUPPORT_DISTANCE_BP
        )
        stop = bisect_right(
            chromosome_summits, promoter.tss + PROMOTER_SUPPORT_DISTANCE_BP
        )
        supporting = values[start:stop]
        eligible_supporting = [
            element for element in supporting if not element.blacklist_overlap
        ]
        percentiles = [
            element.h3k27ac_background_percentile
            for element in eligible_supporting
            if element.h3k27ac_background_percentile is not None
        ]
        percentile = max(percentiles) if percentiles else None
        combined = max(
            (element.combined_activity for element in eligible_supporting), default=0.0
        )
        activity_score = max(
            (
                element.combined_activity * element.h3k27ac_background_percentile
                for element in eligible_supporting
                if element.h3k27ac_background_percentile is not None
            ),
            default=0.0,
        )
        result[promoter.promoter_id] = PromoterActivity(
            dhs_ids=";".join(element.master_dhs_id for element in supporting),
            accessible=int(bool(supporting)),
            atac_signal=max(
                (element.atac_signal for element in eligible_supporting), default=0.0
            ),
            h3k27ac_background_percentile=percentile,
            combined_activity=combined,
            activity_score=activity_score,
            active=int(
                any(element.h3k27ac_active for element in eligible_supporting)
            ),
        )
    return result


BASE_FIELDS = (
    "master_dhs_id",
    "element_chrom",
    "element_start",
    "element_end",
    "element_summit",
    "regulatory_class",
    "blacklist_overlap",
    "nearest_tss_distance_bp",
    "nearest_tss_tie_count",
    "nearest_promoter_ids",
    "nearest_gene_ids",
    "nearest_gene_names",
    "nearest_tss_positions",
    "nearest_tss_strands",
    "nearest_promoter_transcript_ids",
    "evidence_type",
    "nearest_promoter_active_contexts",
    "nearest_promoter_active_context_count",
)

CONTEXT_FIELDS = (
    "element_context_membership",
    "element_atac_signal",
    "element_h3k27ac_background_percentile",
    "element_h3k27ac_active",
    "element_activity_state",
    "element_combined_activity",
    "nearest_promoter_active",
    "active_nearest_promoter_ids",
    "active_nearest_promoter_supporting_element_ids",
    "nearest_promoter_activity_score_max",
    "nearest_promoter_atac_signal_max",
    "nearest_promoter_h3k27ac_background_percentile_max",
    "nearest_promoter_combined_activity_max",
)

CATALOG_STORED_FIELDS = (
    "master_dhs_id",
    "chrom",
    "start",
    "end",
    "summit",
    "context",
    "context_membership",
    "regulatory_class",
    "blacklist_overlap",
    "atac_normalized_cpm_per_kb",
    "h3k27ac_max_500_background_percentile",
    "h3k27ac_active",
    "activity_state",
    "combined_activity_max_500",
)
CATALOG_REQUIRED_FIELDS = set(CATALOG_STORED_FIELDS)
STATIC_CATALOG_FIELDS = (
    "chrom",
    "start",
    "end",
    "summit",
    "regulatory_class",
    "blacklist_overlap",
)


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def _catalog_rows(
    path: Path, contexts: list[str]
) -> tuple[list[str], dict[str, dict[str, dict[str, str]]]]:
    by_element: dict[str, dict[str, dict[str, str]]] = {}
    element_order: list[str] = []
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = CATALOG_REQUIRED_FIELDS - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} lacks catalog fields: {sorted(missing)}")
        for row in reader:
            context = row["context"]
            if context not in contexts:
                raise ValueError(f"{path} contains unconfigured context {context!r}")
            master_id = row["master_dhs_id"]
            if master_id not in by_element:
                by_element[master_id] = {}
                element_order.append(master_id)
            if context in by_element[master_id]:
                raise ValueError(f"{path} repeats {master_id}/{context}")
            by_element[master_id][context] = {
                field: row[field] for field in CATALOG_STORED_FIELDS
            }
    expected = set(contexts)
    for master_id, rows in by_element.items():
        if set(rows) != expected:
            raise ValueError(f"{path} lacks contexts for {master_id}")
        first = rows[contexts[0]]
        static = tuple(first[field] for field in STATIC_CATALOG_FIELDS)
        if any(
            tuple(row[field] for field in STATIC_CATALOG_FIELDS)
            != static
            for row in rows.values()
        ):
            raise ValueError(f"{path} has context-dependent coordinates/class for {master_id}")
    if not by_element:
        raise ValueError(f"{path} is empty")
    return element_order, by_element


def _promoter_index(
    promoters: list[Promoter],
) -> tuple[dict[str, list[int]], dict[str, dict[int, list[Promoter]]]]:
    by_position: dict[str, dict[int, list[Promoter]]] = {}
    for promoter in promoters:
        by_position.setdefault(promoter.chrom, {}).setdefault(promoter.tss, []).append(
            promoter
        )
    positions = {
        chromosome: sorted(chromosome_promoters)
        for chromosome, chromosome_promoters in by_position.items()
    }
    for chromosome_promoters in by_position.values():
        for values in chromosome_promoters.values():
            values.sort(key=lambda promoter: promoter.promoter_id)
    return positions, by_position


def _nearest_promoters(
    *,
    chromosome: str,
    summit: int,
    positions: dict[str, list[int]],
    by_position: dict[str, dict[int, list[Promoter]]],
) -> tuple[int | None, list[Promoter]]:
    chromosome_positions = positions.get(chromosome, [])
    if not chromosome_positions:
        return None, []
    index = bisect_left(chromosome_positions, summit)
    candidate_positions = []
    if index < len(chromosome_positions):
        candidate_positions.append(chromosome_positions[index])
    if index:
        candidate_positions.append(chromosome_positions[index - 1])
    distance = min(abs(position - summit) for position in candidate_positions)
    tied_positions = sorted(
        position
        for position in set(candidate_positions)
        if abs(position - summit) == distance
    )
    tied = [
        promoter
        for position in tied_positions
        for promoter in by_position[chromosome][position]
    ]
    return distance, tied


def build_nearest_tss_enhancer_candidates(
    *,
    catalog_path: Path,
    promoters_path: Path,
    contexts: list[str],
    enhancer_classes: list[str],
    output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Write one row per enhancer with tied nearest TSSs and context activity."""

    if not contexts or len(contexts) != len(set(contexts)):
        raise ValueError("Contexts must be non-empty and unique")
    if not enhancer_classes or len(enhancer_classes) != len(set(enhancer_classes)):
        raise ValueError("Enhancer classes must be non-empty and unique")

    element_order, catalog = _catalog_rows(catalog_path, contexts)
    promoters = read_promoters(promoters_path)
    positions, by_position = _promoter_index(promoters)
    activities = {
        context: promoter_activities(
            promoters,
            read_context_elements(catalog_path, context),
        )
        for context in contexts
    }
    fields = [
        *BASE_FIELDS,
        *[
            f"{context}__{field}"
            for context in contexts
            for field in CONTEXT_FIELDS
        ],
    ]
    output_rows: list[dict[str, Any]] = []
    enhancer_class_set = set(enhancer_classes)
    no_same_chromosome_promoter_count = 0
    tie_element_count = 0
    active_by_context = {context: 0 for context in contexts}

    for master_id in element_order:
        context_rows = catalog[master_id]
        first = context_rows[contexts[0]]
        if first["regulatory_class"] not in enhancer_class_set:
            continue
        distance, nearest = _nearest_promoters(
            chromosome=first["chrom"],
            summit=int(first["summit"]),
            positions=positions,
            by_position=by_position,
        )
        if not nearest:
            no_same_chromosome_promoter_count += 1
        if len(nearest) > 1:
            tie_element_count += 1
        active_contexts = []
        context_values: dict[str, Any] = {}
        for context in contexts:
            element = context_rows[context]
            nearest_activities = [
                (promoter, activities[context][promoter.promoter_id])
                for promoter in nearest
            ]
            active = [
                (promoter, activity)
                for promoter, activity in nearest_activities
                if activity.active
            ]
            if active:
                active_contexts.append(context)
                active_by_context[context] += 1
            active_support = sorted(
                {
                    master
                    for _promoter, activity in active
                    for master in activity.dhs_ids.split(";")
                    if master
                }
            )
            percentile_values = [
                activity.h3k27ac_background_percentile
                for _promoter, activity in nearest_activities
                if activity.h3k27ac_background_percentile is not None
            ]
            prefix = f"{context}__"
            context_values.update(
                {
                    f"{prefix}element_context_membership": element[
                        "context_membership"
                    ],
                    f"{prefix}element_atac_signal": element[
                        "atac_normalized_cpm_per_kb"
                    ],
                    f"{prefix}element_h3k27ac_background_percentile": element[
                        "h3k27ac_max_500_background_percentile"
                    ],
                    f"{prefix}element_h3k27ac_active": element["h3k27ac_active"],
                    f"{prefix}element_activity_state": element["activity_state"],
                    f"{prefix}element_combined_activity": element[
                        "combined_activity_max_500"
                    ],
                    f"{prefix}nearest_promoter_active": int(bool(active)),
                    f"{prefix}active_nearest_promoter_ids": ";".join(
                        promoter.promoter_id for promoter, _activity in active
                    ),
                    f"{prefix}active_nearest_promoter_supporting_element_ids": ";".join(
                        active_support
                    ),
                    f"{prefix}nearest_promoter_activity_score_max": _format(
                        max(
                            (
                                activity.activity_score
                                for _promoter, activity in nearest_activities
                            ),
                        )
                        if nearest_activities
                        else None
                    ),
                    f"{prefix}nearest_promoter_atac_signal_max": _format(
                        max(
                            (
                                activity.atac_signal
                                for _promoter, activity in nearest_activities
                            ),
                        )
                        if nearest_activities
                        else None
                    ),
                    f"{prefix}nearest_promoter_h3k27ac_background_percentile_max": _format(
                        max(percentile_values) if percentile_values else None
                    ),
                    f"{prefix}nearest_promoter_combined_activity_max": _format(
                        max(
                            (
                                activity.combined_activity
                                for _promoter, activity in nearest_activities
                            ),
                        )
                        if nearest_activities
                        else None
                    ),
                }
            )
        output_rows.append(
            {
                "master_dhs_id": master_id,
                "element_chrom": first["chrom"],
                "element_start": first["start"],
                "element_end": first["end"],
                "element_summit": first["summit"],
                "regulatory_class": first["regulatory_class"],
                "blacklist_overlap": first["blacklist_overlap"],
                "nearest_tss_distance_bp": "" if distance is None else distance,
                "nearest_tss_tie_count": len(nearest),
                "nearest_promoter_ids": ";".join(
                    promoter.promoter_id for promoter in nearest
                ),
                "nearest_gene_ids": ";".join(promoter.gene_id for promoter in nearest),
                "nearest_gene_names": ";".join(
                    promoter.gene_name for promoter in nearest
                ),
                "nearest_tss_positions": ";".join(
                    str(promoter.tss) for promoter in nearest
                ),
                "nearest_tss_strands": ";".join(
                    promoter.strand for promoter in nearest
                ),
                "nearest_promoter_transcript_ids": "|".join(
                    promoter.transcript_ids for promoter in nearest
                ),
                "evidence_type": "nearest_annotated_tss",
                "nearest_promoter_active_contexts": ";".join(active_contexts),
                "nearest_promoter_active_context_count": len(active_contexts),
                **context_values,
            }
        )

    _atomic_tsv(output, fields, output_rows)
    metrics = {
        "schema_version": 1,
        "method": "nearest_annotated_tss_all_ties_wide_v2",
        "catalog": str(catalog_path),
        "catalog_sha256": sha256_file(catalog_path),
        "promoters": str(promoters_path),
        "promoters_sha256": sha256_file(promoters_path),
        "contexts": contexts,
        "enhancer_classes": enhancer_classes,
        "promoter_activity_definition": (
            "any blacklist-clear context-member master DHS summit within "
            "inclusive +/-500 bp of the exact annotated TSS with catalog "
            "H3K27ac activity"
        ),
        "candidate_assignment": (
            "minimum absolute enhancer-summit to annotated-TSS distance; all exact "
            "distance ties retained in the same enhancer row"
        ),
        "enhancer_count": len(output_rows),
        "nearest_tss_tie_element_count": tie_element_count,
        "no_same_chromosome_promoter_count": no_same_chromosome_promoter_count,
        "enhancer_with_active_nearest_promoter_by_context": active_by_context,
        "output": str(output),
        "output_sha256": sha256_file(output),
    }
    write_json_atomic(metrics_output, metrics)
    return metrics
