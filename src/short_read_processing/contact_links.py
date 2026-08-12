"""Build context-specific element-promoter graphs and gene-candidate tables."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from contextlib import contextmanager
import csv
from dataclasses import dataclass
import gzip
import io
from itertools import groupby
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Iterable

from .artifacts import sha256_file
from .contact_metadata import verify_reported_checksum
from .contacts import canonical_chromosome, write_json_atomic


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

EDGE_FIELDS = (
    "context",
    "source_node_id",
    "target_node_id",
    "master_dhs_id",
    "promoter_id",
    "gene_id",
    "gene_name",
    "distance_bp",
    "promoter_active",
    "promoter_activity_score",
    "promoter_atac_accessible",
    "promoter_atac_signal",
    "promoter_h3k27ac_background_percentile",
    "contact_strategy",
    "contact_assay",
    "contact_match",
    "contact_resolution_bp",
    "contact_status",
    "same_contact_bin",
    "observed_balanced_contact",
    "powerlaw_expected_contact",
    "contact_weight",
    "observed_over_expected",
    "candidate_link_score",
)

GENE_FIELDS = (
    "context",
    "master_dhs_id",
    "element_chrom",
    "element_start",
    "element_end",
    "element_summit",
    "regulatory_class",
    "element_h3k27ac_background_percentile",
    "element_combined_activity",
    "element_blacklist_overlap",
    "gene_id",
    "gene_name",
    "promoter_count",
    "active_promoter_count",
    "best_promoter_id",
    "best_active_promoter_id",
    "minimum_distance_bp",
    "best_promoter_active",
    "best_promoter_atac_signal",
    "best_promoter_h3k27ac_background_percentile",
    "contact_strategy",
    "contact_assay",
    "contact_match",
    "contact_resolution_bp",
    "best_contact_status",
    "best_observed_balanced_contact",
    "best_powerlaw_expected_contact",
    "best_contact_weight",
    "best_observed_over_expected",
    "candidate_link_score",
    "candidate_gene_rank",
    "active_candidate_gene_rank",
    "contact_gene_rank",
    "nearest_gene_rank",
)

DISTANCE_MODEL_GENE_FIELDS = GENE_FIELDS + (
    "evidence_type",
    "is_primary_candidate",
    "nearest_active_tss_gene_ids",
    "nearest_active_tss_gene_names",
    "nearest_active_tss_distance_bp",
    "nearest_tss_gene_ids",
    "nearest_tss_gene_names",
    "nearest_tss_distance_bp",
)

NEAREST_ACTIVE_PROMOTER_GENE_FIELDS = (
    "context",
    "master_dhs_id",
    "element_chrom",
    "element_start",
    "element_end",
    "element_summit",
    "regulatory_class",
    "element_h3k27ac_background_percentile",
    "element_combined_activity",
    "element_blacklist_overlap",
    "promoter_id",
    "active_promoter_supporting_element_ids",
    "active_promoter_supporting_element_count",
    "active_promoter_associated_element_ids",
    "active_promoter_associated_element_count",
    "gene_id",
    "gene_name",
    "promoter_tss",
    "promoter_start",
    "promoter_end",
    "distance_bp",
    "nearest_active_tss_tie_count",
    "promoter_atac_signal",
    "promoter_h3k27ac_background_percentile",
    "promoter_combined_activity",
    "evidence_type",
)

ENHANCER_REGULATORY_CLASSES = {
    "proximal_enhancer_like",
    "distal_enhancer_like",
}

NODE_FIELDS = (
    "context",
    "node_id",
    "node_type",
    "master_dhs_id",
    "promoter_id",
    "gene_id",
    "gene_name",
    "chrom",
    "start",
    "end",
    "anchor",
    "regulatory_class",
    "atac_accessible",
    "atac_signal",
    "h3k27ac_background_percentile",
    "combined_activity",
    "active",
    "blacklist_overlap",
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


@dataclass(frozen=True)
class ContactEvidence:
    status: str
    same_bin: int
    observed: float | None
    expected: float
    weight: float
    enrichment: float | None


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


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


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


def read_context_elements(path: Path, context: str) -> tuple[list[Element], bool]:
    required = {
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
        blacklist_overlap_defaulted = "blacklist_overlap" not in fieldnames
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
                    blacklist_overlap=(
                        0
                        if blacklist_overlap_defaulted
                        else int(row["blacklist_overlap"])
                    ),
                )
            )
    if not elements:
        raise ValueError(f"No context-member elements for {context} in {path}")
    return elements, blacklist_overlap_defaulted


def promoter_activities(
    promoters: list[Promoter],
    elements: list[Element],
) -> dict[str, PromoterActivity]:
    """Summarize DHSs whose summits lie within 500 bp of each promoter TSS."""

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
                any(element.h3k27ac_active for element in supporting)
            ),
        )
    return result


def powerlaw_expected(distance: int, gamma: float, scale: float, minimum: int) -> float:
    return scale * max(distance, minimum) ** gamma


class ObservedContact:
    """Return balanced contact and explicit distance-model components."""

    def __init__(
        self,
        path: Path,
        *,
        gamma: float,
        scale: float,
        pseudocount_fraction: float,
    ):
        try:
            import cooler  # type: ignore
        except ImportError as error:
            raise RuntimeError("cooler is required for observed contact links") from error
        self.cooler = cooler.Cooler(str(path))
        if self.cooler.binsize is None:
            raise ValueError(f"Variable-bin Cooler is unsupported: {path}")
        self.resolution = int(self.cooler.binsize)
        self.gamma = gamma
        self.scale = scale
        self.pseudocount_fraction = pseudocount_fraction
        self.aliases = {
            canonical_chromosome(chromosome): chromosome
            for chromosome in self.cooler.chromnames
        }
        self._cached_chromosome = ""
        self._cached_matrix = None

    def __call__(
        self, chromosome: str, first_position: int, second_position: int
    ) -> ContactEvidence:
        distance = abs(first_position - second_position)
        expected = powerlaw_expected(
            distance, self.gamma, self.scale, self.resolution
        )
        source_chromosome = self.aliases.get(canonical_chromosome(chromosome))
        same_bin = int(first_position // self.resolution == second_position // self.resolution)
        if source_chromosome is None:
            return ContactEvidence(
                "chromosome_absent_from_matrix",
                same_bin,
                None,
                expected,
                self.pseudocount_fraction * expected,
                None,
            )
        if self._cached_chromosome != source_chromosome:
            self._cached_matrix = self.cooler.matrix(
                balance=True, sparse=True
            ).fetch(source_chromosome).tocsr()
            self._cached_chromosome = source_chromosome
        matrix = self._cached_matrix
        first_bin = first_position // self.resolution
        second_bin = second_position // self.resolution
        if first_bin >= matrix.shape[0] or second_bin >= matrix.shape[1]:
            return ContactEvidence(
                "position_outside_matrix",
                same_bin,
                None,
                expected,
                self.pseudocount_fraction * expected,
                None,
            )
        if same_bin:
            neighbors: list[float] = []
            if first_bin > 0:
                neighbors.append(float(matrix[first_bin - 1, first_bin]))
            if first_bin + 1 < matrix.shape[0]:
                neighbors.append(float(matrix[first_bin, first_bin + 1]))
            finite = [
                value for value in neighbors if math.isfinite(value) and value >= 0
            ]
            if not finite:
                return ContactEvidence(
                    "unbalanced_or_masked_adjacent_bins",
                    same_bin,
                    None,
                    expected,
                    self.pseudocount_fraction * expected,
                    None,
                )
            observed = max(finite)
            status = "adjacent_bin_proxy"
        else:
            observed = float(matrix[first_bin, second_bin])
            if not math.isfinite(observed) or observed < 0:
                return ContactEvidence(
                    "unbalanced_or_masked_bin",
                    same_bin,
                    None,
                    expected,
                    self.pseudocount_fraction * expected,
                    None,
                )
            status = "observed_matrix_pixel"
        return ContactEvidence(
            status,
            same_bin,
            observed,
            expected,
            observed + self.pseudocount_fraction * expected,
            observed / expected if expected > 0 else None,
        )


def _rank(values: list[dict[str, Any]], key) -> dict[str, int]:
    return {
        row["gene_id"]: index
        for index, row in enumerate(sorted(values, key=key), start=1)
    }


def _element_gene_rows(
    element: Element,
    edge_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in edge_rows:
        by_gene[row["gene_id"]].append(row)
    result = []
    for gene_id, rows in by_gene.items():
        best = max(
            rows,
            key=lambda row: (
                float(row["candidate_link_score"]),
                float(row["contact_weight"]),
                -int(row["distance_bp"]),
                row["promoter_id"],
            ),
        )
        active_rows = [row for row in rows if str(row["promoter_active"]) == "1"]
        best_active = (
            max(
                active_rows,
                key=lambda row: (
                    float(row["candidate_link_score"]),
                    float(row["contact_weight"]),
                    -int(row["distance_bp"]),
                    row["promoter_id"],
                ),
            )
            if active_rows
            else None
        )
        result.append(
            {
                "context": best["context"],
                "master_dhs_id": element.master_dhs_id,
                "element_chrom": element.chrom,
                "element_start": element.start,
                "element_end": element.end,
                "element_summit": element.summit,
                "regulatory_class": element.regulatory_class,
                "element_h3k27ac_background_percentile": _format(
                    element.h3k27ac_background_percentile
                ),
                "element_combined_activity": _format(element.combined_activity),
                "element_blacklist_overlap": element.blacklist_overlap,
                "gene_id": gene_id,
                "gene_name": best["gene_name"],
                "promoter_count": len(rows),
                "active_promoter_count": len(active_rows),
                "best_promoter_id": best["promoter_id"],
                "best_active_promoter_id": (
                    best_active["promoter_id"] if best_active else ""
                ),
                "minimum_distance_bp": min(int(row["distance_bp"]) for row in rows),
                "best_promoter_active": best["promoter_active"],
                "best_promoter_atac_signal": best["promoter_atac_signal"],
                "best_promoter_h3k27ac_background_percentile": best[
                    "promoter_h3k27ac_background_percentile"
                ],
                "contact_strategy": best["contact_strategy"],
                "contact_assay": best["contact_assay"],
                "contact_match": best["contact_match"],
                "contact_resolution_bp": best["contact_resolution_bp"],
                "best_contact_status": best["contact_status"],
                "best_observed_balanced_contact": best[
                    "observed_balanced_contact"
                ],
                "best_powerlaw_expected_contact": best[
                    "powerlaw_expected_contact"
                ],
                "best_contact_weight": best["contact_weight"],
                "best_observed_over_expected": best["observed_over_expected"],
                "candidate_link_score": best["candidate_link_score"],
                "_contact_rank_value": max(float(row["contact_weight"]) for row in rows),
                "_active_score": (
                    float(best_active["candidate_link_score"])
                    if best_active
                    else None
                ),
            }
        )
    candidate_ranks = _rank(
        result,
        lambda row: (
            -float(row["candidate_link_score"]),
            int(row["minimum_distance_bp"]),
            row["gene_id"],
        ),
    )
    contact_ranks = _rank(
        result,
        lambda row: (
            -float(row["_contact_rank_value"]),
            int(row["minimum_distance_bp"]),
            row["gene_id"],
        ),
    )
    nearest_ranks = _rank(
        result,
        lambda row: (int(row["minimum_distance_bp"]), row["gene_id"]),
    )
    active = [row for row in result if row["_active_score"] is not None]
    active_ranks = _rank(
        active,
        lambda row: (
            -float(row["_active_score"]),
            int(row["minimum_distance_bp"]),
            row["gene_id"],
        ),
    )
    for row in result:
        gene_id = row["gene_id"]
        row["candidate_gene_rank"] = candidate_ranks[gene_id]
        row["active_candidate_gene_rank"] = active_ranks.get(gene_id, "")
        row["contact_gene_rank"] = contact_ranks[gene_id]
        row["nearest_gene_rank"] = nearest_ranks[gene_id]
        row.pop("_contact_rank_value")
        row.pop("_active_score")
    return sorted(result, key=lambda row: int(row["candidate_gene_rank"]))


def _active_contact_enhancer_gene_rows(
    element: Element,
    edge_rows: list[dict[str, Any]],
    *,
    observed_over_expected_threshold: float,
) -> list[dict[str, Any]]:
    """Project threshold-qualified observed enhancer contacts to genes."""

    if (
        element.regulatory_class not in ENHANCER_REGULATORY_CLASSES
        or not element.h3k27ac_active
    ):
        return []
    qualifying_edges = [
        row
        for row in edge_rows
        if row["observed_over_expected"] != ""
        and float(row["observed_over_expected"])
        >= observed_over_expected_threshold
    ]
    return _element_gene_rows(element, qualifying_edges)


def _focused_candidate_elements(
    *,
    context: str,
    nodes_path: Path,
) -> tuple[dict[str, Element], set[str], int]:
    """Read threshold-qualified enhancer nodes for a focused projection."""

    elements: dict[str, Element] = {}
    all_element_ids: set[str] = set()
    element_node_count = 0
    with _open_text(nodes_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(NODE_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{nodes_path} lacks node fields: {sorted(missing)}")
        for row in reader:
            if row["context"] != context:
                raise ValueError(
                    f"{nodes_path} contains context {row['context']!r}, "
                    f"expected {context!r}"
                )
            if row["node_type"] != "element":
                continue
            element_node_count += 1
            master_dhs_id = row["master_dhs_id"]
            if master_dhs_id in all_element_ids:
                raise ValueError(f"{nodes_path} repeats element {master_dhs_id}")
            all_element_ids.add(master_dhs_id)
            percentile = (
                float(row["h3k27ac_background_percentile"])
                if row["h3k27ac_background_percentile"]
                else None
            )
            element = Element(
                master_dhs_id=master_dhs_id,
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                summit=int(row["anchor"]),
                regulatory_class=row["regulatory_class"],
                activity_state="",
                atac_signal=float(row["atac_signal"]),
                h3k27ac_background_percentile=percentile,
                h3k27ac_active=int(row["active"]),
                combined_activity=float(row["combined_activity"]),
                blacklist_overlap=int(row["blacklist_overlap"]),
            )
            if (
                element.regulatory_class in ENHANCER_REGULATORY_CLASSES
                and element.h3k27ac_active
            ):
                elements[element.master_dhs_id] = element
    return elements, all_element_ids, element_node_count


def build_nearest_active_promoter_gene_candidates(
    *,
    context: str,
    nodes_path: Path,
    output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Assign qualifying enhancers to the nearest supported active TSS."""

    elements, all_element_ids, element_node_count = _focused_candidate_elements(
        context=context,
        nodes_path=nodes_path,
    )

    active_element_ids: set[str] = set()
    active_promoter_associated_element_ids: set[str] = set()
    promoter_rows: list[dict[str, Any]] = []
    seen_promoters: set[str] = set()
    with _open_text(nodes_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(NODE_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{nodes_path} lacks node fields: {sorted(missing)}")
        rows = list(reader)
    for row in rows:
        if row["context"] != context:
            raise ValueError(
                f"{nodes_path} contains context {row['context']!r}, "
                f"expected {context!r}"
            )
        if row["node_type"] == "element" and str(row["active"]) == "1":
            active_element_ids.add(row["master_dhs_id"])
            if row["regulatory_class"] == "promoter_associated":
                active_promoter_associated_element_ids.add(row["master_dhs_id"])

    for row in rows:
        if row["node_type"] != "promoter":
            continue
        promoter_id = row["promoter_id"]
        if promoter_id in seen_promoters:
            raise ValueError(f"{nodes_path} repeats promoter {promoter_id}")
        seen_promoters.add(promoter_id)
        supporting_node_ids = {
            value for value in row["master_dhs_id"].split(";") if value
        }
        unknown = supporting_node_ids - all_element_ids
        if unknown:
            raise ValueError(
                f"{nodes_path} promoter {promoter_id} references unknown elements: "
                f"{sorted(unknown)}"
            )
        if str(row["active"]) != "1":
            continue
        supporting_ids = sorted(supporting_node_ids & active_element_ids)
        if not supporting_ids:
            raise ValueError(
                f"{nodes_path} active promoter {promoter_id} lacks an active "
                "supporting element"
            )
        promoter_associated_ids = sorted(
            supporting_node_ids & active_promoter_associated_element_ids
        )
        promoter_rows.append(
            row
            | {
                "active_promoter_supporting_element_ids": ";".join(
                    supporting_ids
                ),
                "active_promoter_supporting_element_count": len(supporting_ids),
                "active_promoter_associated_element_ids": ";".join(
                    promoter_associated_ids
                ),
                "active_promoter_associated_element_count": len(
                    promoter_associated_ids
                ),
            }
        )

    promoters_by_chromosome: dict[str, dict[int, list[dict[str, Any]]]] = (
        defaultdict(lambda: defaultdict(list))
    )
    for row in promoter_rows:
        promoters_by_chromosome[row["chrom"]][int(row["anchor"])].append(row)
    promoter_tss = {
        chromosome: sorted(by_tss)
        for chromosome, by_tss in promoters_by_chromosome.items()
    }
    for by_tss in promoters_by_chromosome.values():
        for values in by_tss.values():
            values.sort(key=lambda row: row["promoter_id"])

    candidate_count = 0
    candidate_element_ids: set[str] = set()

    def candidate_rows():
        nonlocal candidate_count
        for element in elements.values():
            coordinates = promoter_tss.get(element.chrom, [])
            if not coordinates:
                continue
            index = bisect_left(coordinates, element.summit)
            adjacent = []
            if index < len(coordinates):
                adjacent.append(coordinates[index])
            if index:
                adjacent.append(coordinates[index - 1])
            minimum_distance = min(
                abs(coordinate - element.summit) for coordinate in adjacent
            )
            nearest_coordinates = sorted(
                coordinate
                for coordinate in set(adjacent)
                if abs(coordinate - element.summit) == minimum_distance
            )
            nearest_promoters = [
                row
                for coordinate in nearest_coordinates
                for row in promoters_by_chromosome[element.chrom][coordinate]
            ]
            tie_count = len(nearest_promoters)
            candidate_element_ids.add(element.master_dhs_id)
            candidate_count += tie_count
            for promoter in nearest_promoters:
                yield {
                    "context": context,
                    "master_dhs_id": element.master_dhs_id,
                    "element_chrom": element.chrom,
                    "element_start": element.start,
                    "element_end": element.end,
                    "element_summit": element.summit,
                    "regulatory_class": element.regulatory_class,
                    "element_h3k27ac_background_percentile": _format(
                        element.h3k27ac_background_percentile
                    ),
                    "element_combined_activity": _format(
                        element.combined_activity
                    ),
                    "element_blacklist_overlap": element.blacklist_overlap,
                    "promoter_id": promoter["promoter_id"],
                    "active_promoter_supporting_element_ids": promoter[
                        "active_promoter_supporting_element_ids"
                    ],
                    "active_promoter_supporting_element_count": promoter[
                        "active_promoter_supporting_element_count"
                    ],
                    "active_promoter_associated_element_ids": promoter[
                        "active_promoter_associated_element_ids"
                    ],
                    "active_promoter_associated_element_count": promoter[
                        "active_promoter_associated_element_count"
                    ],
                    "gene_id": promoter["gene_id"],
                    "gene_name": promoter["gene_name"],
                    "promoter_tss": promoter["anchor"],
                    "promoter_start": promoter["start"],
                    "promoter_end": promoter["end"],
                    "distance_bp": minimum_distance,
                    "nearest_active_tss_tie_count": tie_count,
                    "promoter_atac_signal": promoter["atac_signal"],
                    "promoter_h3k27ac_background_percentile": promoter[
                        "h3k27ac_background_percentile"
                    ],
                    "promoter_combined_activity": promoter["combined_activity"],
                    "evidence_type": "nearest_active_promoter_tss",
                }

    _atomic_tsv(output, NEAREST_ACTIVE_PROMOTER_GENE_FIELDS, candidate_rows())
    metrics = {
        "schema_version": 1,
        "context": context,
        "method": "nearest_active_promoter_gene_candidates_v1",
        "evidence_type": "nearest_active_promoter_tss",
        "promoter_activity_definition": (
            "active flag inherited from the context promoter node"
        ),
        "selection": (
            "proximal_or_distal_enhancer_like and catalog H3K27ac active; "
            "nearest same-chromosome active promoter TSS"
        ),
        "distance_limit_bp": None,
        "exact_distance_ties_retained": True,
        "contact_or_powerlaw_used": False,
        "element_node_count": element_node_count,
        "qualifying_element_count": len(elements),
        "active_element_count": len(active_element_ids),
        "active_promoter_associated_element_count": len(
            active_promoter_associated_element_ids
        ),
        "active_promoter_tss_count": len(promoter_rows),
        "element_with_nearest_active_promoter_count": len(candidate_element_ids),
        "qualifying_element_without_active_promoter_on_chromosome_count": (
            len(elements) - len(candidate_element_ids)
        ),
        "nearest_active_promoter_gene_candidate_count": candidate_count,
        "inputs": {
            "nodes": {"path": str(nodes_path), "sha256": sha256_file(nodes_path)},
        },
        "output": {"path": str(output), "sha256": sha256_file(output)},
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def build_active_contact_enhancer_gene_candidates(
    *,
    context: str,
    nodes_path: Path,
    edges_path: Path,
    observed_over_expected_threshold: float,
    output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Derive the focused enhancer--gene table from completed graph tables."""

    if observed_over_expected_threshold < 0:
        raise ValueError("Candidate observed/expected threshold cannot be negative")

    elements, all_element_ids, element_node_count = _focused_candidate_elements(
        context=context,
        nodes_path=nodes_path,
    )

    candidate_count = 0
    candidate_element_ids: set[str] = set()

    def candidate_rows():
        nonlocal candidate_count
        seen: set[str] = set()
        with _open_text(edges_path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(EDGE_FIELDS) - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{edges_path} lacks edge fields: {sorted(missing)}")
            for master_dhs_id, grouped_rows in groupby(
                reader, key=lambda row: row["master_dhs_id"]
            ):
                if master_dhs_id in seen:
                    raise ValueError(
                        f"{edges_path} has non-contiguous edges for {master_dhs_id}"
                    )
                seen.add(master_dhs_id)
                edge_rows = list(grouped_rows)
                if any(row["context"] != context for row in edge_rows):
                    raise ValueError(
                        f"{edges_path} contains an edge outside context {context!r}"
                    )
                if master_dhs_id not in all_element_ids:
                    raise ValueError(
                        f"{edges_path} references unknown element {master_dhs_id}"
                    )
                element = elements.get(master_dhs_id)
                if element is None:
                    continue
                rows = _active_contact_enhancer_gene_rows(
                    element,
                    edge_rows,
                    observed_over_expected_threshold=(
                        observed_over_expected_threshold
                    ),
                )
                if rows:
                    candidate_element_ids.add(master_dhs_id)
                candidate_count += len(rows)
                yield from rows

    _atomic_tsv(output, GENE_FIELDS, candidate_rows())
    metrics = {
        "schema_version": 1,
        "context": context,
        "method": "active_contact_enhancer_gene_candidates_v1",
        "observed_over_expected_threshold": observed_over_expected_threshold,
        "selection": (
            "proximal_or_distal_enhancer_like and catalog H3K27ac active and "
            "observed_over_expected >= contact threshold"
        ),
        "element_node_count": element_node_count,
        "qualifying_element_count": len(elements),
        "element_with_candidate_count": len(candidate_element_ids),
        "active_contact_enhancer_gene_candidate_count": candidate_count,
        "inputs": {
            "nodes": {"path": str(nodes_path), "sha256": sha256_file(nodes_path)},
            "element_promoter_edges": {
                "path": str(edges_path),
                "sha256": sha256_file(edges_path),
            },
        },
        "output": {"path": str(output), "sha256": sha256_file(output)},
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def _nearest_tss_summary(
    edge_rows: list[dict[str, Any]],
) -> tuple[str, str, str]:
    if not edge_rows:
        return "", "", ""
    minimum_distance = min(int(row["distance_bp"]) for row in edge_rows)
    nearest = sorted(
        {
            (row["gene_id"], row["gene_name"])
            for row in edge_rows
            if int(row["distance_bp"]) == minimum_distance
        }
    )
    return (
        ";".join(gene_id for gene_id, _gene_name in nearest),
        ";".join(gene_name for _gene_id, gene_name in nearest),
        str(minimum_distance),
    )


def _active_distance_enhancer_gene_rows(
    element: Element,
    edge_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Project active-promoter distance-model edges to ranked genes."""

    if (
        element.regulatory_class not in ENHANCER_REGULATORY_CLASSES
        or not element.h3k27ac_active
    ):
        return []
    if any(row["contact_strategy"] != "powerlaw" for row in edge_rows):
        raise ValueError("Distance candidates require powerlaw contact edges")
    if any(
        row["observed_balanced_contact"] != ""
        or row["observed_over_expected"] != ""
        for row in edge_rows
    ):
        raise ValueError("Distance candidates cannot contain observed contact values")

    active_edges = [row for row in edge_rows if str(row["promoter_active"]) == "1"]
    nearest_active = _nearest_tss_summary(active_edges)
    nearest_all = _nearest_tss_summary(edge_rows)
    rows = _element_gene_rows(element, active_edges)
    for row in rows:
        row.update(
            {
                "evidence_type": "distance_model_active_promoter",
                "is_primary_candidate": int(row["candidate_gene_rank"] == 1),
                "nearest_active_tss_gene_ids": nearest_active[0],
                "nearest_active_tss_gene_names": nearest_active[1],
                "nearest_active_tss_distance_bp": nearest_active[2],
                "nearest_tss_gene_ids": nearest_all[0],
                "nearest_tss_gene_names": nearest_all[1],
                "nearest_tss_distance_bp": nearest_all[2],
            }
        )
    return rows


def build_active_distance_enhancer_gene_candidates(
    *,
    context: str,
    nodes_path: Path,
    edges_path: Path,
    output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Derive ranked active-promoter candidates for a distance-only context."""

    elements, all_element_ids, element_node_count = _focused_candidate_elements(
        context=context,
        nodes_path=nodes_path,
    )

    candidate_count = 0
    candidate_element_ids: set[str] = set()

    def candidate_rows():
        nonlocal candidate_count
        seen: set[str] = set()
        with _open_text(edges_path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(EDGE_FIELDS) - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{edges_path} lacks edge fields: {sorted(missing)}")
            for master_dhs_id, grouped_rows in groupby(
                reader, key=lambda row: row["master_dhs_id"]
            ):
                if master_dhs_id in seen:
                    raise ValueError(
                        f"{edges_path} has non-contiguous edges for {master_dhs_id}"
                    )
                seen.add(master_dhs_id)
                edge_rows = list(grouped_rows)
                if any(row["context"] != context for row in edge_rows):
                    raise ValueError(
                        f"{edges_path} contains an edge outside context {context!r}"
                    )
                if master_dhs_id not in all_element_ids:
                    raise ValueError(
                        f"{edges_path} references unknown element {master_dhs_id}"
                    )
                element = elements.get(master_dhs_id)
                if element is None:
                    continue
                rows = _active_distance_enhancer_gene_rows(
                    element,
                    edge_rows,
                )
                if rows:
                    candidate_element_ids.add(master_dhs_id)
                candidate_count += len(rows)
                yield from rows

    _atomic_tsv(output, DISTANCE_MODEL_GENE_FIELDS, candidate_rows())
    metrics = {
        "schema_version": 1,
        "context": context,
        "method": "active_distance_enhancer_gene_candidates_v1",
        "evidence_type": "distance_model_active_promoter",
        "promoter_activity_required": True,
        "selection": (
            "proximal_or_distal_enhancer_like and catalog H3K27ac active and "
            "promoter_active = 1"
        ),
        "ranking": "powerlaw contact weight * promoter activity score",
        "nearest_tss_baselines": (
            "closest active promoter TSS and closest annotated promoter TSS "
            "within the upstream maximum-distance window"
        ),
        "observed_contact_values": "not_applicable",
        "element_node_count": element_node_count,
        "qualifying_element_count": len(elements),
        "element_with_candidate_count": len(candidate_element_ids),
        "qualifying_element_without_active_promoter_count": (
            len(elements) - len(candidate_element_ids)
        ),
        "active_distance_enhancer_gene_candidate_count": candidate_count,
        "inputs": {
            "nodes": {"path": str(nodes_path), "sha256": sha256_file(nodes_path)},
            "element_promoter_edges": {
                "path": str(edges_path),
                "sha256": sha256_file(edges_path),
            },
        },
        "output": {"path": str(output), "sha256": sha256_file(output)},
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def build_context_links(
    *,
    context: str,
    context_elements: Path,
    promoters_path: Path,
    powerlaw_path: Path,
    strategy: str,
    contact_assay: str,
    contact_match: str,
    configured_resolution: int,
    maximum_distance: int,
    pseudocount_fraction: float,
    observed_contact_path: Path | None,
    nodes_output: Path,
    edges_output: Path,
    genes_output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Write one weighted graph layer and its element-to-gene projection."""

    if strategy not in {"observed", "powerlaw"}:
        raise ValueError(f"Unsupported contact strategy: {strategy}")
    if maximum_distance < 1 or configured_resolution < 1:
        raise ValueError("Contact distance and resolution must be positive")
    if not 0 <= pseudocount_fraction <= 1:
        raise ValueError("Contact pseudocount fraction must be in [0, 1]")
    promoters = read_promoters(promoters_path)
    elements, blacklist_overlap_defaulted = read_context_elements(
        context_elements, context
    )
    activities = promoter_activities(promoters, elements)
    with powerlaw_path.open(encoding="utf-8") as handle:
        powerlaw = json.load(handle)
    if strategy == "observed":
        if observed_contact_path is None:
            raise ValueError(f"{context}: observed strategy requires a contact map")
        fit = powerlaw["contexts"].get(context)
        if fit is None:
            raise ValueError(f"{context}: observed contact decay fit is absent")
        contact_reader = ObservedContact(
            observed_contact_path,
            gamma=float(fit["gamma"]),
            scale=float(fit["scale"]),
            pseudocount_fraction=pseudocount_fraction,
        )
        resolution = contact_reader.resolution
        if resolution != configured_resolution:
            raise ValueError(
                f"{context}: contact resolution {resolution} differs from configured "
                f"{configured_resolution}"
            )
        gamma = float(fit["gamma"])
        scale = float(fit["scale"])
    else:
        if observed_contact_path is not None:
            raise ValueError(f"{context}: powerlaw strategy cannot receive a contact map")
        atlas_fit = powerlaw["atlas_powerlaw"]
        gamma = float(atlas_fit["gamma"])
        scale = float(atlas_fit["scale"])
        resolution = configured_resolution
        contact_reader = None

    promoters_by_chromosome: dict[str, list[Promoter]] = defaultdict(list)
    promoter_tss: dict[str, list[int]] = {}
    for promoter in promoters:
        promoters_by_chromosome[promoter.chrom].append(promoter)
    for chromosome, values in promoters_by_chromosome.items():
        values.sort(key=lambda promoter: (promoter.tss, promoter.promoter_id))
        promoter_tss[chromosome] = [promoter.tss for promoter in values]

    node_rows: list[dict[str, Any]] = []
    for element in elements:
        node_rows.append(
            {
                "context": context,
                "node_id": f"element:{element.master_dhs_id}",
                "node_type": "element",
                "master_dhs_id": element.master_dhs_id,
                "promoter_id": "",
                "gene_id": "",
                "gene_name": "",
                "chrom": element.chrom,
                "start": element.start,
                "end": element.end,
                "anchor": element.summit,
                "regulatory_class": element.regulatory_class,
                "atac_accessible": 1,
                "atac_signal": _format(element.atac_signal),
                "h3k27ac_background_percentile": _format(
                    element.h3k27ac_background_percentile
                ),
                "combined_activity": _format(element.combined_activity),
                "active": element.h3k27ac_active,
                "blacklist_overlap": element.blacklist_overlap,
            }
        )
    for promoter in promoters:
        activity = activities[promoter.promoter_id]
        node_rows.append(
            {
                "context": context,
                "node_id": f"promoter:{promoter.promoter_id}",
                "node_type": "promoter",
                "master_dhs_id": activity.dhs_ids,
                "promoter_id": promoter.promoter_id,
                "gene_id": promoter.gene_id,
                "gene_name": promoter.gene_name,
                "chrom": promoter.chrom,
                "start": promoter.start,
                "end": promoter.end,
                "anchor": promoter.tss,
                "regulatory_class": "promoter",
                "atac_accessible": activity.accessible,
                "atac_signal": _format(activity.atac_signal),
                "h3k27ac_background_percentile": _format(
                    activity.h3k27ac_background_percentile
                ),
                "combined_activity": _format(activity.combined_activity),
                "active": activity.active,
                "blacklist_overlap": "",
            }
        )
    _atomic_tsv(nodes_output, NODE_FIELDS, node_rows)

    edge_count = 0
    gene_count = 0
    observed_edge_count = 0
    active_edge_count = 0

    def edge_and_gene_rows():
        nonlocal edge_count, gene_count, observed_edge_count, active_edge_count
        for element in elements:
            chromosome_promoters = promoters_by_chromosome.get(element.chrom, [])
            tss_values = promoter_tss.get(element.chrom, [])
            left = bisect_left(tss_values, element.summit - maximum_distance)
            right = bisect_right(tss_values, element.summit + maximum_distance)
            element_edges = []
            for promoter in chromosome_promoters[left:right]:
                distance = abs(element.summit - promoter.tss)
                activity = activities[promoter.promoter_id]
                if contact_reader is None:
                    expected = powerlaw_expected(distance, gamma, scale, resolution)
                    evidence = ContactEvidence(
                        "distance_model",
                        int(element.summit // resolution == promoter.tss // resolution),
                        None,
                        expected,
                        expected,
                        None,
                    )
                else:
                    evidence = contact_reader(
                        element.chrom, element.summit, promoter.tss
                    )
                candidate_score = evidence.weight * activity.activity_score
                row = {
                    "context": context,
                    "source_node_id": f"element:{element.master_dhs_id}",
                    "target_node_id": f"promoter:{promoter.promoter_id}",
                    "master_dhs_id": element.master_dhs_id,
                    "promoter_id": promoter.promoter_id,
                    "gene_id": promoter.gene_id,
                    "gene_name": promoter.gene_name,
                    "distance_bp": distance,
                    "promoter_active": activity.active,
                    "promoter_activity_score": _format(activity.activity_score),
                    "promoter_atac_accessible": activity.accessible,
                    "promoter_atac_signal": _format(activity.atac_signal),
                    "promoter_h3k27ac_background_percentile": _format(
                        activity.h3k27ac_background_percentile
                    ),
                    "contact_strategy": strategy,
                    "contact_assay": contact_assay,
                    "contact_match": contact_match,
                    "contact_resolution_bp": resolution,
                    "contact_status": evidence.status,
                    "same_contact_bin": evidence.same_bin,
                    "observed_balanced_contact": _format(evidence.observed),
                    "powerlaw_expected_contact": _format(evidence.expected),
                    "contact_weight": _format(evidence.weight),
                    "observed_over_expected": _format(evidence.enrichment),
                    "candidate_link_score": _format(candidate_score),
                }
                element_edges.append(row)
                edge_count += 1
                observed_edge_count += int(evidence.observed is not None)
                active_edge_count += activity.active
            for row in element_edges:
                yield "edge", row
            gene_rows = _element_gene_rows(element, element_edges)
            gene_count += len(gene_rows)
            for row in gene_rows:
                yield "gene", row

    # Stream the expensive pair construction once into both atomic outputs.
    edges_output.parent.mkdir(parents=True, exist_ok=True)
    genes_output.parent.mkdir(parents=True, exist_ok=True)
    edge_descriptor, edge_temporary_name = tempfile.mkstemp(
        prefix=f".{edges_output.name}.", suffix=".tmp", dir=edges_output.parent
    )
    gene_descriptor, gene_temporary_name = tempfile.mkstemp(
        prefix=f".{genes_output.name}.", suffix=".tmp", dir=genes_output.parent
    )
    os.close(edge_descriptor)
    os.close(gene_descriptor)
    edge_temporary = Path(edge_temporary_name)
    gene_temporary = Path(gene_temporary_name)
    try:
        with _output_text(edge_temporary, compressed=True) as edge_handle, _output_text(
            gene_temporary, compressed=True
        ) as gene_handle:
            edge_writer = csv.DictWriter(
                edge_handle,
                fieldnames=EDGE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            gene_writer = csv.DictWriter(
                gene_handle,
                fieldnames=GENE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            edge_writer.writeheader()
            gene_writer.writeheader()
            for row_type, row in edge_and_gene_rows():
                if row_type == "edge":
                    edge_writer.writerow(row)
                else:
                    gene_writer.writerow(row)
        os.replace(edge_temporary, edges_output)
        os.replace(gene_temporary, genes_output)
    finally:
        edge_temporary.unlink(missing_ok=True)
        gene_temporary.unlink(missing_ok=True)

    metrics = {
        "schema_version": 1,
        "context": context,
        "contact_strategy": strategy,
        "contact_assay": contact_assay,
        "contact_match": contact_match,
        "contact_resolution_bp": resolution,
        "maximum_distance_bp": maximum_distance,
        "pseudocount_fraction": pseudocount_fraction,
        "blacklist_overlap_annotation": (
            "legacy_catalog_default_zero"
            if blacklist_overlap_defaulted
            else "catalog"
        ),
        "promoter_activity_method": "maximum over context-member master DHSs "
        "whose summits are within 500 bp of the promoter TSS",
        "promoter_support_distance_bp": PROMOTER_SUPPORT_DISTANCE_BP,
        "promoter_support_distance_inclusive": True,
        "link_score": "contact_weight * maximum over blacklist-clear promoter-supporting "
        "DHSs of (combined_activity * h3k27ac_background_percentile)",
        "element_node_count": len(elements),
        "promoter_node_count": len(promoters),
        "active_promoter_count": sum(activity.active for activity in activities.values()),
        "element_promoter_edge_count": edge_count,
        "element_gene_candidate_count": gene_count,
        "observed_edge_count": observed_edge_count,
        "active_promoter_edge_count": active_edge_count,
        "powerlaw_gamma": gamma,
        "powerlaw_scale": scale,
        "inputs": {
            "context_elements": {
                "path": str(context_elements),
                "sha256": sha256_file(context_elements),
            },
            "promoters": {
                "path": str(promoters_path),
                "sha256": sha256_file(promoters_path),
            },
            "powerlaw": {
                "path": str(powerlaw_path),
                "sha256": sha256_file(powerlaw_path),
            },
            "observed_contact": (
                {
                    "path": str(observed_contact_path),
                    "sha256": sha256_file(observed_contact_path),
                }
                if observed_contact_path is not None
                else None
            ),
        },
        "outputs": {
            "nodes": {"path": str(nodes_output), "sha256": sha256_file(nodes_output)},
            "element_promoter_edges": {
                "path": str(edges_output),
                "sha256": sha256_file(edges_output),
            },
            "element_gene_candidates": {
                "path": str(genes_output),
                "sha256": sha256_file(genes_output),
            },
        },
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def aggregate_link_metrics(
    *,
    context_metric_paths: list[Path],
    candidate_metric_paths: list[Path],
    nearest_candidate_metric_paths: list[Path],
    distance_candidate_metric_paths: list[Path],
    source_manifest: Path,
    promoter_metrics: Path,
    contact_metrics: list[Path],
    powerlaw: Path,
    output_metrics: Path,
    output_provenance: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    contexts = []
    for path in context_metric_paths:
        with path.open(encoding="utf-8") as handle:
            contexts.append(json.load(handle))
    contexts.sort(key=lambda row: row["context"])
    candidate_metrics = []
    for path in candidate_metric_paths:
        with path.open(encoding="utf-8") as handle:
            candidate_metrics.append(json.load(handle))
    candidates_by_context = {
        row["context"]: row for row in candidate_metrics
    }
    if len(candidates_by_context) != len(candidate_metrics):
        raise ValueError("Focused candidate metrics repeat a context")
    if set(candidates_by_context) != {row["context"] for row in contexts}:
        raise ValueError("Focused candidate metrics do not match graph contexts")
    nearest_candidate_metrics = []
    for path in nearest_candidate_metric_paths:
        with path.open(encoding="utf-8") as handle:
            nearest_candidate_metrics.append(json.load(handle))
    nearest_candidates_by_context = {
        row["context"]: row for row in nearest_candidate_metrics
    }
    if len(nearest_candidates_by_context) != len(nearest_candidate_metrics):
        raise ValueError("Nearest-promoter candidate metrics repeat a context")
    if set(nearest_candidates_by_context) != {
        row["context"] for row in contexts
    }:
        raise ValueError(
            "Nearest-promoter candidate metrics do not match graph contexts"
        )
    distance_candidate_metrics = []
    for path in distance_candidate_metric_paths:
        with path.open(encoding="utf-8") as handle:
            distance_candidate_metrics.append(json.load(handle))
    distance_candidates_by_context = {
        row["context"]: row for row in distance_candidate_metrics
    }
    if len(distance_candidates_by_context) != len(distance_candidate_metrics):
        raise ValueError("Distance candidate metrics repeat a context")
    powerlaw_contexts = {
        row["context"]
        for row in contexts
        if row["contact_strategy"] == "powerlaw"
    }
    if set(distance_candidates_by_context) != powerlaw_contexts:
        raise ValueError(
            "Distance candidate metrics do not match powerlaw graph contexts"
        )
    metrics = {
        "schema_version": 1,
        "context_count": len(contexts),
        "observed_context_count": sum(
            row["contact_strategy"] == "observed" for row in contexts
        ),
        "powerlaw_context_count": sum(
            row["contact_strategy"] == "powerlaw" for row in contexts
        ),
        "element_promoter_edge_count": sum(
            int(row["element_promoter_edge_count"]) for row in contexts
        ),
        "element_gene_candidate_count": sum(
            int(row["element_gene_candidate_count"]) for row in contexts
        ),
        "active_contact_enhancer_gene_candidate_count": sum(
            int(row["active_contact_enhancer_gene_candidate_count"])
            for row in candidate_metrics
        ),
        "nearest_active_promoter_gene_candidate_count": sum(
            int(row["nearest_active_promoter_gene_candidate_count"])
            for row in nearest_candidate_metrics
        ),
        "active_distance_enhancer_gene_candidate_count": sum(
            int(row["active_distance_enhancer_gene_candidate_count"])
            for row in distance_candidate_metrics
        ),
        "contexts": {
            row["context"]: {
                key: row[key]
                for key in (
                    "contact_strategy",
                    "contact_assay",
                    "contact_match",
                    "contact_resolution_bp",
                    "element_node_count",
                    "promoter_node_count",
                    "active_promoter_count",
                    "element_promoter_edge_count",
                    "element_gene_candidate_count",
                )
            }
            | {
                "active_contact_enhancer_gene_candidate_count": (
                    candidates_by_context[row["context"]][
                        "active_contact_enhancer_gene_candidate_count"
                    ]
                ),
                "nearest_active_promoter_gene_candidate_count": (
                    nearest_candidates_by_context[row["context"]][
                        "nearest_active_promoter_gene_candidate_count"
                    ]
                ),
                "active_distance_enhancer_gene_candidate_count": (
                    distance_candidates_by_context.get(row["context"], {}).get(
                        "active_distance_enhancer_gene_candidate_count", 0
                    )
                ),
            }
            for row in contexts
        },
    }
    write_json_atomic(output_metrics, metrics)
    inputs = [
        source_manifest,
        promoter_metrics,
        powerlaw,
        *contact_metrics,
        *context_metric_paths,
        *candidate_metric_paths,
        *nearest_candidate_metric_paths,
        *distance_candidate_metric_paths,
    ]
    provenance = {
        "schema_version": 1,
        "method": "context_contact_graph_v3",
        "inputs": {
            str(path): {"path": str(path), "sha256": sha256_file(path)}
            for path in inputs
        },
        "outputs": {
            "metrics": {
                "path": str(output_metrics),
                "sha256": sha256_file(output_metrics),
            }
        },
    }
    write_json_atomic(output_provenance, provenance)
    return metrics, provenance
