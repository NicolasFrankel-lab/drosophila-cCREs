"""Build contact-independent nearest-TSS enhancer candidate tables."""

from __future__ import annotations

from bisect import bisect_left
import csv
import gzip
from pathlib import Path
from typing import Any

from .artifacts import sha256_file
from .contact_links import (
    Promoter,
    _atomic_tsv,
    promoter_activities,
    read_context_elements,
    read_promoters,
)
from .contacts import write_json_atomic


BASE_FIELDS = (
    "master_dhs_id",
    "element_chrom",
    "element_start",
    "element_end",
    "element_summit",
    "regulatory_class",
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
    "atac_normalized_cpm_per_kb",
    "h3k27ac_max_500_background_percentile",
    "h3k27ac_active",
    "activity_state",
    "combined_activity_max_500",
)
CATALOG_REQUIRED_FIELDS = set(CATALOG_STORED_FIELDS)


def _open_catalog(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8", newline="")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8", newline="")
    )


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def _catalog_rows(
    path: Path, contexts: list[str]
) -> tuple[list[str], dict[str, dict[str, dict[str, str]]]]:
    by_element: dict[str, dict[str, dict[str, str]]] = {}
    element_order: list[str] = []
    with _open_catalog(path) as handle:
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
        static = tuple(first[field] for field in ("chrom", "start", "end", "summit", "regulatory_class"))
        if any(
            tuple(row[field] for field in ("chrom", "start", "end", "summit", "regulatory_class"))
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
            read_context_elements(catalog_path, context)[0],
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
        "method": "nearest_annotated_tss_all_ties_wide_v1",
        "catalog": str(catalog_path),
        "catalog_sha256": sha256_file(catalog_path),
        "promoters": str(promoters_path),
        "promoters_sha256": sha256_file(promoters_path),
        "contexts": contexts,
        "enhancer_classes": enhancer_classes,
        "promoter_activity_definition": (
            "any context-member master DHS summit within inclusive +/-500 bp of "
            "the exact annotated TSS with catalog H3K27ac activity"
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
