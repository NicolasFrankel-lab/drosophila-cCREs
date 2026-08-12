import csv
import gzip
from pathlib import Path

from short_read_processing.contact_links import PROMOTER_FIELDS, _atomic_tsv
from short_read_processing.nearest_tss_links import (
    build_nearest_tss_enhancer_candidates,
)


CATALOG_FIELDS = (
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
    "blacklist_overlap",
)


def _read_tsv(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _catalog_row(
    master_id,
    summit,
    context,
    *,
    member,
    regulatory_class,
    posterior,
):
    return {
        "master_dhs_id": master_id,
        "chrom": "chr2L",
        "start": summit - 50,
        "end": summit + 50,
        "summit": summit,
        "context": context,
        "context_membership": int(member),
        "regulatory_class": regulatory_class,
        "atac_normalized_cpm_per_kb": 4 if member else 0,
        "h3k27ac_max_500_background_percentile": posterior if member else "",
        "h3k27ac_active": int(member),
        "activity_state": (
            "active_h3k27ac_background_gt_0.60" if member else "inactive_in_context"
        ),
        "combined_activity_max_500": 3 if member else 0,
        "blacklist_overlap": 0,
    }


def test_wide_nearest_tss_table_keeps_one_enhancer_row_and_ties(tmp_path):
    promoters = tmp_path / "promoters.tsv.gz"
    _atomic_tsv(
        promoters,
        PROMOTER_FIELDS,
        [
            {
                "promoter_id": "DM6PROM00000001",
                "gene_id": "gene_a",
                "gene_name": "A",
                "chrom": "chr2L",
                "tss": 1000,
                "strand": "+",
                "start": 500,
                "end": 1500,
                "transcript_ids": "tx_a",
                "transcript_n": 1,
                "selection_rule": "all_distinct_gene_tss",
            },
            {
                "promoter_id": "DM6PROM00000002",
                "gene_id": "gene_b",
                "gene_name": "B",
                "chrom": "chr2L",
                "tss": 2200,
                "strand": "-",
                "start": 1700,
                "end": 2700,
                "transcript_ids": "tx_b1;tx_b2",
                "transcript_n": 2,
                "selection_rule": "all_distinct_gene_tss",
            },
        ],
    )
    catalog = tmp_path / "catalog.tsv.gz"
    rows = []
    for context in ("c1", "c2"):
        rows.extend(
            [
                _catalog_row(
                    "enh_tie",
                    1600,
                    context,
                    member=context == "c1",
                    regulatory_class="proximal_enhancer_like",
                    posterior=0.8,
                ),
                _catalog_row(
                    "enh_b",
                    3000,
                    context,
                    member=True,
                    regulatory_class="proximal_enhancer_like",
                    posterior=0.7,
                ),
                _catalog_row(
                    "support_a",
                    1000,
                    context,
                    member=context == "c1",
                    regulatory_class="promoter_associated",
                    posterior=0.9,
                ),
                _catalog_row(
                    "support_b",
                    2200,
                    context,
                    member=context == "c2",
                    regulatory_class="promoter_associated",
                    posterior=0.8,
                ),
            ]
        )
    _atomic_tsv(catalog, CATALOG_FIELDS, rows)

    output = tmp_path / "wide.tsv.gz"
    metrics = build_nearest_tss_enhancer_candidates(
        catalog_path=catalog,
        promoters_path=promoters,
        contexts=["c1", "c2"],
        enhancer_classes=["proximal_enhancer_like", "distal_enhancer_like"],
        output=output,
        metrics_output=tmp_path / "wide.json",
    )

    by_id = {row["master_dhs_id"]: row for row in _read_tsv(output)}
    assert set(by_id) == {"enh_tie", "enh_b"}
    tied = by_id["enh_tie"]
    assert tied["nearest_tss_distance_bp"] == "600"
    assert tied["nearest_tss_tie_count"] == "2"
    assert tied["nearest_gene_ids"] == "gene_a;gene_b"
    assert tied["evidence_type"] == "nearest_annotated_tss"
    assert tied["nearest_promoter_active_contexts"] == "c1;c2"
    assert tied["c1__active_nearest_promoter_ids"] == "DM6PROM00000001"
    assert tied["c2__active_nearest_promoter_ids"] == "DM6PROM00000002"
    assert tied["c1__active_nearest_promoter_supporting_element_ids"] == "support_a"
    assert tied["c2__active_nearest_promoter_supporting_element_ids"] == "support_b"
    assert tied["c2__element_context_membership"] == "0"
    assert by_id["enh_b"]["nearest_gene_ids"] == "gene_b"
    assert by_id["enh_b"]["c1__nearest_promoter_active"] == "0"
    assert by_id["enh_b"]["c2__nearest_promoter_active"] == "1"
    assert metrics["enhancer_count"] == 2
    assert metrics["nearest_tss_tie_element_count"] == 1


def test_wide_nearest_tss_output_is_byte_deterministic(tmp_path):
    promoters = tmp_path / "promoters.tsv.gz"
    _atomic_tsv(
        promoters,
        PROMOTER_FIELDS,
        [
            {
                "promoter_id": "DM6PROM00000001",
                "gene_id": "gene_a",
                "gene_name": "A",
                "chrom": "chr2L",
                "tss": 1000,
                "strand": "+",
                "start": 500,
                "end": 1500,
                "transcript_ids": "tx_a",
                "transcript_n": 1,
                "selection_rule": "all_distinct_gene_tss",
            }
        ],
    )
    catalog = tmp_path / "catalog.tsv.gz"
    _atomic_tsv(
        catalog,
        CATALOG_FIELDS,
        [
            _catalog_row(
                "enh",
                1700,
                "c1",
                member=True,
                regulatory_class="proximal_enhancer_like",
                posterior=0.8,
            )
        ],
    )
    outputs = []
    for index in (1, 2):
        output = tmp_path / f"wide-{index}.tsv.gz"
        build_nearest_tss_enhancer_candidates(
            catalog_path=catalog,
            promoters_path=promoters,
            contexts=["c1"],
            enhancer_classes=["proximal_enhancer_like", "distal_enhancer_like"],
            output=output,
            metrics_output=tmp_path / f"wide-{index}.json",
        )
        outputs.append(output.read_bytes())

    assert outputs[0] == outputs[1]
