import json
from pathlib import Path

import pytest

from short_read_processing.activity import _tsv_content, sha256_file
from short_read_processing.activity_tmm import TMM_FACTOR_FIELDS
from short_read_processing.integrated_report import (
    _qc_rows,
    _qc_source_roots,
    _selected_qc_source,
    build_integrated_qc_report,
    discover_report_source_files,
)
from short_read_processing.regulatory_elements import SUMMARY_FIELDS


def test_report_qc_matches_selected_bam_when_library_ids_are_reused(tmp_path):
    records = []
    for name, digest, mapped in [("old", "a" * 64, 100), ("new", "b" * 64, 900)]:
        root = tmp_path / name
        (root / "qc").mkdir(parents=True)
        manifest = root / "provenance/manifests/final-bams.tsv"
        manifest.parent.mkdir(parents=True)
        manifest.write_text(f"library_id\tbam_sha256\ne13_h3\t{digest}\n")
        (root / "qc/metrics.json").write_text(json.dumps([
            {"sample": "e13_h3", "mapped_alignments": mapped},
            {"sample": "unselected", "mapped_alignments": 50},
        ]))
        records.extend(discover_report_source_files([root]))
    libraries = {"e13_h3": {"bam_sha256": "b" * 64}}
    roots = _qc_source_roots(records, libraries)
    assert _qc_rows(records, roots) == [{"sample": "e13_h3", "mapped_alignments": 900}]
    assert roots == {"e13_h3": {(tmp_path / "new").resolve()}}
    assert _qc_rows(list(reversed(records)), roots) == _qc_rows(records, roots)
    assert not _selected_qc_source("e13_h3", {"source_root": tmp_path / "old"}, roots)
    assert _selected_qc_source("e13_h3", {"source_root": tmp_path / "new"}, roots)
    missing = _qc_source_roots(records, {"e13_h3": {"bam_sha256": "c" * 64}})
    assert _qc_rows(records, missing) == []

    # A conflicting metric within the matched source remains an error.
    (tmp_path / "new/qc/metrics.json").write_text(json.dumps([
        {"sample": "e13_h3", "mapped_alignments": 900},
        {"sample": "e13_h3", "mapped_alignments": 901},
    ]))
    with pytest.raises(ValueError, match="Conflicting upstream QC metrics"):
        _qc_rows(records, roots)


def test_integrated_report_discovers_sources_and_writes_deterministic_outputs(tmp_path):
    upstream = tmp_path / "upstream"
    (upstream / "qc").mkdir(parents=True)
    qc = upstream / "qc" / "metrics.json"
    qc.write_text(
        json.dumps(
            [
                {
                    "sample": "atac_1",
                    "assay": "atac",
                    "layout": "paired",
                    "mapped_alignments": 1000,
                    "frip": 0.42,
                    "peak_count": 50,
                }
            ]
        )
    )
    discovered = discover_report_source_files([upstream])
    assert [record["kind"] for record in discovered] == ["qc_metrics"]

    master = tmp_path / "master.json"
    master.write_text(
        json.dumps(
            {
                "master_dhs_count": 100,
                "source_peak_count": 160,
                "multi_context_master_dhs_count": 40,
                "master_width_min": 50,
                "master_width_median": 200,
                "master_width_max": 400,
                "context_peak_counts": {"ctx": 80},
            }
        )
    )
    activity_metrics = tmp_path / "activity.json"
    activity_metrics.write_text(
        json.dumps(
            {
                "master_dhs_count": 100,
                "context_count": 1,
                "library_count": 4,
                "normalization_method": "tmm_background_10kb_v1",
            }
        )
    )
    catalog_metrics = tmp_path / "catalog.json"
    catalog_metrics.write_text(
        json.dumps(
            {
                "master_dhs_count": 100,
                "context_count": 1,
                "catalog_row_count": 100,
                "activity_definition": {"threshold": 0.6},
                "active_blacklist_overlap_element_count_by_context": {"ctx": 2},
            }
        )
    )
    factors = tmp_path / "factors.tsv"
    factor_rows = []
    for library_id, assay in (
        ("atac_1", "atac"),
        ("atac_2", "atac"),
        ("h3_1", "h3k27ac"),
        ("h3_2", "h3k27ac"),
    ):
        factor_rows.append(
            {
                "library_id": library_id,
                "assay": assay,
                "context": "ctx",
                "total_units": 1000,
                "feature_count": 20,
                "tmm_normalization_factor": 1.0,
                "effective_library_size": 1000,
                "normalization_method": "tmm_background_10kb_v1",
            }
        )
    factors.write_text(_tsv_content(TMM_FACTOR_FIELDS, factor_rows))

    summary = tmp_path / "summary.tsv"
    summary.write_text(
        _tsv_content(
            SUMMARY_FIELDS,
            [
                {
                    "context": "ctx",
                    "regulatory_class": "distal_enhancer_like",
                    "h3k27ac_active": 1,
                    "activity_state": "active_h3k27ac_background_gt_0.60",
                    "element_n": 10,
                }
            ],
        )
    )
    nearest_candidates = tmp_path / "nearest.tsv.gz"
    nearest_candidates.write_bytes(b"candidates")
    nearest_metrics = tmp_path / "nearest.json"
    nearest_metrics.write_text(
        json.dumps(
            {
                "enhancer_count": 10,
                "nearest_tss_tie_element_count": 1,
                "no_same_chromosome_promoter_count": 0,
                "enhancer_with_active_nearest_promoter_by_context": {"ctx": 8},
            }
        )
    )

    config = {
        "project": "atlas",
        "run_id": "catalog-v1",
        "assay": "activity",
        "reference": {"name": "dm6"},
        "activity": {
            "contexts": ["ctx"],
            "master": {"stats_json": str(master)},
            "libraries": [
                {
                    "id": library_id,
                    "context": "ctx",
                    "assay": assay,
                    "layout": "paired",
                    "qc_status": "accepted",
                }
                for library_id, assay in (
                    ("atac_1", "atac"),
                    ("atac_2", "atac"),
                    ("h3_1", "h3k27ac"),
                    ("h3_2", "h3k27ac"),
                )
            ],
        },
        "provenance": {},
    }
    current = {
        "activity_metrics": activity_metrics,
        "regulatory_element_metrics": catalog_metrics,
        "normalization_factors": factors,
        "regulatory_element_summary": summary,
        "nearest_tss_enhancer_candidates": nearest_candidates,
        "nearest_tss_enhancer_metrics": nearest_metrics,
    }
    output_html = tmp_path / "report.html"
    output_pdf = tmp_path / "report.pdf"
    output_metrics = tmp_path / "report.json"

    def fake_pdf_renderer(html_text: str, path: Path) -> str:
        assert "Integrated atlas regulatory-catalog QC report" in html_text
        path.write_bytes(b"%PDF-1.4\n% test\n")
        return "test-renderer"

    result = build_integrated_qc_report(
        config=config,
        source_files=[
            {
                **discovered[0],
                "path": str(discovered[0]["path"]),
                "source_root": str(discovered[0]["source_root"]),
            }
        ],
        current_files=current,
        output_html=output_html,
        output_pdf=output_pdf,
        output_metrics=output_metrics,
        pdf_renderer=fake_pdf_renderer,
    )
    first_hashes = tuple(
        sha256_file(path) for path in (output_html, output_pdf, output_metrics)
    )
    build_integrated_qc_report(
        config=config,
        source_files=[
            {
                **discovered[0],
                "path": str(discovered[0]["path"]),
                "source_root": str(discovered[0]["source_root"]),
            }
        ],
        current_files=current,
        output_html=output_html,
        output_pdf=output_pdf,
        output_metrics=output_metrics,
        pdf_renderer=fake_pdf_renderer,
    )

    assert "contact_links" not in result
    assert "Nearest-TSS enhancer annotation" in output_html.read_text()
    assert "contact" not in output_html.read_text().lower()
    assert "&gt; 0.60 background percentile" in output_html.read_text()
    assert "Blacklist-flagged active calls" in output_html.read_text()
    assert first_hashes == tuple(
        sha256_file(path) for path in (output_html, output_pdf, output_metrics)
    )
