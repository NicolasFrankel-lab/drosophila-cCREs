import json
from pathlib import Path

from short_read_processing.activity import _tsv_content, sha256_file
from short_read_processing.activity_tmm import TMM_FACTOR_FIELDS
from short_read_processing.integrated_report import (
    build_integrated_qc_report,
    discover_report_source_files,
)
from short_read_processing.regulatory_elements import SUMMARY_FIELDS


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
            {"master_dhs_count": 100, "context_count": 1, "catalog_row_count": 100}
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
    contact_metrics = tmp_path / "contact-metrics.json"
    contact_metrics.write_text(
        json.dumps(
            {
                "observed_context_count": 1,
                "powerlaw_context_count": 0,
                "element_promoter_edge_count": 250,
                "element_gene_candidate_count": 200,
                "active_contact_enhancer_gene_candidate_count": 50,
                "nearest_active_promoter_gene_candidate_count": 40,
                "active_distance_enhancer_gene_candidate_count": 0,
                "contexts": {
                    "ctx": {
                        "contact_strategy": "observed",
                        "contact_assay": "Micro-C",
                        "contact_match": "test_match",
                        "contact_resolution_bp": 5000,
                        "active_promoter_count": 20,
                        "element_promoter_edge_count": 250,
                        "element_gene_candidate_count": 200,
                        "active_contact_enhancer_gene_candidate_count": 50,
                        "nearest_active_promoter_gene_candidate_count": 40,
                        "active_distance_enhancer_gene_candidate_count": 0,
                    }
                },
            }
        )
    )
    contact_provenance = tmp_path / "contact-provenance.json"
    contact_provenance.write_text(json.dumps({"method": "test"}))
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
        "contact_graph_metrics": contact_metrics,
        "contact_graph_provenance": contact_provenance,
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

    assert result["contact_links"]["observed_context_count"] == 1
    assert "Default nearest-TSS enhancer candidates" in output_html.read_text()
    assert "Optional context contact links and candidate genes" in output_html.read_text()
    assert "Nearest-active-promoter candidates" in output_html.read_text()
    assert "Focused distance candidates" in output_html.read_text()
    assert "&gt; 0.60 background percentile" in output_html.read_text()
    assert first_hashes == tuple(
        sha256_file(path) for path in (output_html, output_pdf, output_metrics)
    )

    nearest_only_files = {
        name: path
        for name, path in current.items()
        if not name.startswith("contact_graph_")
    }
    nearest_only_html = tmp_path / "nearest-only.html"
    nearest_only = build_integrated_qc_report(
        config=config,
        source_files=[
            {
                **discovered[0],
                "path": str(discovered[0]["path"]),
                "source_root": str(discovered[0]["source_root"]),
            }
        ],
        current_files=nearest_only_files,
        output_html=nearest_only_html,
        output_pdf=tmp_path / "nearest-only.pdf",
        output_metrics=tmp_path / "nearest-only.json",
        pdf_renderer=fake_pdf_renderer,
    )
    assert nearest_only["contact_links"] is None
    assert "Default nearest-TSS enhancer candidates" in nearest_only_html.read_text()
    assert "Optional context contact links" not in nearest_only_html.read_text()
