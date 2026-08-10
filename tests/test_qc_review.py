import csv
from pathlib import Path

import pytest

from short_read_processing.artifacts import FINAL_BAM_FILTERING_CONTRACT
from short_read_processing.qc_review import (
    build_review_table,
    parse_cross_correlation,
    parse_fragment_histogram,
    parse_tss_profile,
)


def test_qc_metric_parsers(tmp_path):
    profile = tmp_path / "profile.tsv"
    values = [1.0] * 400
    values[199] = 4.0
    profile.write_text(
        "bin labels\n"
        "bins\t\t" + "\t".join(str(index) for index in range(1, 401)) + "\n"
        "sample\tgenes\t" + "\t".join(str(value) for value in values) + "\n"
    )
    assert parse_tss_profile(profile) == pytest.approx(4.0)

    histogram = tmp_path / "histogram.tsv"
    histogram.write_text("fragment_length\tcount\n100\t3\n200\t1\n")
    assert parse_fragment_histogram(histogram) == {
        "atac_fragments_counted": "4",
        "atac_median_fragment_length_bp": "100",
        "atac_fraction_lt150bp": "0.75",
        "atac_fraction_180_250bp": "0.25",
    }

    cross_correlation = tmp_path / "cross_correlation.txt"
    cross_correlation.write_text(
        "h3.final.bam\t1000\t165,300\t0.8,0.7\t55\t0.6\t1500\t0.5\t"
        "1.2\t1.5\t1\n"
    )
    parsed = parse_cross_correlation(cross_correlation)
    assert parsed["suggested_fragment_length_bp"] == "165"
    assert parsed["phantompeak_nsc"] == "1.2"
    assert parsed["phantompeak_rsc"] == "1.5"


def test_build_review_table_combines_assay_specific_metrics(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    manifest = tmp_path / "final-bams.tsv"
    manifest_rows = []
    for library_id, assay, layout in (
        ("atac", "atac", "paired"),
        ("h3", "h3k27ac", "single"),
    ):
        bam = input_dir / f"{library_id}.bam"
        bai = input_dir / f"{library_id}.bam.bai"
        bam.write_bytes(b"bam")
        bai.write_bytes(b"bai")
        manifest_rows.append(
            f"{library_id}\t{assay}\teye\ttreatment\t{layout}\t{bam}\t{bai}\t"
            f"dm6\t{FINAL_BAM_FILTERING_CONTRACT}\t{'a' * 64}\t{'b' * 64}\t"
            "pending_review\n"
        )
    manifest.write_text(
        "library_id\tassay\tcontext\trole\tlayout\tbam\tbai\tgenome\t"
        "filtering_contract\tbam_sha256\tbai_sha256\tqc_status\n"
        + "".join(manifest_rows)
    )

    metrics = tmp_path / "metrics.tsv"
    metrics.write_text(
        "sample\tassay\tlayout\trole\ttotal_alignments\tmapped_alignments\t"
        "properly_paired\tfrip_numerator\tfrip_denominator\tfrip\tpeak_count\n"
        "atac\tatac\tpaired\ttreatment\t200\t190\t180\t25\t90\t0.277778\t10\n"
        "h3\tchip_histone\tsingle\ttreatment\t100\t95\t\t10\t95\t0.105263\t5\n"
    )
    profile = tmp_path / "atac.profile.tsv"
    values = [1.0] * 400
    values[199] = 4.0
    profile.write_text(
        "bin labels\n"
        "bins\t\t" + "\t".join(str(index) for index in range(1, 401)) + "\n"
        "sample\tgenes\t" + "\t".join(str(value) for value in values) + "\n"
    )
    histogram = tmp_path / "atac.histogram.tsv"
    histogram.write_text("fragment_length\tcount\n100\t3\n200\t1\n")
    cross_correlation = tmp_path / "h3.cross_correlation.txt"
    cross_correlation.write_text(
        "h3.final.bam\t1000\t165\t0.8\t55\t0.6\t1500\t0.5\t1.2\t1.5\t1\n"
    )
    cross_plot = tmp_path / "h3.cross_correlation.pdf"
    cross_plot.write_bytes(b"pdf")
    multiqc = tmp_path / "multiqc_report.html"
    multiqc.write_text("report")
    output = tmp_path / "qc" / "library-review.tsv"

    assert build_review_table(
        final_bam_manifest=manifest,
        metrics_tsv=metrics,
        library_qc={
            "atac": {
                "tss_profile": str(profile),
                "fragment_histogram": str(histogram),
            },
            "h3": {
                "cross_correlation": str(cross_correlation),
                "cross_correlation_plot": str(cross_plot),
            },
        },
        multiqc_report=multiqc,
        output=output,
    ) == 2

    with output.open(newline="") as handle:
        rows = {
            row["library_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert rows["atac"]["atac_tss_enrichment"] == "4"
    assert rows["atac"]["properly_paired_fraction"] == "0.9"
    assert rows["atac"]["qc_decision"] == ""
    assert rows["h3"]["phantompeak_nsc"] == "1.2"
    assert rows["h3"]["estimated_fragment_length_bp"] == "165"
