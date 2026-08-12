import csv
import gzip
import json
from pathlib import Path

import pytest

from short_read_processing.activity_tmm import TMM_BACKGROUND_METHOD, TMM_FACTOR_FIELDS
from short_read_processing.catalog_visualization import (
    build_catalog_beds,
    centered_extension_flanks,
    mean_unionbedg_rows,
    read_track_factors,
)


def test_catalog_beds_preserve_membership_and_continuous_annotations(tmp_path):
    master = tmp_path / "master.bed"
    summits = tmp_path / "summits.bed"
    matrix = tmp_path / "matrix.tsv"
    master.write_text(
        "chr1\t10\t40\tDHS0000001\t0\t.\n"
        "chr1\t50\t90\tDHS0000002\t0\t.\n"
    )
    summits.write_text(
        "chr1\t20\t21\tDHS0000001\t0\t.\n"
        "chr1\t70\t71\tDHS0000002\t0\t.\n"
    )
    matrix.write_text(
        "master_dhs_id\tchrom\tstart\tend\tsummit\tcontext_n\tctx\n"
        "DHS0000001\tchr1\t10\t40\t20\t1\t1\n"
        "DHS0000002\tchr1\t50\t90\t70\t0\t0\n"
    )
    elements = tmp_path / "ctx.elements.tsv.gz"
    fields = [
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
    ]
    with gzip.open(elements, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "master_dhs_id": "DHS0000001",
                "chrom": "chr1",
                "start": 10,
                "end": 40,
                "summit": 20,
                "context": "ctx",
                "context_membership": 1,
                "regulatory_class": "distal_enhancer_like",
                "nearest_tss_distance_bp": 1250,
                "blacklist_overlap": 1,
                "h3k27ac_max_500_background_percentile": 0.25,
                "h3k27ac_active": 0,
                "activity_state": "excluded_blacklist",
            }
        )
    output_master = tmp_path / "bed/master_dhs.bed"
    output_dhs = tmp_path / "bed/ctx.dhs.bed"
    output_elements = tmp_path / "bed/ctx.elements.bed"
    output_manifest = tmp_path / "bed/bed_tracks.json"

    result = build_catalog_beds(
        master_bed=master,
        summit_bed=summits,
        context_matrix=matrix,
        element_paths={"ctx": elements},
        output_master_bed=output_master,
        output_context_dhs={"ctx": output_dhs},
        output_element_beds={"ctx": output_elements},
        output_manifest=output_manifest,
    )

    assert output_master.read_text() == master.read_text()
    assert output_dhs.read_text().splitlines() == [
        "chr1\t10\t40\tDHS0000001\t0\t.\t20\t21\t0,145,130"
    ]
    element_fields = output_elements.read_text().strip().split("\t")
    assert element_fields[3] == (
        "DHS0000001|distal_enhancer_like|tss_distance=1250|"
        "max3_background_percentile=0.25|inactive|excluded_blacklist|blacklist_overlap"
    )
    assert element_fields[4] == "250"
    assert element_fields[8] == "0,0,0"
    assert result["context_metrics"]["ctx"]["context_element_count"] == 1
    assert json.loads(output_manifest.read_text())["status"] == "ok"


def test_track_factors_and_unionbedg_mean_are_exact(tmp_path):
    assert centered_extension_flanks(150) == (75, 74)
    assert sum(centered_extension_flanks(150)) + 1 == 150
    with pytest.raises(ValueError, match="positive"):
        centered_extension_flanks(0)

    factors = tmp_path / "factors.tsv"
    with factors.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TMM_FACTOR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for library_id, effective_size in (("rep1", 100.0), ("rep2", 200.0)):
            writer.writerow(
                {
                    "library_id": library_id,
                    "assay": "atac",
                    "context": "ctx",
                    "total_units": 100,
                    "feature_count": 10,
                    "tmm_normalization_factor": effective_size / 100,
                    "effective_library_size": effective_size,
                    "normalization_method": TMM_BACKGROUND_METHOD,
                }
            )

    parsed = read_track_factors(
        factor_path=factors,
        library_ids=["rep2", "rep1"],
        assay="atac",
        context="ctx",
    )
    assert list(parsed) == ["rep2", "rep1"]
    assert parsed["rep1"]["scale"] == 10_000
    assert parsed["rep2"]["scale"] == 5_000

    rows = list(
        mean_unionbedg_rows(
            ["chr1\t0\t10\t2\t0\n", "chr1\t10\t20\t2\t4\n"],
            chromosome_sizes=[("chr1", 20)],
            library_n=2,
        )
    )
    assert rows == [("chr1", 0, 10, 1.0), ("chr1", 10, 20, 3.0)]

    single_library_rows = list(
        mean_unionbedg_rows(
            ["chr1\t0\t10\t2.5\n"],
            chromosome_sizes=[("chr1", 20)],
            library_n=1,
        )
    )
    assert single_library_rows == [("chr1", 0, 10, 2.5)]

    with pytest.raises(ValueError, match="unsorted or overlapping"):
        list(
            mean_unionbedg_rows(
                ["chr1\t10\t20\t1\n", "chr1\t5\t8\t1\n"],
                chromosome_sizes=[("chr1", 20)],
                library_n=1,
            )
        )
