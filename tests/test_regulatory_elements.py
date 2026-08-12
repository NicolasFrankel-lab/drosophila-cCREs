import csv
import gzip
import math
from pathlib import Path

import pytest

from short_read_processing.activity import _tsv_content, write_deterministic_gzip
from short_read_processing.activity_tmm import (
    TMM_ACTIVITY_FIELDS,
    TMM_BACKGROUND_METHOD,
    TMM_FACTOR_FIELDS,
)
from short_read_processing.regulatory_elements import (
    CATALOG_FIELDS,
    H3K27AC_QUANTIFICATION_FIELDS,
    WINDOW_COUNT_FIELDS,
    WINDOW_FIELDS,
    _is_h3k27ac_active,
    _read_blacklist,
    blacklist_overlap,
    build_h3k27ac_quantification,
    build_regulatory_catalog,
    nearest_tss,
    regulatory_class,
    write_window_counts,
    write_window_definitions,
)


def test_h3k27ac_activity_uses_strict_060_cutoff_and_excludes_blacklist():
    assert not _is_h3k27ac_active(
        member=True, blacklisted=False, background_percentile=0.60
    )
    assert _is_h3k27ac_active(
        member=True, blacklisted=False, background_percentile=0.600001
    )
    assert not _is_h3k27ac_active(
        member=True, blacklisted=True, background_percentile=1.0
    )
    assert not _is_h3k27ac_active(
        member=False, blacklisted=False, background_percentile=1.0
    )


def _master_files(tmp_path: Path) -> tuple[Path, Path, Path]:
    master = tmp_path / "master.bed"
    summit = tmp_path / "summits.bed"
    chrom_sizes = tmp_path / "chrom.sizes"
    master.write_text(
        "chr1\t50\t151\tDHS0000001\t0\t.\n"
        "chr1\t400\t501\tDHS0000002\t0\t.\n"
    )
    summit.write_text(
        "chr1\t100\t101\tDHS0000001\t0\t.\n"
        "chr1\t450\t451\tDHS0000002\t0\t.\n"
    )
    chrom_sizes.write_text("chr1\t1000\n")
    return master, summit, chrom_sizes


def test_window_definitions_clip_boundaries_and_retain_empty_flanks(tmp_path):
    master, summit, chrom_sizes = _master_files(tmp_path)
    table = tmp_path / "windows.tsv.gz"
    bed = tmp_path / "windows.bed"

    metrics = write_window_definitions(
        master_bed=master,
        summit_bed=summit,
        chrom_sizes=chrom_sizes,
        output_table=table,
        output_bed=bed,
    )

    assert metrics == {
        "window_count": 6,
        "countable_window_count": 5,
        "zero_width_window_count": 1,
    }
    with gzip.open(table, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    first = {row["window"]: row for row in rows if row["master_dhs_id"] == "DHS0000001"}
    assert (int(first["left_500"]["start"]), int(first["left_500"]["end"])) == (0, 0)
    assert (int(first["center_500"]["start"]), int(first["center_500"]["end"])) == (0, 350)
    assert (int(first["right_500"]["start"]), int(first["right_500"]["end"])) == (350, 850)
    assert len(bed.read_text().splitlines()) == 5


def test_window_count_formatter_fills_zero_width_window(tmp_path):
    master, summit, chrom_sizes = _master_files(tmp_path)
    table = tmp_path / "windows.tsv.gz"
    bed = tmp_path / "windows.bed"
    write_window_definitions(
        master_bed=master,
        summit_bed=summit,
        chrom_sizes=chrom_sizes,
        output_table=table,
        output_bed=bed,
    )
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text(
        "".join(
            f"{line}\t{index}\n"
            for index, line in enumerate(bed.read_text().splitlines(), start=1)
        )
    )
    total = tmp_path / "total.txt"
    total.write_text("100\n")
    output = tmp_path / "counts.tsv.gz"

    write_window_counts(
        window_table=table,
        coverage_path=coverage,
        total_units_path=total,
        library_id="h3_1",
        context="ctx",
        output_path=output,
    )

    with gzip.open(output, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    empty = next(row for row in rows if row["window_id"] == "DHS0000001|left_500")
    assert empty["raw_count"] == "0"
    assert empty["width_bp"] == "0"


def test_nearest_tss_classes_use_inclusive_500_and_1000_boundaries():
    tss = {"chr1": ([100, 500], {100: ["a", "b"], 500: ["c"]})}
    assert nearest_tss("chr1", 350, tss) == (150, "c")
    assert nearest_tss("chrX", 350, tss) == (None, "")
    assert regulatory_class(500) == "promoter_associated"
    assert regulatory_class(501) == "proximal_enhancer_like"
    assert regulatory_class(1000) == "proximal_enhancer_like"
    assert regulatory_class(1001) == "distal_enhancer_like"
    assert regulatory_class(None) == "unclassified_no_tss_on_contig"


def test_blacklist_overlap_uses_unioned_interval_bases(tmp_path):
    path = tmp_path / "blacklist.bed"
    path.write_text(
        "chr1\t60\t90\n"
        "chr1\t80\t120\n"
        "chr1\t140\t160\n"
    )
    blacklist = _read_blacklist(path)

    overlap_bp, overlap_fraction = blacklist_overlap("chr1", 50, 151, blacklist)

    assert overlap_bp == 71
    assert overlap_fraction == pytest.approx(71 / 101)
    assert blacklist_overlap("chr2", 50, 151, blacklist) == (0, 0.0)


def _write_factor_table(path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TMM_FACTOR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for library_id in ("h3_1", "h3_2"):
            writer.writerow(
                {
                    "library_id": library_id,
                    "assay": "h3k27ac",
                    "context": "ctx",
                    "total_units": 1000,
                    "feature_count": 10,
                    "tmm_normalization_factor": 1,
                    "effective_library_size": 1000,
                    "normalization_method": TMM_BACKGROUND_METHOD,
                }
            )


def _write_activity_table(path: Path) -> None:
    rows = []
    for master_id, start, end, summit in (
        ("DHS0000001", 50, 151, 100),
        ("DHS0000002", 400, 501, 450),
    ):
        row = {field: 0 for field in TMM_ACTIVITY_FIELDS}
        row.update(
            {
                "master_dhs_id": master_id,
                "chrom": "chr1",
                "start": start,
                "end": end,
                "summit": summit,
                "width_bp": end - start,
                "context": "ctx",
                "normalization_method": TMM_BACKGROUND_METHOD,
                "atac_normalized_cpm_per_kb": 9.0,
            }
        )
        rows.append(row)
    write_deterministic_gzip(path, _tsv_content(TMM_ACTIVITY_FIELDS, rows))


def _write_window_count_table(
    path: Path,
    definitions: list[dict[str, str]],
    *,
    library_id: str,
    counts: dict[tuple[str, str], int],
) -> None:
    rows = []
    for definition in definitions:
        rows.append(
            {
                **definition,
                "library_id": library_id,
                "context": "ctx",
                "raw_count": counts[(definition["master_dhs_id"], definition["window"])],
                "total_units": 1000,
            }
        )
    write_deterministic_gzip(path, _tsv_content(WINDOW_COUNT_FIELDS, rows))


def test_catalog_averages_replicates_before_max_and_uses_center_tie_break(tmp_path):
    master, summit, chrom_sizes = _master_files(tmp_path)
    windows = tmp_path / "windows.tsv.gz"
    write_window_definitions(
        master_bed=master,
        summit_bed=summit,
        chrom_sizes=chrom_sizes,
        output_table=windows,
        output_bed=tmp_path / "windows.bed",
    )
    with gzip.open(windows, "rt", newline="") as handle:
        definitions = list(csv.DictReader(handle, delimiter="\t"))
    count_paths = {"h3_1": tmp_path / "h3_1.tsv.gz", "h3_2": tmp_path / "h3_2.tsv.gz"}
    first_counts = {}
    second_counts = {}
    for definition in definitions:
        key = (definition["master_dhs_id"], definition["window"])
        first_counts[key] = 0
        second_counts[key] = 0
    # On DHS2, center and left context means tie at 50 counts. Center must win.
    # DHS2's left window is boundary-clipped to 200 bp; 40 counts there have
    # the same per-kb context mean as 100 counts in the 500-bp center window.
    first_counts[("DHS0000002", "left_500")] = 40
    second_counts[("DHS0000002", "center_500")] = 100
    for library_id, values in (("h3_1", first_counts), ("h3_2", second_counts)):
        _write_window_count_table(
            count_paths[library_id], definitions, library_id=library_id, counts=values
        )
    factors = tmp_path / "factors.tsv"
    _write_factor_table(factors)
    activity = tmp_path / "activity.tsv.gz"
    _write_activity_table(activity)
    matrix = tmp_path / "matrix.tsv"
    matrix.write_text(
        "master_dhs_id\tchrom\tstart\tend\tsummit\tcontext_n\tctx\n"
        "DHS0000001\tchr1\t50\t151\t100\t0\t0\n"
        "DHS0000002\tchr1\t400\t501\t450\t1\t1\n"
    )
    tss = tmp_path / "tss.bed"
    tss.write_text("chr1\t100\t101\ttss1\t0\t+\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t425\t475\n")
    catalog = tmp_path / "catalog.tsv.gz"
    wide = tmp_path / "wide.tsv.gz"
    context_elements = tmp_path / "ctx.elements.tsv.gz"
    h3k27ac_signal = tmp_path / "h3k27ac_signal.tsv.gz"

    build_h3k27ac_quantification(
        master_bed=master,
        summit_bed=summit,
        context_matrix=matrix,
        blacklist_bed=blacklist,
        window_table=windows,
        window_count_paths=count_paths,
        factor_table=factors,
        contexts=["ctx"],
        output_signal=h3k27ac_signal,
        output_metrics=tmp_path / "h3k27ac_metrics.json",
        output_provenance=tmp_path / "h3k27ac_provenance.json",
    )

    with gzip.open(h3k27ac_signal, "rt", newline="") as handle:
        signal_reader = csv.DictReader(handle, delimiter="\t")
        assert signal_reader.fieldnames == H3K27AC_QUANTIFICATION_FIELDS
        signal_rows = list(signal_reader)
    second_signal = next(
        row for row in signal_rows if row["master_dhs_id"] == "DHS0000002"
    )
    assert second_signal["h3k27ac_max_500_window"] == "center_500"
    assert second_signal["h3k27ac_max_500_raw_count_sum"] == "100"
    assert second_signal["h3k27ac_max_500_background_n"] == "1"
    assert float(
        second_signal["h3k27ac_max_500_background_percentile"]
    ) == pytest.approx(1.0)
    assert float(
        second_signal[
            "h3k27ac_max_500_background_upper_tail_probability"
        ]
    ) == pytest.approx(0.5)

    build_regulatory_catalog(
        master_bed=master,
        summit_bed=summit,
        context_matrix=matrix,
        tss_bed=tss,
        blacklist_bed=blacklist,
        h3k27ac_quantification=h3k27ac_signal,
        activity_table=activity,
        contexts=["ctx"],
        output_catalog=catalog,
        output_wide=wide,
        output_element_paths={"ctx": context_elements},
        output_summary=tmp_path / "summary.tsv",
        output_metrics=tmp_path / "metrics.json",
        output_provenance=tmp_path / "provenance.json",
    )

    with gzip.open(catalog, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == CATALOG_FIELDS
        rows = list(reader)
    second = next(row for row in rows if row["master_dhs_id"] == "DHS0000002")
    assert second["h3k27ac_max_500_window"] == "center_500"
    assert second["h3k27ac_active"] == "0"
    assert second["activity_state"] == "excluded_blacklist"
    assert second["blacklist_overlap"] == "1"
    assert second["blacklist_overlap_bp"] == "50"
    assert float(second["blacklist_overlap_fraction"]) == pytest.approx(50 / 101)
    assert float(second["h3k27ac_left_500_normalized_cpm_per_kb"]) == pytest.approx(100000)
    assert float(second["h3k27ac_center_500_normalized_cpm_per_kb"]) == pytest.approx(100000)
    assert second["h3k27ac_max_500_raw_count_sum"] == "100"
    assert float(
        second["h3k27ac_max_500_background_percentile"]
    ) == pytest.approx(1.0)
    assert float(second["combined_activity_max_500"]) == pytest.approx(
        math.sqrt(9 * 100000)
    )
    with gzip.open(wide, "rt", newline="") as handle:
        wide_reader = csv.DictReader(handle, delimiter="\t")
        assert "ctx__h3k27ac_max_500_background_percentile" in wide_reader.fieldnames
        assert "ctx__h3k27ac_active" in wide_reader.fieldnames
        assert len(list(wide_reader)) == 2
    with gzip.open(context_elements, "rt", newline="") as handle:
        element_reader = csv.DictReader(handle, delimiter="\t")
        assert element_reader.fieldnames == CATALOG_FIELDS
        assert len(list(element_reader)) == 1
