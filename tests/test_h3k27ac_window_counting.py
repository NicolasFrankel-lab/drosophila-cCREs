import csv
import gzip
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
from snakemake.utils import format as format_shell

from short_read_processing.regulatory_elements import write_window_definitions


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(
    not shutil.which("bedtools") or not shutil.which("pigz"),
    reason="bedtools and pigz are required",
)
@pytest.mark.parametrize("truncated", [False, True])
def test_window_counting_with_fragments_after_last_query_contig(tmp_path, truncated):
    master = tmp_path / "master.bed"
    summits = tmp_path / "summits.bed"
    sizes = tmp_path / "chrom.sizes"
    master.write_text("chr1\t400\t600\tDHS0000001\t0\t.\n")
    summits.write_text("chr1\t500\t501\tDHS0000001\t0\t.\n")
    sizes.write_text("chr1\t2000\nchr2\t200000\n")
    windows = tmp_path / "windows.tsv.gz"
    windows_bed = tmp_path / "windows.bed"
    write_window_definitions(
        master_bed=master,
        summit_bed=summits,
        chrom_sizes=sizes,
        output_table=windows,
        output_bed=windows_bed,
    )
    units = tmp_path / "units.bed.gz"
    with gzip.open(units, "wt") as handle:
        handle.write("chr1\t450\t451\tr0\t0\t.\n")
        # Exceed the pipe buffer after the last queried chromosome. The former
        # pigz | bedtools command produced correct counts but failed with 141.
        for index in range(100000):
            handle.write(f"chr2\t{index}\t{index + 1}\tr{index + 1}\t0\t.\n")
    if truncated:
        units.write_bytes(units.read_bytes()[:-8])
    total = tmp_path / "total.txt"
    total.write_text("100001\n")
    output = tmp_path / "counts.tsv.gz"
    log = tmp_path / "count.log"

    # Execute the production rule's shell, including decompression validation
    # and the real count formatter, without unrelated upstream workflow jobs.
    rules = (REPO_ROOT / "workflow/rules/activity.smk").read_text()
    block = rules.split("    rule count_activity_h3k27ac_windows:", 1)[1]
    shell = block.split('r"""', 1)[1].split('"""', 1)[0]
    command = format_shell(
        shell,
        input=SimpleNamespace(
            units=units,
            total=total,
            windows=windows,
            windows_bed=windows_bed,
            chrom_sizes=sizes,
            script=REPO_ROOT / "src/format_regulatory_window_counts.py",
        ),
        output=SimpleNamespace(counts=output),
        params=SimpleNamespace(context="ctx"),
        wildcards=SimpleNamespace(library="h3_1"),
        log=log,
    )
    result = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", command],
        capture_output=True,
        text=True,
        timeout=30,
    )
    if truncated:
        assert result.returncode != 0
        assert not output.exists()
        return
    assert result.returncode == 0, result.stderr + log.read_text()
    with gzip.open(output, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["window"]: int(row["raw_count"]) for row in rows} == {
        "left_500": 0,
        "center_500": 1,
        "right_500": 0,
    }
