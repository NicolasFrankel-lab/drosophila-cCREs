from pathlib import Path

from short_read_processing.macs3 import (
    atac_qpois_callpeak_arguments,
    callpeak_arguments,
    callpeak_bedgraph_paths,
)


def test_callpeak_argv_and_bedgraph_outputs():
    config = {
        "command": "callpeak",
        "format": "BAM",
        "qvalue": 0.01,
        "broad": False,
        "nomodel": True,
        "shift": -75,
        "extsize": 150,
        "write_bedgraph": True,
        "spmr": True,
    }
    arguments = callpeak_arguments(
        config,
        treatment_bam=Path("treatment.bam"),
        name="atac_rep1",
        genome_size="dm",
        output_dir=Path("peaks"),
    )

    assert arguments[:2] == ["macs3", "callpeak"]
    assert arguments[arguments.index("-q") + 1] == "0.01"
    assert "-B" in arguments
    assert "--SPMR" in arguments
    assert "-c" not in arguments
    assert arguments[arguments.index("--shift") + 1] == "-75"
    assert arguments[arguments.index("--extsize") + 1] == "150"
    assert callpeak_bedgraph_paths(Path("peaks"), "atac_rep1") == (
        Path("peaks/atac_rep1_treat_pileup.bdg"),
        Path("peaks/atac_rep1_control_lambda.bdg"),
    )


def test_h3k27ac_broad_callpeak_is_ip_only():
    config = {
        "command": "callpeak",
        "format": "BAMPE",
        "qvalue": 0.01,
        "broad": True,
        "broad_cutoff": 0.1,
        "nomodel": False,
        "shift": None,
        "extsize": None,
        "write_bedgraph": True,
        "spmr": True,
    }
    arguments = callpeak_arguments(
        config,
        treatment_bam=Path("h3k27ac.bam"),
        name="h3k27ac_rep1",
        genome_size="dm",
        output_dir=Path("peaks"),
    )

    assert "-c" not in arguments
    assert "--broad" in arguments
    assert arguments[arguments.index("--broad-cutoff") + 1] == "0.1"


def test_atac_qpois_uses_two_ended_bed_without_spmr():
    arguments = atac_qpois_callpeak_arguments(
        {
            "command": "callpeak",
            "mode": "tn5_qpois",
            "format": "BED",
            "qvalue": 0.1,
            "broad": False,
            "nomodel": True,
            "shift": -75,
            "extsize": 150,
            "write_bedgraph": True,
            "spmr": False,
        },
        insertion_bed=Path("insertions.bed.gz"),
        name="atac_rep1",
        genome_size="dm",
        output_dir=Path("peaks"),
    )

    assert arguments[arguments.index("-f") + 1] == "BED"
    assert arguments[arguments.index("-t") + 1] == "insertions.bed.gz"
    assert "--keep-dup" in arguments
    assert "-B" in arguments
    assert "--SPMR" not in arguments
