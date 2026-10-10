import shlex

from short_read_processing.macs3 import (
    atac_qpois_callpeak_arguments,
    callpeak_arguments,
)
from short_read_processing.workflow_config import aria2_checksum


def raw_fastq(wildcards):
    return RAW_FASTQ_BY_UNIT[wildcards.unit]


def trimmed_fastq(wildcards):
    return TRIMMED_FASTQ_BY_UNIT[wildcards.unit]


def lane_input(wildcards, mate):
    return LANES[(wildcards.sample, wildcards.lane)][mate]


def trim_parameter(wildcards, section, name):
    return SAMPLES[wildcards.sample]["parameters"][section][name]


def adapter_value(wildcards, mate):
    trimming = SAMPLES[wildcards.sample]["parameters"]["trimming"]
    preset = trimming["adapter_preset"]
    if preset == "custom":
        return f"file:{trimming['adapter_fasta']}"
    if preset == "nextera":
        return NEXTERA_ADAPTER
    return TRUSEQ_R1_ADAPTER if mate == "r1" else TRUSEQ_R2_ADAPTER


def trimmed_lane_reads(wildcards):
    key = (wildcards.sample, wildcards.lane)
    reads = [TRIMMED_R1[key]]
    if SAMPLES[wildcards.sample]["layout"] == "paired":
        reads.append(TRIMMED_R2[key])
    return reads


def bowtie_lane_arguments(wildcards):
    key = (wildcards.sample, wildcards.lane)
    r1 = shlex.quote(TRIMMED_R1[key])
    if SAMPLES[wildcards.sample]["layout"] == "paired":
        r2 = shlex.quote(TRIMMED_R2[key])
        return f"-1 {r1} -2 {r2}"
    return f"-U {r1}"


def sample_lane_bams(wildcards):
    return [
        LANE_BAMS[(sample, lane)]
        for sample, lane in LANES
        if sample == wildcards.sample
    ]


def alignment_lane_from_unit(wildcards):
    return LANE_BY_UNIT[str(wildcards.unit)]


def chunked_lane_reads(wildcards):
    sample, lane = alignment_lane_from_unit(wildcards)
    reads = [TRIMMED_R1[(sample, lane)]]
    if SAMPLES[sample]["layout"] == "paired":
        reads.append(TRIMMED_R2[(sample, lane)])
    return reads


def chunked_lane_layout(wildcards):
    sample, _lane = alignment_lane_from_unit(wildcards)
    return SAMPLES[sample]["layout"]


def alignment_chunk_manifest(wildcards):
    checkpoint_result = checkpoints.split_lane_fastq.get(unit=str(wildcards.unit))
    chunk_root = Path(str(checkpoint_result.output.manifest)).parent
    manifest_path = chunk_root / "manifest.json"
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("unit") != str(wildcards.unit):
        raise ValueError(f"Alignment chunk manifest unit mismatch: {manifest_path}")
    if int(manifest.get("records_per_chunk", 0)) != ALIGNMENT_CHUNK_PAIRS:
        raise ValueError(f"Alignment chunk size mismatch: {manifest_path}")
    return chunk_root, manifest


def alignment_chunk_entry(wildcards):
    chunk_root, manifest = alignment_chunk_manifest(wildcards)
    chunk_id = str(wildcards.chunk)
    for entry in manifest["chunks"]:
        if entry["id"] == chunk_id:
            return chunk_root, entry
    raise ValueError(
        f"Alignment chunk {chunk_id!r} is absent from {chunk_root / 'manifest.json'}"
    )


def alignment_chunk_reads(wildcards):
    chunk_root, entry = alignment_chunk_entry(wildcards)
    reads = [str(chunk_root / entry["r1"])]
    if entry.get("r2"):
        reads.append(str(chunk_root / entry["r2"]))
    return reads


def alignment_chunk_bowtie_arguments(wildcards):
    reads = [shlex.quote(path) for path in alignment_chunk_reads(wildcards)]
    sample, _lane = alignment_lane_from_unit(wildcards)
    if SAMPLES[sample]["layout"] == "paired":
        return f"-1 {reads[0]} -2 {reads[1]}"
    return f"-U {reads[0]}"


def alignment_chunk_bams(wildcards):
    _chunk_root, manifest = alignment_chunk_manifest(wildcards)
    return [
        (
            f"{ALIGNMENT_CHUNK_ROOT}/bams/"
            f"{wildcards.unit}.{entry['id']}.coordsort.bam"
        )
        for entry in manifest["chunks"]
    ]


def alignment_chunk_logs(wildcards):
    _chunk_root, manifest = alignment_chunk_manifest(wildcards)
    return [
        (
            f"{ALIGNMENT_CHUNK_ROOT}/logs/"
            f"{wildcards.unit}.{entry['id']}.bowtie2.log"
        )
        for entry in manifest["chunks"]
    ]


def alignment_chunk_sample(wildcards):
    sample, _lane = alignment_lane_from_unit(wildcards)
    return sample


def alignment_chunk_layout_arguments(wildcards):
    sample, _lane = alignment_lane_from_unit(wildcards)
    if SAMPLES[sample]["layout"] == "paired":
        maximum = SAMPLES[sample]["parameters"]["alignment"][
            "maximum_fragment_length"
        ]
        return f"--no-mixed --no-discordant -X {maximum}"
    return ""


def alignment_chunk_preset(wildcards):
    sample, _lane = alignment_lane_from_unit(wildcards)
    return f"--{SAMPLES[sample]['parameters']['alignment']['preset']}"


def final_bam_validation_input(wildcards):
    validation = EXTERNAL_BAM_VALIDATIONS.get(wildcards.sample)
    return [validation] if validation else []


def worker_threads(wildcards, threads):
    return max(1, threads - 1)


def bowtie_layout_arguments(wildcards):
    sample = SAMPLES[wildcards.sample]
    if sample["layout"] == "paired":
        maximum = sample["parameters"]["alignment"]["maximum_fragment_length"]
        return f"--no-mixed --no-discordant -X {maximum}"
    return ""


def bowtie_preset(wildcards):
    preset = SAMPLES[wildcards.sample]["parameters"]["alignment"]["preset"]
    return f"--{preset}"


def required_flags(wildcards):
    return "-f 2" if SAMPLES[wildcards.sample]["layout"] == "paired" else ""


def excluded_flags(wildcards):
    excluded = 4 | 256 | 512 | 2048
    if SAMPLES[wildcards.sample]["layout"] == "paired":
        excluded |= 8
    if SAMPLES[wildcards.sample]["parameters"]["filtering"]["remove_duplicates"]:
        excluded |= 1024
    return excluded


def callpeak_argv(wildcards):
    sample = wildcards.sample
    arguments = callpeak_arguments(
        SAMPLES[sample]["peak_caller"],
        treatment_bam=Path(FINAL_BAMS[sample]),
        name=sample,
        genome_size=REFERENCE["macs3_genome_size"],
        output_dir=Path(f"{RESULT_ROOT}/peaks/{sample}"),
    )
    return shlex.join(arguments)


def atac_insertion_bam(wildcards):
    if wildcards.sample in ATAC_SHORT_BAMS:
        return ATAC_SHORT_BAMS[wildcards.sample]
    return FINAL_BAMS[wildcards.sample]


def atac_insertion_bai(wildcards):
    if wildcards.sample in ATAC_SHORT_BAIS:
        return ATAC_SHORT_BAIS[wildcards.sample]
    return FINAL_BAIS[wildcards.sample]


def atac_qpois_replicate_callpeak_argv(wildcards):
    sample = wildcards.sample
    arguments = atac_qpois_callpeak_arguments(
        SAMPLES[sample]["peak_caller"],
        insertion_bed=Path(ATAC_INSERTIONS[sample]).resolve(),
        name=sample,
        genome_size=REFERENCE["macs3_genome_size"],
        output_dir=Path("."),
    )
    return shlex.join(arguments)


def atac_condition_peak_config(condition):
    sample = ATAC_CONDITIONS[condition].samples[0]
    return SAMPLES[sample]["peak_caller"]


def atac_qpois_condition_callpeak_argv(wildcards):
    condition = wildcards.condition
    arguments = atac_qpois_callpeak_arguments(
        atac_condition_peak_config(condition),
        insertion_bed=Path(ATAC_CONDITION_INSERTIONS[condition]).resolve(),
        name=condition,
        genome_size=REFERENCE["macs3_genome_size"],
        output_dir=Path("."),
    )
    return shlex.join(arguments)


def atac_condition_insertion_inputs(wildcards):
    return [
        ATAC_INSERTIONS[sample]
        for sample in ATAC_CONDITIONS[wildcards.condition].samples
    ]


def atac_condition_insertion_count_inputs(wildcards):
    return [
        ATAC_INSERTION_COUNTS[sample]
        for sample in ATAC_CONDITIONS[wildcards.condition].samples
    ]


def atac_condition_peak_validation_inputs(wildcards):
    return [
        EXTERNAL_QC_PEAK_VALIDATIONS[sample]
        for sample in ATAC_CONDITIONS[wildcards.condition].samples
        if sample in EXTERNAL_QC_PEAK_VALIDATIONS
    ]


def atac_consensus_replicate_arguments(wildcards):
    condition = wildcards.condition
    return " ".join(
        "--replicate {}".format(
            shlex.quote(f"{sample}={ATAC_REPLICATE_REFINED[sample]}")
        )
        for sample in ATAC_CONDITIONS[condition].samples
    )


def atac_master_arguments(_wildcards):
    arguments = []
    for condition in ATAC_CONDITION_IDS:
        arguments.extend(
            ["--context-peaks", f"{condition}={ATAC_CONDITION_CONSENSUS[condition]}"]
        )
        arguments.extend(
            ["--context-signal", f"{condition}={ATAC_CONDITION_SIGNAL_BW[condition]}"]
        )
    return shlex.join(arguments)
