import json


rule resolved_config_provenance:
    output:
        config=RESOLVED_CONFIG
    params:
        config=json.dumps(config, sort_keys=True)
    conda:
        "../envs/reporting.yaml"
    log:
        f"{RESULT_ROOT}/logs/provenance/resolved_config.log"
    script:
        "../scripts/write_provenance.py"


rule atac_shift_bam:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        bai=lambda wc: FINAL_BAIS[wc.sample],
        validated=final_bam_validation_input
    output:
        bam=temp(f"{WORK_ROOT}/atac_shift/{{sample}}.shifted.bam"),
        bai=temp(f"{WORK_ROOT}/atac_shift/{{sample}}.shifted.bam.bai")
    wildcard_constraints:
        sample=ATAC_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/tracks/{{sample}}.alignmentSieve.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        temporary=$(mktemp -d $(dirname {output.bam:q})/.{wildcards.sample}.shift.XXXXXX)
        trap 'rm -rf "$temporary"' EXIT
        alignmentSieve --ATACshift -b {input.bam:q} -o "$temporary/unsorted.bam" \
          --numberOfProcessors {threads} > {log:q} 2>&1
        samtools sort -@ {threads} -o "$temporary/shifted.bam" \
          "$temporary/unsorted.bam" 2>> {log:q}
        samtools index -@ {threads} "$temporary/shifted.bam" \
          "$temporary/shifted.bam.bai" 2>> {log:q}
        samtools quickcheck -v "$temporary/shifted.bam" 2>> {log:q}
        mv "$temporary/shifted.bam" {output.bam:q}
        mv "$temporary/shifted.bam.bai" {output.bai:q}
        """


rule atac_shifted_bigwig:
    input:
        bam=f"{WORK_ROOT}/atac_shift/{{sample}}.shifted.bam",
        bai=f"{WORK_ROOT}/atac_shift/{{sample}}.shifted.bam.bai"
    output:
        bw=temp(f"{RESULT_ROOT}/tracks/{{sample}}.Tn5-shifted.CPM.bw")
    wildcard_constraints:
        sample=ATAC_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/tracks/{{sample}}.shifted.bamCoverage.log"
    shell:
        "mkdir -p $(dirname {output.bw:q}) $(dirname {log:q}) && "
        "bamCoverage -b {input.bam:q} -o {output.bw:q} --outFileFormat bigwig "
        "--normalizeUsing CPM --binSize 10 --numberOfProcessors {threads} > {log:q} 2>&1"


rule atac_tss_matrix:
    input:
        bw=f"{RESULT_ROOT}/tracks/{{sample}}.Tn5-shifted.CPM.bw",
        tss=str(REFERENCE["tss_bed"])
    output:
        matrix=temp(f"{RESULT_ROOT}/qc/tss/{{sample}}.matrix.gz"),
        table=temp(f"{RESULT_ROOT}/qc/tss/{{sample}}.matrix.tsv")
    wildcard_constraints:
        sample=ATAC_RE
    threads: 4
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/{{sample}}.computeMatrix.log"
    shell:
        "mkdir -p $(dirname {output.matrix:q}) $(dirname {log:q}) && "
        "computeMatrix reference-point --referencePoint TSS -S {input.bw:q} -R {input.tss:q} "
        "-b 2000 -a 2000 --skipZeros --numberOfProcessors {threads} "
        "-o {output.matrix:q} --outFileNameMatrix {output.table:q} > {log:q} 2>&1"


rule atac_tss_profile:
    input:
        matrix=f"{RESULT_ROOT}/qc/tss/{{sample}}.matrix.gz"
    output:
        plot=f"{RESULT_ROOT}/qc/tss/{{sample}}.profile.png",
        table=f"{RESULT_ROOT}/qc/tss/{{sample}}.profile.tsv"
    wildcard_constraints:
        sample=ATAC_RE
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/{{sample}}.plotProfile.log"
    shell:
        "plotProfile -m {input.matrix:q} -out {output.plot:q} --perGroup "
        "--outFileNameData {output.table:q} > {log:q} 2>&1"


rule atac_fragment_histogram:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        validated=final_bam_validation_input
    output:
        histogram=f"{RESULT_ROOT}/qc/fragments/{{sample}}.histogram.tsv"
    wildcard_constraints:
        sample=ATAC_RE
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/{{sample}}.fragment_histogram.log"
    shell:
        r"""
        mkdir -p $(dirname {output.histogram:q})
        printf 'fragment_length\tcount\n' > {output.histogram:q}
        samtools view {input.bam:q} 2> {log:q} \
          | awk '$9 > 0 {{print $9}}' \
          | sort -n \
          | uniq -c \
          | awk 'BEGIN{{OFS="\t"}} {{print $2,$1}}' >> {output.histogram:q}
        """


rule frip:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        peaks=lambda wc: PEAKS[wc.sample],
        validated=final_bam_validation_input
    output:
        tsv=f"{RESULT_ROOT}/qc/frip/{{sample}}.tsv",
        json=f"{RESULT_ROOT}/qc/frip/{{sample}}.json"
    params:
        layout=lambda wc: SAMPLES[wc.sample]["layout"]
    threads: 4
    resources:
        mem_mb=4000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/{{sample}}.frip.log"
    script:
        "../scripts/frip.py"


rule chip_cross_correlation:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        validated=final_bam_validation_input
    output:
        metrics=f"{RESULT_ROOT}/qc/chip/{{sample}}.cross_correlation.txt",
        plot=f"{RESULT_ROOT}/qc/chip/{{sample}}.cross_correlation.pdf"
    wildcard_constraints:
        sample=CHIP_TREATMENT_RE
    threads: 4
    resources:
        mem_mb=6000
    conda:
        "../envs/chip_qc.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/{{sample}}.cross_correlation.log"
    shell:
        "mkdir -p $(dirname {output.metrics:q}) $(dirname {log:q}) && "
        "run_spp.R -c={input.bam:q} -p={threads} -savp={output.plot:q} "
        "-out={output.metrics:q} > {log:q} 2>&1"


METRICS_MANIFEST = {
    sample: {
        "assay": ASSAY,
        "layout": SAMPLES[sample]["layout"],
        "role": SAMPLES[sample]["role"],
        "flagstat": FLAGSTATS[sample],
        "frip": FRIP_JSON.get(sample),
        "peaks": PEAKS.get(sample),
    }
    for sample in SAMPLE_IDS
}


rule metrics_summary:
    input:
        flagstats=list(FLAGSTATS.values()),
        frip=list(FRIP_JSON.values()),
        peaks=list(PEAKS.values())
    output:
        tsv=METRICS_TSV,
        json=METRICS_JSON
    params:
        manifest=json.dumps(METRICS_MANIFEST, sort_keys=True)
    conda:
        "../envs/reporting.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/metrics_summary.log"
    script:
        "../scripts/aggregate_metrics.py"


rule multiqc:
    input:
        raw_fastqc=RAW_FASTQC,
        trimmed_fastqc=TRIMMED_FASTQC,
        cutadapt=CUTADAPT_JSON,
        alignment=ALIGNMENT_QC_FILES,
        frip=list(FRIP_TSV.values()),
        metrics=METRICS_TSV,
        cross_correlation=list(CROSSCORRELATION.values())
    output:
        report=MULTIQC_REPORT
    params:
        outdir=lambda wc, output: str(Path(output.report).parent),
        scan=lambda wc, output: str(Path(output.report).parents[2])
    conda:
        "../envs/reporting.yaml"
    log:
        f"{RESULT_ROOT}/logs/qc/multiqc.log"
    shell:
        "mkdir -p {params.outdir:q} $(dirname {log:q}) && "
        "multiqc {params.scan:q} --outdir {params.outdir:q} "
        "--filename multiqc_report.html --force > {log:q} 2>&1"


if SAMPLE_IDS:
    rule library_qc_review_table:
        input:
            manifest=FINAL_BAM_EXPORT_MANIFEST,
            metrics=METRICS_TSV,
            multiqc=MULTIQC_REPORT,
            tss_profiles=list(TSS_PROFILE_TABLES.values()),
            fragment_histograms=list(FRAGMENT_HISTOGRAMS.values()),
            cross_correlation=list(CROSSCORRELATION.values()),
            cross_correlation_plots=list(CROSSCORRELATION_PLOTS.values())
        output:
            table=LIBRARY_QC_REVIEW
        params:
            library_qc=json.dumps(LIBRARY_QC_REVIEW_INPUTS, sort_keys=True)
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/qc/library-review.log"
        script:
            "../scripts/build_library_qc_review.py"
