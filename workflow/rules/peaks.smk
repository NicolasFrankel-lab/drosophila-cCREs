rule callpeak_h3k27ac:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        validated=final_bam_validation_input
    output:
        peaks=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.broadPeak",
        gapped=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.gappedPeak",
        xls=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.xls",
        treat_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_treat_pileup.bdg",
        control_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        command=callpeak_argv
    wildcard_constraints:
        sample=CHIP_TREATMENT_RE
    conda:
        "../envs/peaks.yaml"
    resources:
        mem_mb=8000
    log:
        f"{RESULT_ROOT}/logs/peaks/{{sample}}.callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"
