"""TOBIAS ATAC footprinting per context over the master DHS registry."""


if FOOTPRINTS:
    rule validate_footprint_bam:
        input:
            bam=lambda wc: FOOTPRINT_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: FOOTPRINT_LIBRARIES[wc.library]["bai"],
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=f"{RESULT_ROOT}/provenance/footprints/bams/{{library}}.validated.json"
        params:
            expected=lambda wc: FOOTPRINT_LIBRARIES[wc.library]
        wildcard_constraints:
            library=FOOTPRINT_LIBRARY_RE
        threads: 2
        resources:
            mem_mb=2000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/footprints/validation/{{library}}.log"
        script:
            "../scripts/validate_external_bam.py"


    # Replicates of a context are pooled before bias correction for more cut sites.
    rule merge_footprint_context_bams:
        input:
            bams=lambda wc: [FOOTPRINT_LIBRARIES[library]["bam"] for library in FOOTPRINT_CONTEXT_LIBRARIES[wc.context]],
            validations=lambda wc: [FOOTPRINT_BAM_VALIDATIONS[library] for library in FOOTPRINT_CONTEXT_LIBRARIES[wc.context]]
        output:
            bam=temp(f"{FOOTPRINT_WORK}/{{context}}.merged.bam"),
            bai=temp(f"{FOOTPRINT_WORK}/{{context}}.merged.bam.bai")
        wildcard_constraints:
            context=FOOTPRINT_CONTEXT_RE
        threads: 8
        resources:
            mem_mb=4000
        conda:
            "../envs/footprinting.yaml"
        log:
            f"{RESULT_ROOT}/logs/footprints/{{context}}.merge.log"
        shell:
            r"""
            mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.bam:q})/.{wildcards.context}.merge.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            samtools merge -f -@ {threads} -o "$temporary/merged.bam" {input.bams:q} 2> {log:q}
            samtools index -@ {threads} "$temporary/merged.bam" "$temporary/merged.bam.bai" 2>> {log:q}
            samtools quickcheck -v "$temporary/merged.bam" 2>> {log:q}
            mv "$temporary/merged.bam" {output.bam:q}
            mv "$temporary/merged.bam.bai" {output.bai:q}
            """


    rule tobias_atacorrect:
        input:
            bam=f"{FOOTPRINT_WORK}/{{context}}.merged.bam",
            bai=f"{FOOTPRINT_WORK}/{{context}}.merged.bam.bai",
            fasta=str(REFERENCE["fasta"]),
            fasta_index=FASTA_INDEX,
            peaks=ATAC_MASTER_BED,
            blacklist=str(REFERENCE["blacklist_bed"]),
            master_validation=EXTERNAL_MASTER_VALIDATION
        output:
            corrected=f"{FOOTPRINT_ROOT}/{{context}}.corrected.bw",
            uncorrected=f"{FOOTPRINT_ROOT}/{{context}}.uncorrected.bw",
            report=f"{FOOTPRINT_ROOT}/qc/{{context}}.atacorrect.pdf"
        params:
            read_shift=" ".join(str(value) for value in FOOTPRINTS["atacorrect"]["read_shift"]),
            k_flank=FOOTPRINTS["atacorrect"]["k_flank"],
            bg_shift=FOOTPRINTS["atacorrect"]["bg_shift"],
            window=FOOTPRINTS["atacorrect"]["window"],
            score_mat=FOOTPRINTS["atacorrect"]["score_mat"]
        wildcard_constraints:
            context=FOOTPRINT_CONTEXT_RE
        threads: 8
        resources:
            mem_mb=12000
        conda:
            "../envs/footprinting.yaml"
        log:
            f"{RESULT_ROOT}/logs/footprints/{{context}}.atacorrect.log"
        shell:
            r"""
            mkdir -p $(dirname {output.report:q}) $(dirname {log:q})
            temporary=$(mktemp -d {FOOTPRINT_ROOT}/.{wildcards.context}.atacorrect.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            TOBIAS ATACorrect --bam {input.bam:q} --genome {input.fasta:q} \
              --peaks {input.peaks:q} --blacklist {input.blacklist:q} \
              --read_shift {params.read_shift} --k_flank {params.k_flank} \
              --bg_shift {params.bg_shift} --window {params.window} \
              --score_mat {params.score_mat} --track-off bias expected \
              --prefix {wildcards.context} --outdir "$temporary" --cores {threads} \
              > {log:q} 2>&1
            mv "$temporary/{wildcards.context}_corrected.bw" {output.corrected:q}
            mv "$temporary/{wildcards.context}_uncorrected.bw" {output.uncorrected:q}
            mv "$temporary/{wildcards.context}_atacorrect.pdf" {output.report:q}
            """


    rule tobias_scorebigwig:
        input:
            signal=f"{FOOTPRINT_ROOT}/{{context}}.corrected.bw",
            regions=ATAC_MASTER_BED
        output:
            scores=f"{FOOTPRINT_ROOT}/{{context}}.footprints.bw"
        params:
            score=FOOTPRINTS["scorebigwig"]["score"],
            fp_min=FOOTPRINTS["scorebigwig"]["fp_min"],
            fp_max=FOOTPRINTS["scorebigwig"]["fp_max"],
            flank_min=FOOTPRINTS["scorebigwig"]["flank_min"],
            flank_max=FOOTPRINTS["scorebigwig"]["flank_max"],
            extend_bp=FOOTPRINTS["scorebigwig"]["extend"]
        wildcard_constraints:
            context=FOOTPRINT_CONTEXT_RE
        threads: 8
        resources:
            mem_mb=8000
        conda:
            "../envs/footprinting.yaml"
        log:
            f"{RESULT_ROOT}/logs/footprints/{{context}}.scorebigwig.log"
        shell:
            r"""
            mkdir -p $(dirname {log:q})
            temporary=$(mktemp -d {FOOTPRINT_ROOT}/.{wildcards.context}.scores.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            TOBIAS ScoreBigwig --signal {input.signal:q} --regions {input.regions:q} \
              --score {params.score} --fp-min {params.fp_min} --fp-max {params.fp_max} \
              --flank-min {params.flank_min} --flank-max {params.flank_max} \
              --extend {params.extend_bp} --output "$temporary/scores.bw" --cores {threads} \
              > {log:q} 2>&1
            mv "$temporary/scores.bw" {output.scores:q}
            """
