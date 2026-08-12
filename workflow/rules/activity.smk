"""Master-DHS background-TMM quantification and regulatory catalog."""


if ACTIVITY:
    rule validate_activity_bam:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=f"{RESULT_ROOT}/provenance/activity/bams/{{library}}.validated.json"
        params:
            expected=lambda wc: ACTIVITY_LIBRARIES[wc.library]
        wildcard_constraints:
            library=ACTIVITY_LIBRARY_RE
        threads: 2
        resources:
            mem_mb=2000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/validation/{{library}}.log"
        script:
            "../scripts/validate_external_bam.py"


    rule validate_activity_master:
        input:
            master_bed=str(ACTIVITY["master"]["master_bed"]),
            summits_bed=str(ACTIVITY["master"]["summits_bed"]),
            membership_tsv=str(ACTIVITY["master"]["membership_tsv"]),
            context_matrix_tsv=str(ACTIVITY["master"]["context_matrix_tsv"]),
            stats_json=str(ACTIVITY["master"]["stats_json"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            validation=ACTIVITY_MASTER_VALIDATION
        params:
            expected=ACTIVITY["master"]
        resources:
            mem_mb=2000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/validation/master.log"
        script:
            "../scripts/validate_external_master.py"


    rule prepare_activity_atac_insertions:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            validated=lambda wc: ACTIVITY_BAM_VALIDATIONS[wc.library],
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            sorter=str(REPO_ROOT / "src" / "sort_bed_by_reference.sh")
        output:
            bed=f"{ACTIVITY_QUANTIFICATION_ROOT}/units/{{library}}.bed.gz",
            unit_count=f"{ACTIVITY_QUANTIFICATION_ROOT}/units/{{library}}.count.txt"
        params:
            maximum=int(ACTIVITY["atac_fragment_maximum"])
        wildcard_constraints:
            library=ACTIVITY_ATAC_LIBRARY_RE
        threads: 4
        resources:
            mem_mb=6000
        conda:
            "../envs/atac_qc.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/units/{{library}}.atac.log"
        shell:
            r"""
            mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.bed:q})/.{wildcards.library}.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            samtools view -@ {threads} -b -f 2 -F 3852 \
              -e 'tlen != 0 && tlen > -{params.maximum} && tlen < {params.maximum}' \
              -o "$temporary/short.bam" {input.bam:q} 2> {log:q}
            samtools index -@ {threads} "$temporary/short.bam" \
              "$temporary/short.bam.bai" 2>> {log:q}
            alignmentSieve --ATACshift -b "$temporary/short.bam" \
              -o "$temporary/shifted.bam" --numberOfProcessors {threads} \
              >> {log:q} 2>&1
            bedtools bamtobed -i "$temporary/shifted.bam" \
              | awk -v countfile="$temporary/count" 'BEGIN {{OFS="\t"}}
                  $6 == "+" {{print $1,$2,$2+1,$4,$5,$6; n++; next}}
                  $6 == "-" {{print $1,$3-1,$3,$4,$5,$6; n++; next}}
                  END {{print n+0 > countfile}}' \
              | bash {input.sorter:q} {input.chrom_sizes:q} \
                  "$temporary/sort" 2>> {log:q} \
              | pigz -p 2 -c > "$temporary/units.bed.gz"
            test "$(samtools view -c "$temporary/shifted.bam")" \
              -eq "$(cat "$temporary/count")"
            pigz -t "$temporary/units.bed.gz"
            mv "$temporary/units.bed.gz" {output.bed:q}
            mv "$temporary/count" {output.unit_count:q}
            """


    rule prepare_activity_h3k27ac_fragments:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            validated=lambda wc: ACTIVITY_BAM_VALIDATIONS[wc.library],
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            sorter=str(REPO_ROOT / "src" / "sort_bed_by_reference.sh")
        output:
            bed=f"{ACTIVITY_QUANTIFICATION_ROOT}/units/{{library}}.bed.gz",
            unit_count=f"{ACTIVITY_QUANTIFICATION_ROOT}/units/{{library}}.count.txt"
        wildcard_constraints:
            library=ACTIVITY_H3K27AC_PAIRED_LIBRARY_RE
        threads: 4
        resources:
            mem_mb=4000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/units/{{library}}.h3k27ac.log"
        shell:
            r"""
            mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.bed:q})/.{wildcards.library}.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            samtools view -@ {threads} -f 2 -F 3852 -e 'tlen > 0' \
              {input.bam:q} 2> {log:q} \
              | awk -v countfile="$temporary/count" 'BEGIN {{OFS="\t"}}
                  {{start=$4-1; print $3,start,start+$9,$1,0,"."; n++}}
                  END {{print n+0 > countfile}}' \
              | bash {input.sorter:q} {input.chrom_sizes:q} \
                  "$temporary/sort" 2>> {log:q} \
              | pigz -p 2 -c > "$temporary/units.bed.gz"
            pigz -t "$temporary/units.bed.gz"
            mv "$temporary/units.bed.gz" {output.bed:q}
            mv "$temporary/count" {output.unit_count:q}
            """


    rule prepare_activity_h3k27ac_single_fragments:
        input:
            bam=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bam"],
            bai=lambda wc: ACTIVITY_LIBRARIES[wc.library]["bai"],
            validated=lambda wc: ACTIVITY_BAM_VALIDATIONS[wc.library],
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            sorter=str(REPO_ROOT / "src" / "sort_bed_by_reference.sh"),
            implementation=str(REPO_ROOT / "src" / "extend_single_end_fragments.py")
        output:
            bed=f"{ACTIVITY_QUANTIFICATION_ROOT}/units/{{library}}.bed.gz",
            unit_count=f"{ACTIVITY_QUANTIFICATION_ROOT}/units/{{library}}.count.txt"
        params:
            fragment_length=lambda wc: int(
                ACTIVITY_LIBRARIES[wc.library]["estimated_fragment_length_bp"]
            )
        wildcard_constraints:
            library=ACTIVITY_H3K27AC_SINGLE_LIBRARY_RE
        threads: 4
        resources:
            mem_mb=4000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/units/{{library}}.h3k27ac-single.log"
        shell:
            r"""
            mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.bed:q})/.{wildcards.library}.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            samtools view -@ {threads} -b -F 3844 {input.bam:q} 2> {log:q} \
              | bedtools bamtobed -i - 2>> {log:q} \
              | python {input.implementation:q} \
                  --chrom-sizes {input.chrom_sizes:q} \
                  --fragment-length {params.fragment_length} \
                  --count-output "$temporary/count" 2>> {log:q} \
              | bash {input.sorter:q} {input.chrom_sizes:q} \
                  "$temporary/sort" 2>> {log:q} \
              | pigz -p 2 -c > "$temporary/units.bed.gz"
            pigz -t "$temporary/units.bed.gz"
            test "$(pigz -dc "$temporary/units.bed.gz" | wc -l)" \
              -eq "$(cat "$temporary/count")"
            mv "$temporary/units.bed.gz" {output.bed:q}
            mv "$temporary/count" {output.unit_count:q}
            """


    rule count_activity_library:
        input:
            units=lambda wc: ACTIVITY_UNIT_BEDS[wc.library],
            total=lambda wc: ACTIVITY_UNIT_COUNTS[wc.library],
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            master_validated=ACTIVITY_MASTER_VALIDATION,
            script=str(REPO_ROOT / "src" / "count_activity_units.py"),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "activity.py")
        output:
            signal=f"{ACTIVITY_QUANTIFICATION_ROOT}/libraries/{{library}}.signal.tsv.gz",
            summary=f"{ACTIVITY_QUANTIFICATION_ROOT}/libraries/{{library}}.summary.json"
        params:
            assay=lambda wc: ACTIVITY_LIBRARIES[wc.library]["assay"],
            context=lambda wc: ACTIVITY_LIBRARIES[wc.library]["context"]
        wildcard_constraints:
            library=ACTIVITY_LIBRARY_RE
        resources:
            mem_mb=4000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/counts/{{library}}.log"
        shell:
            "mkdir -p $(dirname {output.signal:q}) $(dirname {log:q}) && "
            "python {input.script:q} --master-bed {input.master:q} "
            "--summit-bed {input.summits:q} --units-bed {input.units:q} "
            "--chrom-sizes {input.chrom_sizes:q} "
            "--total-units {input.total:q} --library-id {wildcards.library:q} "
            "--assay {params.assay:q} --cohort atlas "
            "--context {params.context:q} --output {output.signal:q} "
            "--summary {output.summary:q} > {log:q} 2>&1"


    rule build_activity_background_bins:
        input:
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            autosomes=str(REFERENCE["autosomes_file"]),
            script=str(REPO_ROOT / "src" / "make_activity_background_bins.py"),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "activity_tmm.py")
        output:
            bins=ACTIVITY_BACKGROUND_BINS
        params:
            bin_width=10000
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/background-bins.log"
        shell:
            "mkdir -p $(dirname {output.bins:q}) $(dirname {log:q}) && "
            "python {input.script:q} --chrom-sizes {input.chrom_sizes:q} "
            "--autosomes {input.autosomes:q} --bin-width {params.bin_width} "
            "--output {output.bins:q} > {log:q} 2>&1"


    rule count_activity_background_library:
        input:
            units=lambda wc: ACTIVITY_UNIT_BEDS[wc.library],
            bins=ACTIVITY_BACKGROUND_BINS,
            script=str(REPO_ROOT / "src" / "count_activity_background.sh")
        output:
            counts=temp(f"{ACTIVITY_WORK}/background_10kb/{{library}}.counts.tsv.gz")
        wildcard_constraints:
            library=ACTIVITY_LIBRARY_RE
        threads: 2
        resources:
            mem_mb=2000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/background/{{library}}.log"
        shell:
            "mkdir -p $(dirname {log:q}) && "
            "bash {input.script:q} {input.units:q} {input.bins:q} "
            "{output.counts:q} {threads} > {log:q} 2>&1"


    rule build_activity_tmm_inputs:
        input:
            signals=list(ACTIVITY_LIBRARY_SIGNALS.values()),
            background=list(ACTIVITY_BACKGROUND_COUNTS.values()),
            script=str(REPO_ROOT / "src" / "build_activity_tmm_inputs.py"),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "activity_tmm.py")
        output:
            counts=ACTIVITY_TMM_COUNTS,
            metadata=ACTIVITY_TMM_METADATA
        params:
            method=TMM_BACKGROUND_METHOD,
            signal_arguments=" ".join(
                "--signal {}".format(
                    shlex.quote(f"{library}={ACTIVITY_LIBRARY_SIGNALS[library]}")
                )
                for library in ACTIVITY_LIBRARY_IDS
            ),
            background_arguments=" ".join(
                "--background-count {}".format(
                    shlex.quote(f"{library}={ACTIVITY_BACKGROUND_COUNTS[library]}")
                )
                for library in ACTIVITY_LIBRARY_IDS
            )
        resources:
            mem_mb=4000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/build-tmm-inputs.log"
        shell:
            "mkdir -p $(dirname {output.counts:q}) $(dirname {log:q}) && "
            "python {input.script:q} --method {params.method:q} "
            "{params.signal_arguments} {params.background_arguments} "
            "--output-counts {output.counts:q} "
            "--output-metadata {output.metadata:q} > {log:q} 2>&1"


    rule calculate_activity_tmm_factors:
        input:
            counts=ACTIVITY_TMM_COUNTS,
            metadata=ACTIVITY_TMM_METADATA,
            script=str(REPO_ROOT / "src" / "calculate_activity_tmm.R")
        output:
            factors=ACTIVITY_TMM_FACTORS,
            receipt=ACTIVITY_TMM_RECEIPT
        params:
            method=TMM_BACKGROUND_METHOD
        resources:
            mem_mb=4000
        conda:
            "../envs/activity_tmm.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/calculate-tmm-factors.log"
        shell:
            "mkdir -p $(dirname {output.factors:q}) $(dirname {log:q}) && "
            "Rscript {input.script:q} --method {params.method:q} "
            "--counts {input.counts:q} --metadata {input.metadata:q} "
            "--output-factors {output.factors:q} "
            "--output-receipt {output.receipt:q} > {log:q} 2>&1"


    rule build_activity_tmm_table:
        input:
            factors=ACTIVITY_TMM_FACTORS,
            receipt=ACTIVITY_TMM_RECEIPT,
            counts=ACTIVITY_TMM_COUNTS,
            metadata=ACTIVITY_TMM_METADATA,
            signals=list(ACTIVITY_LIBRARY_SIGNALS.values()),
            master_validated=ACTIVITY_MASTER_VALIDATION,
            bam_validations=list(ACTIVITY_BAM_VALIDATIONS.values()),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "activity_tmm.py")
        output:
            context_signal=ACTIVITY_CONTEXT_SIGNAL,
            activity=ACTIVITY_TABLE,
            metrics=ACTIVITY_METRICS,
            provenance=ACTIVITY_PROVENANCE
        params:
            signal_paths=ACTIVITY_LIBRARY_SIGNALS,
            contexts=ACTIVITY_CONTEXTS,
            activity_config=ACTIVITY,
            workflow_provenance=config.get("provenance", {})
        resources:
            mem_mb=8000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/build-table.log"
        script:
            "../scripts/build_activity_tmm_table.py"


    rule build_activity_h3k27ac_windows:
        input:
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            master_validated=ACTIVITY_MASTER_VALIDATION,
            script=str(REPO_ROOT / "src" / "make_regulatory_windows.py"),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "regulatory_elements.py")
        output:
            table=ACTIVITY_H3K27AC_WINDOWS,
            bed=ACTIVITY_H3K27AC_WINDOWS_BED
        resources:
            mem_mb=2000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/h3k27ac-windows.log"
        shell:
            "mkdir -p $(dirname {output.table:q}) $(dirname {log:q}) && "
            "python {input.script:q} --master-bed {input.master:q} "
            "--summit-bed {input.summits:q} "
            "--chrom-sizes {input.chrom_sizes:q} "
            "--output-table {output.table:q} --output-bed {output.bed:q} "
            "> {log:q} 2>&1"


    rule count_activity_h3k27ac_windows:
        input:
            units=lambda wc: ACTIVITY_UNIT_BEDS[wc.library],
            total=lambda wc: ACTIVITY_UNIT_COUNTS[wc.library],
            windows=ACTIVITY_H3K27AC_WINDOWS,
            windows_bed=ACTIVITY_H3K27AC_WINDOWS_BED,
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            script=str(REPO_ROOT / "src" / "format_regulatory_window_counts.py"),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "regulatory_elements.py")
        output:
            counts=f"{ACTIVITY_H3K27AC_ROOT}/libraries/{{library}}.window_counts.tsv.gz"
        params:
            context=lambda wc: ACTIVITY_LIBRARIES[wc.library]["context"]
        wildcard_constraints:
            library=wildcard_regex(ACTIVITY_H3K27AC_LIBRARIES)
        threads: 2
        resources:
            mem_mb=4000
        conda:
            "../envs/alignment.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/{{library}}.h3k27ac-windows.log"
        shell:
            r"""
            mkdir -p $(dirname {output.counts:q}) $(dirname {log:q})
            temporary=$(mktemp -d $(dirname {output.counts:q})/.{wildcards.library}.XXXXXX)
            trap 'rm -rf "$temporary"' EXIT
            pigz -dc {input.units:q} \
              | bedtools coverage -a {input.windows_bed:q} -b stdin \
                  -counts -sorted -g {input.chrom_sizes:q} \
                  > "$temporary/coverage.tsv" 2> {log:q}
            python {input.script:q} --window-table {input.windows:q} \
              --coverage "$temporary/coverage.tsv" \
              --total-units {input.total:q} \
              --library-id {wildcards.library:q} \
              --context {params.context:q} --output {output.counts:q} \
              >> {log:q} 2>&1
            pigz -t {output.counts:q}
            """


    rule build_activity_h3k27ac_signal:
        input:
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            context_matrix=str(ACTIVITY["master"]["context_matrix_tsv"]),
            blacklist=str(REFERENCE["blacklist_bed"]),
            windows=ACTIVITY_H3K27AC_WINDOWS,
            counts=list(ACTIVITY_H3K27AC_WINDOW_COUNTS.values()),
            factors=ACTIVITY_TMM_FACTORS,
            master_validated=ACTIVITY_MASTER_VALIDATION,
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "regulatory_elements.py")
        output:
            signal=ACTIVITY_H3K27AC_SIGNAL,
            metrics=ACTIVITY_H3K27AC_METRICS,
            provenance=ACTIVITY_H3K27AC_PROVENANCE
        params:
            window_counts=ACTIVITY_H3K27AC_WINDOW_COUNTS,
            contexts=ACTIVITY_CONTEXTS
        resources:
            mem_mb=12000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/quantification/h3k27ac-signal.log"
        script:
            "../scripts/build_activity_h3k27ac_signal.py"


    rule build_regulatory_element_catalog:
        input:
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            context_matrix=str(ACTIVITY["master"]["context_matrix_tsv"]),
            tss=str(REFERENCE["tss_bed"]),
            blacklist=str(REFERENCE["blacklist_bed"]),
            h3k27ac_signal=ACTIVITY_H3K27AC_SIGNAL,
            activity=ACTIVITY_TABLE,
            master_validated=ACTIVITY_MASTER_VALIDATION,
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "regulatory_elements.py")
        output:
            catalog=ACTIVITY_REGULATORY_CATALOG,
            wide=ACTIVITY_REGULATORY_WIDE,
            elements=list(ACTIVITY_REGULATORY_ELEMENTS.values()),
            summary=ACTIVITY_REGULATORY_SUMMARY,
            metrics=ACTIVITY_REGULATORY_METRICS,
            provenance=ACTIVITY_REGULATORY_PROVENANCE
        params:
            element_paths=ACTIVITY_REGULATORY_ELEMENTS,
            contexts=ACTIVITY_CONTEXTS
        resources:
            mem_mb=12000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/catalog/catalog.log"
        script:
            "../scripts/build_regulatory_element_catalog.py"


    rule build_catalog_bed_tracks:
        input:
            master=str(ACTIVITY["master"]["master_bed"]),
            summits=str(ACTIVITY["master"]["summits_bed"]),
            context_matrix=str(ACTIVITY["master"]["context_matrix_tsv"]),
            elements=list(ACTIVITY_REGULATORY_ELEMENTS.values()),
            catalog_metrics=ACTIVITY_REGULATORY_METRICS,
            implementation=str(
                REPO_ROOT
                / "src"
                / "short_read_processing"
                / "catalog_visualization.py"
            )
        output:
            master=ACTIVITY_CATALOG_MASTER_BED,
            context_dhs=list(ACTIVITY_CATALOG_CONTEXT_DHS.values()),
            element_beds=list(ACTIVITY_CATALOG_ELEMENT_BEDS.values()),
            manifest=ACTIVITY_CATALOG_BED_MANIFEST
        params:
            element_paths=ACTIVITY_REGULATORY_ELEMENTS,
            context_dhs_paths=ACTIVITY_CATALOG_CONTEXT_DHS,
            element_bed_paths=ACTIVITY_CATALOG_ELEMENT_BEDS
        resources:
            mem_mb=3000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/catalog/bed-tracks.log"
        script:
            "../scripts/build_catalog_bed_tracks.py"


    rule build_context_mean_bigwig:
        input:
            units=lambda wc: [
                ACTIVITY_UNIT_BEDS[library]
                for library in ACTIVITY_CONTEXT_ASSAY_LIBRARIES[
                    (wc.context, wc.assay)
                ]
            ],
            factors=ACTIVITY_TMM_FACTORS,
            chrom_sizes=str(REFERENCE["chrom_sizes"]),
            implementation=str(
                REPO_ROOT
                / "src"
                / "short_read_processing"
                / "catalog_visualization.py"
            )
        output:
            bigwig=(
                f"{ACTIVITY_CATALOG_TRACK_ROOT}/"
                "{context}.{assay}.mean.background_tmm.bw"
            ),
            metrics=(
                f"{ACTIVITY_CATALOG_TRACK_ROOT}/"
                "{context}.{assay}.mean.background_tmm.json"
            )
        params:
            context=lambda wc: wc.context,
            assay=lambda wc: wc.assay,
            atac_extension_bp=int(ACTIVITY["atac_browser_extension_bp"]),
            unit_paths=lambda wc: {
                library: ACTIVITY_UNIT_BEDS[library]
                for library in ACTIVITY_CONTEXT_ASSAY_LIBRARIES[
                    (wc.context, wc.assay)
                ]
            }
        wildcard_constraints:
            context=ACTIVITY_CONTEXT_RE,
            assay="atac|h3k27ac"
        threads: 2
        resources:
            mem_mb=5000
        conda:
            "../envs/catalog_tracks.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/catalog/tracks/{{context}}.{{assay}}.log"
        script:
            "../scripts/build_context_mean_bigwig.py"


    rule build_context_igv_session:
        input:
            atac=lambda wc: ACTIVITY_CONTEXT_MEAN_BIGWIGS[(wc.context, "atac")],
            h3k27ac=lambda wc: ACTIVITY_CONTEXT_MEAN_BIGWIGS[
                (wc.context, "h3k27ac")
            ],
            context_dhs=lambda wc: ACTIVITY_CATALOG_CONTEXT_DHS[wc.context],
            master_dhs=ACTIVITY_CATALOG_MASTER_BED,
            elements=lambda wc: ACTIVITY_CATALOG_ELEMENT_BEDS[wc.context],
            implementation=str(REPO_ROOT / "src" / "build_igv_session.py")
        output:
            session=f"{ACTIVITY_CATALOG_IGV_ROOT}/{{context}}.xml"
        params:
            context=lambda wc: wc.context,
            genome=str(REFERENCE["name"])
        wildcard_constraints:
            context=ACTIVITY_CONTEXT_RE
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/catalog/igv/{{context}}.log"
        script:
            "../scripts/build_catalog_igv_session.py"


    rule build_all_contexts_igv_session:
        input:
            atac=[
                ACTIVITY_CONTEXT_MEAN_BIGWIGS[(context, "atac")]
                for context in ACTIVITY_CONTEXTS
            ],
            h3k27ac=[
                ACTIVITY_CONTEXT_MEAN_BIGWIGS[(context, "h3k27ac")]
                for context in ACTIVITY_CONTEXTS
            ],
            context_dhs=list(ACTIVITY_CATALOG_CONTEXT_DHS.values()),
            master_dhs=ACTIVITY_CATALOG_MASTER_BED,
            elements=list(ACTIVITY_CATALOG_ELEMENT_BEDS.values()),
            implementation=str(REPO_ROOT / "src" / "build_igv_session.py")
        output:
            session=ACTIVITY_ALL_CONTEXTS_IGV_SESSION
        params:
            contexts=ACTIVITY_CONTEXTS,
            genome=str(REFERENCE["name"]),
            atac_paths={
                context: ACTIVITY_CONTEXT_MEAN_BIGWIGS[(context, "atac")]
                for context in ACTIVITY_CONTEXTS
            },
            h3k27ac_paths={
                context: ACTIVITY_CONTEXT_MEAN_BIGWIGS[(context, "h3k27ac")]
                for context in ACTIVITY_CONTEXTS
            },
            context_dhs_paths=ACTIVITY_CATALOG_CONTEXT_DHS,
            element_paths=ACTIVITY_CATALOG_ELEMENT_BEDS
        resources:
            mem_mb=1000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/catalog/igv/all-contexts.log"
        script:
            "../scripts/build_all_contexts_catalog_igv_session.py"


    rule build_integrated_qc_report:
        input:
            sources=ACTIVITY_REPORT_SOURCE_FILES,
            activity_metrics=ACTIVITY_METRICS,
            activity_provenance=ACTIVITY_PROVENANCE,
            factors=ACTIVITY_TMM_FACTORS,
            tmm_receipt=ACTIVITY_TMM_RECEIPT,
            context_signal=ACTIVITY_CONTEXT_SIGNAL,
            activity_table=ACTIVITY_TABLE,
            master_metrics=str(ACTIVITY["master"]["stats_json"]),
            catalog=ACTIVITY_REGULATORY_CATALOG,
            wide=ACTIVITY_REGULATORY_WIDE,
            elements=list(ACTIVITY_REGULATORY_ELEMENTS.values()),
            summary=ACTIVITY_REGULATORY_SUMMARY,
            catalog_metrics=ACTIVITY_REGULATORY_METRICS,
            catalog_provenance=ACTIVITY_REGULATORY_PROVENANCE,
            catalog_master_bed=ACTIVITY_CATALOG_MASTER_BED,
            context_dhs_beds=list(ACTIVITY_CATALOG_CONTEXT_DHS.values()),
            element_beds=list(ACTIVITY_CATALOG_ELEMENT_BEDS.values()),
            bed_track_manifest=ACTIVITY_CATALOG_BED_MANIFEST,
            mean_bigwigs=list(ACTIVITY_CONTEXT_MEAN_BIGWIGS.values()),
            mean_bigwig_metrics=list(
                ACTIVITY_CONTEXT_MEAN_BIGWIG_METRICS.values()
            ),
            igv_sessions=list(ACTIVITY_CONTEXT_IGV_SESSIONS.values()),
            all_contexts_igv_session=ACTIVITY_ALL_CONTEXTS_IGV_SESSION,
            nearest_tss_candidates=NEAREST_TSS_ENHANCER_CANDIDATES,
            nearest_tss_metrics=NEAREST_TSS_ENHANCER_METRICS,
            contact_link_metrics=(CONTACT_LINK_METRICS if CONTACTS else []),
            contact_link_provenance=(CONTACT_LINK_PROVENANCE if CONTACTS else []),
            implementation=str(REPO_ROOT / "src" / "short_read_processing" / "integrated_report.py")
        output:
            html=ACTIVITY_REPORT_HTML,
            pdf=ACTIVITY_REPORT_PDF,
            metrics=ACTIVITY_REPORT_METRICS
        params:
            workflow_config=config,
            report_sources=config["report"]["source_files"],
            current_files={
                "activity_metrics": ACTIVITY_METRICS,
                "activity_provenance": ACTIVITY_PROVENANCE,
                "normalization_factors": ACTIVITY_TMM_FACTORS,
                "tmm_software": ACTIVITY_TMM_RECEIPT,
                "context_signal": ACTIVITY_CONTEXT_SIGNAL,
                "master_dhs_activity": ACTIVITY_TABLE,
                "master_elements_long": ACTIVITY_REGULATORY_CATALOG,
                "master_elements_wide": ACTIVITY_REGULATORY_WIDE,
                "regulatory_element_summary": ACTIVITY_REGULATORY_SUMMARY,
                "regulatory_element_metrics": ACTIVITY_REGULATORY_METRICS,
                "regulatory_element_provenance": ACTIVITY_REGULATORY_PROVENANCE,
                "catalog_master_dhs_bed": ACTIVITY_CATALOG_MASTER_BED,
                "catalog_bed_manifest": ACTIVITY_CATALOG_BED_MANIFEST,
                **{
                    f"context_elements_{context}": path
                    for context, path in ACTIVITY_REGULATORY_ELEMENTS.items()
                },
                **{
                    f"context_dhs_bed_{context}": path
                    for context, path in ACTIVITY_CATALOG_CONTEXT_DHS.items()
                },
                **{
                    f"context_elements_bed_{context}": path
                    for context, path in ACTIVITY_CATALOG_ELEMENT_BEDS.items()
                },
                **{
                    f"mean_bigwig_{context}_{assay}": path
                    for (context, assay), path in ACTIVITY_CONTEXT_MEAN_BIGWIGS.items()
                },
                **{
                    f"mean_bigwig_metrics_{context}_{assay}": path
                    for (context, assay), path in ACTIVITY_CONTEXT_MEAN_BIGWIG_METRICS.items()
                },
                **{
                    f"igv_session_{context}": path
                    for context, path in ACTIVITY_CONTEXT_IGV_SESSIONS.items()
                },
                "igv_session_all_contexts": ACTIVITY_ALL_CONTEXTS_IGV_SESSION,
                "nearest_tss_enhancer_candidates": NEAREST_TSS_ENHANCER_CANDIDATES,
                "nearest_tss_enhancer_metrics": NEAREST_TSS_ENHANCER_METRICS,
                **(
                    {
                        "contact_graph_metrics": CONTACT_LINK_METRICS,
                        "contact_graph_provenance": CONTACT_LINK_PROVENANCE,
                    }
                    if CONTACTS
                    else {}
                ),
            }
        resources:
            mem_mb=4000
        conda:
            "../envs/final_report.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/report/integrated-qc.log"
        script:
            "../scripts/build_integrated_qc_report.py"
