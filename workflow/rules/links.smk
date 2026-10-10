"""Nearest-annotated-TSS enhancer links and promoter activity."""


if NEAREST_TSS_LINKS:
    rule build_nearest_tss_promoters:
        input:
            annotation=str(NEAREST_TSS_LINKS["promoter_annotation"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            promoters=NEAREST_TSS_PROMOTERS,
            metrics=NEAREST_TSS_PROMOTER_METRICS
        params:
            promoter_width=int(NEAREST_TSS_LINKS["promoter_width_bp"]),
            promoter_id_prefix=f"{str(REFERENCE['name']).upper()}PROM",
            annotation_checksum=NEAREST_TSS_LINKS[
                "promoter_annotation_checksum"
            ]
        resources:
            mem_mb=4000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/nearest-tss/promoters.log"
        script:
            "../scripts/build_nearest_tss_promoters.py"


    rule build_nearest_tss_enhancer_candidates:
        input:
            catalog=ACTIVITY_REGULATORY_CATALOG,
            promoters=NEAREST_TSS_PROMOTERS
        output:
            candidates=NEAREST_TSS_ENHANCER_CANDIDATES,
            metrics=NEAREST_TSS_ENHANCER_METRICS
        params:
            contexts=ACTIVITY_CONTEXTS,
            enhancer_classes=NEAREST_TSS_LINKS["enhancer_classes"]
        resources:
            mem_mb=8000
        conda:
            "../envs/reporting.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/nearest-tss/enhancers.log"
        script:
            "../scripts/build_nearest_tss_enhancer_candidates.py"
