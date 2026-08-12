"""Download, normalize, model, and link context-resolved chromatin contacts."""


if CONTACTS:
    rule download_contact_mcool:
        output:
            "data/raw/contacts/{source_id}.mcool"
        params:
            url=lambda wc: CONTACT_SOURCE_BY_ID[wc.source_id]["url"],
            checksum=lambda wc: CONTACT_DOWNLOAD_CHECKSUM_ARGUMENTS[wc.source_id]
        wildcard_constraints:
            source_id=wildcard_regex(CONTACT_SOURCE_IDS_BY_FORMAT["mcool"])
        resources:
            mem_mb=1000,
            contact_download_slots=1
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/contacts/download-{{source_id}}.log"
        shell:
            "mkdir -p $(dirname {output:q}) $(dirname {log:q}) && "
            "aria2c --continue=true --allow-overwrite=true "
            "--auto-file-renaming=false --file-allocation=none "
            "--max-connection-per-server=1 --split=1 --min-split-size=5M "
            "--check-integrity=true {params.checksum} "
            "--dir=$(dirname {output:q}) --out=$(basename {output:q}) {params.url:q} "
            "> {log:q} 2>&1"


    rule download_contact_cool_gz:
        output:
            "data/raw/contacts/{source_id}.cool.gz"
        params:
            url=lambda wc: CONTACT_SOURCE_BY_ID[wc.source_id]["url"],
            checksum=lambda wc: CONTACT_DOWNLOAD_CHECKSUM_ARGUMENTS[wc.source_id]
        wildcard_constraints:
            source_id=wildcard_regex(CONTACT_SOURCE_IDS_BY_FORMAT["cool.gz"])
        resources:
            mem_mb=1000,
            contact_download_slots=1
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/contacts/download-{{source_id}}.log"
        shell:
            "mkdir -p $(dirname {output:q}) $(dirname {log:q}) && "
            "aria2c --continue=true --allow-overwrite=true "
            "--auto-file-renaming=false --file-allocation=none "
            "--max-connection-per-server=1 --split=1 --min-split-size=5M "
            "--check-integrity=true {params.checksum} "
            "--dir=$(dirname {output:q}) --out=$(basename {output:q}) {params.url:q} "
            "> {log:q} 2>&1"


    rule download_contact_h5:
        output:
            "data/raw/contacts/{source_id}.h5"
        params:
            url=lambda wc: CONTACT_SOURCE_BY_ID[wc.source_id]["url"],
            checksum=lambda wc: CONTACT_DOWNLOAD_CHECKSUM_ARGUMENTS[wc.source_id]
        wildcard_constraints:
            source_id=wildcard_regex(CONTACT_SOURCE_IDS_BY_FORMAT["h5"])
        resources:
            mem_mb=1000,
            contact_download_slots=1
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/contacts/download-{{source_id}}.log"
        shell:
            "mkdir -p $(dirname {output:q}) $(dirname {log:q}) && "
            "aria2c --continue=true --allow-overwrite=true "
            "--auto-file-renaming=false --file-allocation=none "
            "--max-connection-per-server=1 --split=1 --min-split-size=5M "
            "--check-integrity=true {params.checksum} "
            "--dir=$(dirname {output:q}) --out=$(basename {output:q}) {params.url:q} "
            "> {log:q} 2>&1"


    rule standardize_context_contacts:
        input:
            manifest=CONTACT_SOURCE_MANIFEST,
            sources=lambda wc: [
                row["local_path"]
                for row in CONTACT_SOURCE_ROWS_BY_CONTEXT[wc.context]
            ]
        output:
            cool=f"{CONTACT_NORMALIZED_ROOT}/{{context}}.balanced.cool",
            metrics=f"{CONTACT_NORMALIZED_ROOT}/{{context}}.metrics.json"
        params:
            repository_root=str(REPO_ROOT),
            resolution=lambda wc: int(
                CONTACT_CONTEXT_BY_ID[wc.context]["resolution_bp"]
            ),
            workdir=lambda wc: f"{CONTACT_WORK_ROOT}/{wc.context}"
        wildcard_constraints:
            context=wildcard_regex(CONTACT_OBSERVED_CONTEXTS)
        threads: 1
        resources:
            mem_mb=16000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/contacts/{{context}}.normalize.log"
        script:
            "../scripts/standardize_contacts.py"


    rule fit_contact_powerlaw:
        input:
            contacts=list(CONTACT_NORMALIZED.values())
        output:
            powerlaw=CONTACT_POWERLAW
        params:
            contacts=CONTACT_NORMALIZED,
            canonical_chromosomes=CONTACTS["canonical_chromosomes"],
            maximum_distance=int(CONTACTS["maximum_distance_bp"])
        threads: 1
        resources:
            mem_mb=16000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/contacts/fit-powerlaw.log"
        script:
            "../scripts/fit_contact_powerlaw.py"


    rule build_contact_promoters:
        input:
            annotation=str(CONTACTS["promoter_annotation"]),
            chrom_sizes=str(REFERENCE["chrom_sizes"])
        output:
            promoters=CONTACT_PROMOTERS,
            metrics=CONTACT_PROMOTER_METRICS
        params:
            canonical_chromosomes=CONTACTS["canonical_chromosomes"],
            promoter_width=int(CONTACTS["promoter_width_bp"]),
            annotation_checksum=CONTACTS["promoter_annotation_checksum"]
        resources:
            mem_mb=4000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/promoters.log"
        script:
            "../scripts/build_contact_promoters.py"


    rule build_context_contact_links:
        input:
            elements=ACTIVITY_REGULATORY_CATALOG,
            promoters=CONTACT_PROMOTERS,
            powerlaw=CONTACT_POWERLAW,
            observed=lambda wc: (
                [CONTACT_NORMALIZED[wc.context]]
                if wc.context in CONTACT_NORMALIZED
                else []
            )
        output:
            nodes=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.nodes.tsv.gz",
            edges=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.element_promoter_edges.tsv.gz",
            genes=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.element_gene_candidates.tsv.gz",
            metrics=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.metrics.json"
        params:
            strategy=lambda wc: CONTACT_CONTEXT_BY_ID[wc.context]["strategy"],
            assay=lambda wc: CONTACT_CONTEXT_BY_ID[wc.context]["assay"],
            match=lambda wc: CONTACT_CONTEXT_BY_ID[wc.context]["match"],
            resolution=lambda wc: int(
                CONTACT_CONTEXT_BY_ID[wc.context]["resolution_bp"]
            ),
            maximum_distance=int(CONTACTS["maximum_distance_bp"]),
            pseudocount_fraction=float(CONTACTS["pseudocount_fraction"]),
        wildcard_constraints:
            context=wildcard_regex(CONTACT_CONTEXTS)
        threads: 1
        resources:
            mem_mb=12000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/{{context}}.links.log"
        script:
            "../scripts/build_context_contact_links.py"


    rule build_active_contact_enhancer_gene_candidates:
        input:
            nodes=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.nodes.tsv.gz",
            edges=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.element_promoter_edges.tsv.gz"
        output:
            candidates=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.active_contact_enhancer_gene_candidates.tsv.gz",
            metrics=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.active_contact_enhancer_gene_candidates.metrics.json"
        params:
            observed_over_expected_threshold=float(
                CONTACTS.get("candidate_observed_over_expected_threshold", 1.0)
            )
        wildcard_constraints:
            context=wildcard_regex(CONTACT_CONTEXTS)
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/{{context}}.active-contact-candidates.log"
        script:
            "../scripts/build_active_contact_enhancer_gene_candidates.py"


    rule build_nearest_active_promoter_gene_candidates:
        input:
            nodes=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.nodes.tsv.gz"
        output:
            candidates=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.nearest_active_promoter_gene_candidates.tsv.gz",
            metrics=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.nearest_active_promoter_gene_candidates.metrics.json"
        wildcard_constraints:
            context=wildcard_regex(CONTACT_CONTEXTS)
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/{{context}}.nearest-active-promoter.log"
        script:
            "../scripts/build_nearest_active_promoter_gene_candidates.py"


    rule build_active_distance_enhancer_gene_candidates:
        input:
            nodes=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.nodes.tsv.gz",
            edges=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.element_promoter_edges.tsv.gz"
        output:
            candidates=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.active_distance_enhancer_gene_candidates.tsv.gz",
            metrics=f"{CONTACT_LINK_ROOT}/contexts/{{context}}.active_distance_enhancer_gene_candidates.metrics.json"
        wildcard_constraints:
            context=wildcard_regex(CONTACT_POWERLAW_CONTEXTS)
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/{{context}}.active-distance-candidates.log"
        script:
            "../scripts/build_active_distance_enhancer_gene_candidates.py"


    rule aggregate_contact_links:
        input:
            manifest=CONTACT_SOURCE_MANIFEST,
            promoter_metrics=CONTACT_PROMOTER_METRICS,
            contact_metrics=list(CONTACT_NORMALIZED_METRICS.values()),
            powerlaw=CONTACT_POWERLAW,
            context_metrics=list(CONTACT_LINK_CONTEXT_METRICS.values()),
            candidate_metrics=list(
                CONTACT_LINK_CONTEXT_ACTIVE_CONTACT_METRICS.values()
            ),
            nearest_candidate_metrics=list(
                CONTACT_LINK_CONTEXT_NEAREST_ACTIVE_PROMOTER_METRICS.values()
            ),
            distance_candidate_metrics=list(
                CONTACT_LINK_CONTEXT_ACTIVE_DISTANCE_METRICS.values()
            )
        output:
            metrics=CONTACT_LINK_METRICS,
            provenance=CONTACT_LINK_PROVENANCE
        resources:
            mem_mb=2000
        conda:
            "../envs/contacts.yaml"
        log:
            f"{RESULT_ROOT}/logs/activity/links/aggregate.log"
        script:
            "../scripts/aggregate_contact_links.py"
