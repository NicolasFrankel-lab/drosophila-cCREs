"""Derive one context's focused enhancer--gene candidates from graph tables."""

from pathlib import Path

from short_read_processing.contact_links import (
    build_active_contact_enhancer_gene_candidates,
)


metrics = build_active_contact_enhancer_gene_candidates(
    context=str(snakemake.wildcards.context),
    nodes_path=Path(str(snakemake.input.nodes)),
    edges_path=Path(str(snakemake.input.edges)),
    observed_over_expected_threshold=float(
        snakemake.params.observed_over_expected_threshold
    ),
    output=Path(str(snakemake.output.candidates)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Retained {metrics['active_contact_enhancer_gene_candidate_count']} "
    f"active-contact enhancer-gene candidates for {metrics['context']}\n",
    encoding="utf-8",
)
