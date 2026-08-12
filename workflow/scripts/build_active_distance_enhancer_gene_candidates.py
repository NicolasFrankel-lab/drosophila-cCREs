"""Derive active-promoter candidates for one distance-model context."""

from pathlib import Path

from short_read_processing.contact_links import (
    build_active_distance_enhancer_gene_candidates,
)


metrics = build_active_distance_enhancer_gene_candidates(
    context=str(snakemake.wildcards.context),
    nodes_path=Path(str(snakemake.input.nodes)),
    edges_path=Path(str(snakemake.input.edges)),
    output=Path(str(snakemake.output.candidates)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Retained {metrics['active_distance_enhancer_gene_candidate_count']} "
    f"active-promoter distance candidates for {metrics['context']}\n",
    encoding="utf-8",
)
