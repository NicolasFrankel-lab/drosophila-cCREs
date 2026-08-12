"""Assign qualifying enhancers to the nearest supported active promoter TSS."""

from pathlib import Path

from short_read_processing.contact_links import (
    build_nearest_active_promoter_gene_candidates,
)


metrics = build_nearest_active_promoter_gene_candidates(
    context=str(snakemake.wildcards.context),
    nodes_path=Path(str(snakemake.input.nodes)),
    output=Path(str(snakemake.output.candidates)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Assigned {metrics['element_with_nearest_active_promoter_count']} "
    f"enhancers to nearest active promoter TSSs for {metrics['context']}\n",
    encoding="utf-8",
)
