"""Build one context's weighted element-promoter graph and gene projection."""

from pathlib import Path

from short_read_processing.contact_links import build_context_links


observed = list(snakemake.input.observed)
metrics = build_context_links(
    context=str(snakemake.wildcards.context),
    context_elements=Path(str(snakemake.input.elements)),
    promoters_path=Path(str(snakemake.input.promoters)),
    powerlaw_path=Path(str(snakemake.input.powerlaw)),
    strategy=str(snakemake.params.strategy),
    contact_assay=str(snakemake.params.assay),
    contact_match=str(snakemake.params.match),
    configured_resolution=int(snakemake.params.resolution),
    maximum_distance=int(snakemake.params.maximum_distance),
    pseudocount_fraction=float(snakemake.params.pseudocount_fraction),
    observed_contact_path=Path(str(observed[0])) if observed else None,
    nodes_output=Path(str(snakemake.output.nodes)),
    edges_output=Path(str(snakemake.output.edges)),
    genes_output=Path(str(snakemake.output.genes)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built {metrics['element_promoter_edge_count']} element-promoter edges and "
    f"{metrics['element_gene_candidate_count']} element-gene candidates for "
    f"{metrics['context']}\n",
    encoding="utf-8",
)
