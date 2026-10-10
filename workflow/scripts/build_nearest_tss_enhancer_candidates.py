"""Build the wide nearest-TSS enhancer table."""

from pathlib import Path

from short_read_processing.nearest_tss_links import (
    build_nearest_tss_enhancer_candidates,
)


metrics = build_nearest_tss_enhancer_candidates(
    catalog_path=Path(str(snakemake.input.catalog)),
    promoters_path=Path(str(snakemake.input.promoters)),
    contexts=list(snakemake.params.contexts),
    enhancer_classes=list(snakemake.params.enhancer_classes),
    output=Path(str(snakemake.output.candidates)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Assigned {metrics['enhancer_count']} enhancers to nearest annotated TSSs\n",
    encoding="utf-8",
)
