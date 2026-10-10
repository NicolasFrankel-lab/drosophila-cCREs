"""Build the GTF promoter registry for nearest-TSS links."""

from pathlib import Path

from short_read_processing.nearest_tss_links import build_promoter_table


metrics = build_promoter_table(
    annotation=Path(str(snakemake.input.annotation)),
    chrom_sizes=Path(str(snakemake.input.chrom_sizes)),
    canonical_chromosomes=None,
    promoter_width=int(snakemake.params.promoter_width),
    promoter_id_prefix=str(snakemake.params.promoter_id_prefix),
    annotation_checksum=str(snakemake.params.annotation_checksum),
    output=Path(str(snakemake.output.promoters)),
    metrics_output=Path(str(snakemake.output.metrics)),
)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(
    f"Built {metrics['promoter_count']} nearest-TSS promoter nodes\n",
    encoding="utf-8",
)
