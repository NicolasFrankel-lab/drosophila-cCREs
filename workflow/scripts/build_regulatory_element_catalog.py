"""Build the background-TMM, max-three-window regulatory-element catalog."""

from pathlib import Path

from short_read_processing.regulatory_elements import build_regulatory_catalog


metrics = build_regulatory_catalog(
    master_bed=Path(str(snakemake.input.master)),
    summit_bed=Path(str(snakemake.input.summits)),
    context_matrix=Path(str(snakemake.input.context_matrix)),
    tss_bed=Path(str(snakemake.input.tss)),
    blacklist_bed=Path(str(snakemake.input.blacklist)),
    h3k27ac_quantification=Path(str(snakemake.input.h3k27ac_signal)),
    activity_table=Path(str(snakemake.input.activity)),
    contexts=list(snakemake.params.contexts),
    output_catalog=Path(str(snakemake.output.catalog)),
    output_wide=Path(str(snakemake.output.wide)),
    output_element_paths={
        str(context): Path(str(path))
        for context, path in dict(snakemake.params.element_paths).items()
    },
    output_summary=Path(str(snakemake.output.summary)),
    output_metrics=Path(str(snakemake.output.metrics)),
    output_provenance=Path(str(snakemake.output.provenance)),
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    f"Built {metrics['catalog_row_count']} regulatory-element context rows\n",
    encoding="utf-8",
)
