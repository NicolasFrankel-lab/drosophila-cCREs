"""Build method-neutral summit-window H3K27ac quantification outputs."""

from pathlib import Path

from short_read_processing.regulatory_elements import build_h3k27ac_quantification


metrics = build_h3k27ac_quantification(
    master_bed=Path(str(snakemake.input.master)),
    summit_bed=Path(str(snakemake.input.summits)),
    context_matrix=Path(str(snakemake.input.context_matrix)),
    blacklist_bed=Path(str(snakemake.input.blacklist)),
    window_table=Path(str(snakemake.input.windows)),
    window_count_paths={
        str(library): Path(str(path))
        for library, path in dict(snakemake.params.window_counts).items()
    },
    factor_table=Path(str(snakemake.input.factors)),
    contexts=list(snakemake.params.contexts),
    output_signal=Path(str(snakemake.output.signal)),
    output_metrics=Path(str(snakemake.output.metrics)),
    output_provenance=Path(str(snakemake.output.provenance)),
)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).write_text(
    f"Built summit-window H3K27ac signal for {metrics['row_count']} DHS-context rows\n",
    encoding="utf-8",
)
