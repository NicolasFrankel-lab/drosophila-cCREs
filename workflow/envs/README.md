# Rule environments

These files define isolated tool groups for Snakemake rules. With the local
profile, Snakemake creates resolved environments under `.snakemake/conda` in the
repository root. Do not create them as named environments in Mamba's global
environment directory.

- `read_qc.yaml`: raw/trimmed FASTQ QC and adapter trimming
- `alignment.yaml`: short-read alignment and BAM filtering
- `peaks.yaml`: MACS3 callpeak and bdgcmp/qpois
- `atac_qc.yaml`: Tn5 insertion preparation, qpois refinement, BigWigs, and ATAC QC
- `chip_qc.yaml`: R-based H3K27ac cross-correlation QC only
- `reporting.yaml`: MultiQC aggregation and standard-library Python catalog, nearest-TSS and validation scripts
- `activity_tmm.yaml`: edgeR TMM factors for 10-kb autosomal background counts
- `catalog_tracks.yaml`: bedtools/coreutils/pyBigWig catalog signal-track construction
- `final_report.yaml`: deterministic HTML-to-PDF rendering for integrated QC
- `footprinting.yaml`: samtools and TOBIAS for context-pooled ATAC footprinting
