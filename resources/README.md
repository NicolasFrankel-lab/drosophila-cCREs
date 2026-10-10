# Pipeline resources

This directory holds the versioned inputs of the eight-context dm6 atlas and
the adapter sequences used for trimming.

## Atlas inputs

- `atlas_samples.tsv`: the pipeline input. Four columns (`accession`,
  `library_id`, `assay`, `context`) and 35 public accessions forming 34
  biological libraries: 20 ATAC-seq and 14 IP-only H3K27ac libraries across the
  contexts `ab`, `e13`, `e5`, `ead`, `hid`, `lb`, `o` and `wid`. `e5_atac_rep1`
  has two technical runs. `o_h3k27ac_cutrun_rep1` is CUT&RUN; the other
  H3K27ac libraries are ChIP-seq.
- `atlas_qc_decisions.tsv`: one QC decision per biological library
  (`library_id`, `assay`, `context`, `qc_decision`,
  `estimated_fragment_length_bp`, `notes`). All 34 libraries pass. Single-end
  H3K27ac rows carry the reviewed fragment length used to extend reads; it is
  blank for paired-end and ATAC libraries. Apply it to each per-assay QC
  manifest with `src/review_final_bam_manifest.py --assay`.

Detailed dataset selection and QC evidence belong to the separate
atlas-analysis repository. The workflow writes its generated manifests below
`results/<project>/<run_id>/provenance/manifests/`; reviewed manifests derived
from the decision table belong below the ignored `data/reviewed/` namespace.

## Adapter sequences

The 158 sequence records in `adapters/adapters.fa` come from
`BBMap_39.96.tar.gz:bbmap/resources/adapters.fa` in the parent workspace.

- BBMap archive SHA-256: `e173bdd0d3ca047f378c71dad568a148596c1690bf36abca93e918569c9fb382`
- Upstream extracted FASTA SHA-256: `85abe9d3e40dc37c968f7e4c1227e05976a4ed0583d1dd442d375aa7516f13a9`
- Repository FASTA SHA-256: `74e19a3b2b09f8fa84bf2e59877f025591bfe64fb50fcb8ee6169631ccb58468`
- Records: 158

The complete file contains TruSeq, Nextera, and other sequences. It must not be
passed wholesale to Cutadapt as though all 158 entries were equivalent 3-prime
adapters.

The reviewed `adapters/nextera.fa` and `adapters/truseq.fa` subsets contain the
read-through sequences associated with the named workflow presets. The full
BBMap collection is retained for provenance and custom adapter review.
