# AGENTS.md

Repository-specific guidance for coding agents and maintainers.

## Priorities

1. Scientific and computational correctness
2. Reproducibility and idempotency
3. Minimal focused changes
4. Verification
5. Clear user-facing behavior
6. Security and data hygiene

Do not silently change scientific defaults, thresholds, output formats, file
names, or sample relationships.

## Repository contract

- The public data input is an accession CSV/TSV conforming to
  `schemas/sample-sheet.schema.yaml`.
- The canonical user-facing columns are `accession`, `library_id`, `assay`, and
  `context`. Genome is a project-level CLI value; assay defaults and accession
  layout are resolved rather than repeated per row.
- Users do not provide raw FASTQ paths to the primary entry point.
- `src/run_pipeline.py` is the canonical accession-to-results command.
- `workflow/Snakefile` is the canonical processing workflow. Put shared rules
  under `workflow/rules/`; avoid creating parallel experimental Snakefiles for
  production behavior.
- `dm6` and `hg38` reference preparation belongs in the workflow DAG.
- Raw FASTQs are immutable inputs. Never overwrite or modify them in place.
- Technical runs share a `library_id`, run lane-level QC/alignment separately,
  and merge before biological-library duplicate marking.
- Biological replicates have distinct `library_id` values and remain separate
  through replicate peak calling. ATAC libraries with the same `context` are
  pooled automatically.

## Scientific invariants

- ATAC uses the two-ended Tn5/MACS3-qpois branch.
- Paired-end ATAC qpois retains proper pairs with `0 < abs(TLEN) < 150`, applies
  Tn5 offsets, and converts both mates to one-base insertion records. MACS3 uses
  `-f BED --nomodel --shift -75 --extsize 150 --keep-dup all -B -q 0.10`.
- ATAC qpois is computed with `macs3 bdgcmp -m qpois` from the unscaled
  treatment pileup and local lambda. Do not add `--SPMR` to this branch.
- Qpois refinement evaluates exponents 2-325 and retains 50-400 bp components;
  broad components split as the exponent rises. Refined boundaries are not an
  independent FDR-controlled peak call.
- The sample-sheet `context` column controls ATAC replicate pooling. Replicate
  calls are evidence; a context peak must be covered by the configured fraction
  in the configured number of biological libraries.
- The final ATAC endpoint is a summit-aware master DHS registry across the
  replicate-supported contexts. Summits come from pooled context signal;
  clustering requires reciprocal summit containment, a complete-linkage summit
  span of at most 150 bp by default, and at most one peak per context. Do not
  replace this with a simple interval union or fixed-width resizing. Adjacent
  clusters with representative summits less than 50 bp apart merge only when
  their context sets are disjoint; a shared context is evidence for two sites.
- H3K27ac integration through context-resolved activity annotation is part of
  this processing DAG. H3K27ac activity uses the maximum background-TMM
  normalized CPM/kb among the three 500-bp summit windows. Its context-specific
  background is blacklist-clear master DHSs that do not contribute to that
  context, and activity requires a blacklist-clear context member with a
  midrank background percentile strictly greater than 0.60. The default links branch is contact-independent and
  assigns each proximal/distal enhancer to its nearest annotated TSS, retaining
  exact ties in one enhancer row and annotating promoter activity by context.
  For the canonical nine-context dm6 atlas, contact-map normalization and
  context-specific element--promoter candidate links are an explicit optional
  branch selected with `--with-contacts`; ABC candidate resizing/scoring
  remains outside it.
- Contact replicates are coarsened only by exact integer count summation and
  merged as raw counts before ICE balancing. Balancing starts with a
  200-iteration limit and restarts from uniform weights with a 2,000-iteration
  limit only if the first attempt does not converge. The seven observed
  contexts retain assay/match caveats; `e13` and `hid` use an explicitly
  labeled atlas-wide distance model rather than a substituted contact map.
- Contact links connect context-member master DHSs to distinct GTF gene/TSS
  promoters within 1 Mb. Candidate scores combine contact weight with promoter
  ATAC/H3K27ac activity and are rankings, not calibrated probabilities or
  causal target assignments. Promoter activity is supported by context-member
  DHSs whose summits lie within an inclusive 500 bp of the exact TSS; it does
  not use DHS-boundary overlap with the promoter BED interval. The focused
  enhancer--gene table retains only
  proximal/distal enhancer-like elements with the catalog H3K27ac activity call
  and observed promoter edges with observed/expected contact at least 1 before
  gene collapse; its distance-modeled contexts are header-only. Distance-only
  contexts additionally receive a separate active-promoter distance projection
  whose observed contact and observed/expected fields remain blank.
  It is derived downstream from completed node/edge tables and must not trigger
  contact normalization or exhaustive link reconstruction when backfilled.
- TF ChIP defaults to narrow MACS3 peaks; histone ChIP defaults to broad peaks.
- ChIP controls are explicit `control_library` relationships. IP-only ChIP is
  represented by omitting that optional column and runs MACS3 without `-c`.
- ChIP MACS3 `callpeak` outputs include `-B --SPMR` treatment-pileup and
  control-lambda bedGraphs. The ATAC qpois branch deliberately omits `--SPMR`.
- FRiP denominators must be documented and fragment-aware for paired-end data.

If a task changes any invariant, update the schema, validators, tests, README,
and provenance together and call out the behavior change explicitly.

## Environment and dependencies

- Keep the orchestration environment at repository-root `.venv`.
- Set `MAMBA_ROOT_PREFIX="$PWD/.micromamba"` when practical.
- Set `CONDA_PKGS_DIRS="$PWD/.micromamba/pkgs"` and
  `CONDA_ENVS_PATH="$PWD/.micromamba/envs"` so package caches and named-prefix
  resolution remain inside the repository.
- Set `CONDARC="$PWD/.condarc"` to use the repository's strict channel policy.
- Keep Snakemake rule environments below `.snakemake/conda` using the local
  profile.
- Add orchestration dependencies to `environment.yml` only when required by
  repository commands or tests.
- Add bioinformatics dependencies to the smallest applicable file under
  `workflow/envs/`.
- R/Bioconductor dependencies belong only in `workflow/envs/chip_qc.yaml` unless
  a new R-dependent rule is explicitly introduced.

Do not create named environments outside the repository.

## Idempotency and file safety

- Every command must be safe to rerun with identical inputs and parameters.
- Prefer temporary/staging outputs followed by atomic replacement.
- Do not replace unchanged manifests or resolved configs merely to update a
  timestamp.
- Preserve aria2 resume state and validate reported checksums.
- Promote SRA-derived FASTQs only after conversion and compression complete.
- Let Snakemake own output completeness and keep `rerun-incomplete: true` in
  execution profiles.
- Never write a new parameter attempt into an existing scientific run
  namespace; use a new `run_id`.

## Parallelism and resources

- Preserve independent lane and sample jobs so Snakemake can schedule them
  concurrently.
- Declare realistic `threads` and `mem_mb` on every compute-heavy rule.
- Avoid nested oversubscription: one rule's subprocesses share its declared
  thread allocation.
- Make aggregate concurrency a CLI/profile concern (`--cores`, `--jobs`, and
  resource limits), not a hard-coded workflow limit.

## SLURM

- All site-specific SLURM scripts and profiles belong under ignored `slurm/`.
- Do not commit cluster hostnames, accounts, partitions, credentials, or home
  directory paths.
- Never run download, alignment, peak calling, or environment installation on
  the `cranex` login node. Submit compute work through SLURM.
- Long cluster runs must write logs under `logs/` and be safe to resubmit.

## Data and Git hygiene

Do not commit:

- `.venv/`, `.micromamba/`, `.snakemake/`, caches, or locks;
- raw FASTQs, references, BAM/CRAM files, BigWigs, bedGraphs, or generated
  results;
- SLURM scripts or site profiles;
- secrets, credentials, private metadata, or large temporary fixtures.

The versioned `resources/` directory is intentionally limited to adapter FASTAs,
the two canonical ATAC/H3K27ac atlas input tables, and the reviewed contact
source manifest. Detailed selection/QC metadata belongs in the atlas-analysis
repository. Do not add random download selections or generated smoke-test data
here.

Preserve unrelated user changes in a dirty worktree. Use `rg`/`rg --files` for
search and `apply_patch` for focused file edits.

## Verification

Run the smallest relevant checks during development, then the full practical
suite before handoff:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export CONDA_PKGS_DIRS="$PWD/.micromamba/pkgs"
export CONDA_ENVS_PATH="$PWD/.micromamba/envs"
export CONDARC="$PWD/.condarc"
export XDG_CACHE_HOME="$PWD/.cache"

mamba run --prefix "$PWD/.venv" pytest -q

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --cores 8 --dry-run

git diff --check
```

The loopback aria2 integration test may require network/socket permission in a
sandbox. Report skipped tests and environmental failures honestly.

For schema/configuration changes, include positive and negative validation
tests. For workflow changes, include a Snakemake dry-run assertion for every
new rule branch. For qpois refinement changes, test interval boundaries,
exponent semantics, empty inputs, and deterministic output ordering.

## Documentation

Update `README.md` whenever setup, inputs, commands, defaults, dependencies,
outputs, or restart behavior changes. Update `resources/README.md` when atlas
or adapter provenance changes. Regenerate `docs/workflow-dag.svg` after
changing rule dependencies.

Final handoffs must state what changed, what was tested, what passed or failed,
what was not tested, and any remaining scientific uncertainty.
