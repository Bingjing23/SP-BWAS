# SP-BWAS

This repository contains the code and analysis notes for an AD-related brain-wide association study (BWAS) workflow.

## What is here

- `SSTAT_BingJing/16_BingJing.Rmd`: R Markdown script for building UK Biobank-derived phenotypes and running BWAS/BLUP-style analyses with `OSCA`.
- `AlzDisease_LMM/`: AD BWAS result files and related outputs.

## Notes

- Data files are not included in version control.
- The scripts assume the original UK Biobank and REML/LDSC files already exist on the same server paths used in the Rmd.
- `brainMapR` is used downstream for summary-statistic based analyses.

## Quick start

1. Open [`SSTAT_BingJing/16_BingJing.Rmd`](SSTAT_BingJing/16_BingJing.Rmd).
2. Update the file paths if needed for your environment.
3. Run the phenotype construction and BWAS steps on the compute environment that has `OSCA` and the required UK Biobank inputs.

## Downstream brainMapR analysis

The intended downstream analysis compares Alzheimer BWAS summary statistics in
`AlzDisease_LMM/` against UKB-derived risk-factor BWAS summary statistics in
`SSTAT_BingJing/` using:

```r
brainMapR::sumR2_regression_bivariate(
  inputPath = c("AlzDisease_LMM/", "outputs/pilot/derived_inputs/"),
  bwasFile = c(
    "BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat",
    "sstat_FS_All_moda_total_hyperTension.Probe.linear"
  ),
  bwasSampleSize = c(3542, "NMISS"),
  refPanel = c("AVERAGE"),
  outputPath = "outputs/pilot/brainMapR_ADvsHC_hyperTension/"
)
```

The official pilot using the inputs below has run successfully:

- AD-side input: `AlzDisease_LMM/BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat`
- Trait-side input: `outputs/pilot/derived_inputs/sstat_FS_All_moda_total_hyperTension.Probe.linear`
- AD sample size: fixed numeric value `3542`.
- Trait sample size: use the `NMISS` column from the UKB trait file.
- Reference panel: `AVERAGE`.
- Output directory: `outputs/pilot/brainMapR_ADvsHC_hyperTension/`.
- Log directory: `logs/`.

Raw input files should not be modified. UKB trait files use `Voxel` as the
spatial identifier, while brainMapR expects `Probe`. Create derived corrected
copies under `outputs/pilot/derived_inputs/` and rename only the header column
`Voxel` to `Probe`.

AD meta-analysis files in `AlzDisease_LMM/` do not carry `NMISS` because the
sample size is not preserved in the meta-analysis summary statistics. Supply
fixed numeric sample sizes for the AD-side input, and keep using `"NMISS"` for
UKB BWAS files that contain an `NMISS` column. The currently confirmed AD
sample sizes are:

| AD trait | Sample size |
| --- | ---: |
| ADvsHC | 3542 |
| MCIvsHC | 3976 |
| Conversion1year | 1257 |
| Conversion2years | 1199 |
| Conversion3years | 1031 |
| Conversion4years | 1285 |
| Conversion5years | 1197 |
| MMSE | 6981 |

Server environment used for brainMapR:

```bash
module purge
module load R/4.5.0

export R_LIBS_USER=/mnt/backedup/home/bingjinZ/Rlibs/R-4.5.0
export RGL_USE_NULL=TRUE
export LD_LIBRARY_PATH=/software/conda-envs/pkgs/bzip2-1.0.8-h4bc722e_7/lib:${LD_LIBRARY_PATH:-}
```

Use the narrow `bzip2` library path above for `rgl` on the QIMR HPC. Do not use
the broader `/software/conda-envs/lib` path because it can cause
`magick`/ImageMagick dynamic-library conflicts.

Official package route:

```r
install.packages("devtools", lib = Sys.getenv("R_LIBS_USER"),
                 repos = "https://cloud.r-project.org")
devtools::install_github("jean997/GFA", force = TRUE)
devtools::install_github("baptisteCD/brainMapR", force = TRUE)
library(brainMapR)
library(GFA)
```

The validated package versions were:

- `brainMapR_1.1.0.9000`
- `GFA_1.0.0.0449`

Do not use the CRAN `GFA_1.0.5`, the older RDM
`brainMapR_0.8.0.9000.tar.gz`, or a local replacement for `GFA::ldsc_rg()` for
the official analysis.

Before running pilot or batch jobs, verify:

```bash
Rscript -e 'library(brainMapR); library(GFA); stopifnot("sumR2_regression_bivariate" %in% ls("package:brainMapR")); stopifnot("ldsc_rg" %in% ls("package:GFA")); stopifnot("snp_ldsc" %in% ls("package:GFA")); sessionInfo()'
```

Submit long brainMapR jobs through PBS on a compute node rather than running
long `Rscript` jobs on the login node.

## Analysis plan

The analysis proceeds in three stages: pilot, small batch, then full batch.
This keeps the workflow reproducible while avoiding a large failed run caused
by a simple input-format or environment issue.

The scientific comparison is map-level rather than SNP-level: each job compares
one AD-related BWAS map from `AlzDisease_LMM/` with one UKB risk-factor BWAS map
from `SSTAT_BingJing/` using
`brainMapR::sumR2_regression_bivariate()`. The target output is a pairwise
comparison table or matrix whose rows are AD-related BWAS maps and whose
columns are UKB risk-factor BWAS maps.

Current design decisions:

- Use the official `brainMapR_1.1.0.9000` plus `jean997/GFA` route only.
- Use `AVERAGE` as the primary reference panel for pilot and small batch.
  `UKB` can be added later as a sensitivity reference panel after the workflow
  is stable.
- Preserve all raw input files.
- Use fixed numeric sample sizes for AD meta-analysis BWAS files. The currently
  included AD sample sizes were confirmed by Baptiste because these
  meta-analysis summary statistics do not carry `NMISS`.
- Use the `NMISS` column for UKB risk-factor BWAS files.
- Convert UKB `Voxel` headers to `Probe` only in derived copies under
  `outputs/batch/derived_inputs/`.
- Some UKB files were written with an extra row-name column before `Chr`; the
  derived-input script removes that extra column when detected. Header-only
  UKB files are excluded from regenerated batch manifests.
- Include only AD traits with confirmed sample sizes in default batch designs.
- Exclude obvious non-analysis UKB variables such as `eid`, `sex`, and pilot
  duplicate files from the default risk-factor batch.

The current small batch is designed as a controlled scale-up from the successful
pilot:

```text
AD maps:
  ADvsHC
  MCIvsHC

Risk-factor maps:
  hyperTension
  T2D
  hyperCholesterolemia
  ldl_direct
  majorDepression
  smoking_ever
```

This gives 12 pairwise jobs. These traits were chosen because they include the
pilot trait plus vascular/metabolic, psychiatric, and lifestyle candidates that
are directly relevant to AD risk. This follows the meeting strategy of testing
a small number of trait pairs first, then rolling the same workflow out to all
selected traits. If the small batch completes cleanly, the same
manifest-driven workflow can be extended to the full default design:

```text
8 AD maps with confirmed sample sizes x 68 UKB risk-factor maps = 544 jobs
```

The full grid is a default engineering design and can be adjusted if a smaller
core trait set is preferred.

Do not interpret small-batch results as final biological conclusions. The small
batch is primarily a workflow validation and early signal-checking step. Full
interpretation should wait until the full batch is complete and outputs have
been collected into comparable matrices, including gray-matter correlation and
sumR2-style outputs where available.

## Batch workflow

The batch workflow is manifest-driven. Do not hand-code pairwise comparisons
one by one.

From the project root on the HPC:

```bash
module purge
module load R/4.5.0

export R_LIBS_USER=/mnt/backedup/home/bingjinZ/Rlibs/R-4.5.0
export RGL_USE_NULL=TRUE
export LD_LIBRARY_PATH=/software/conda-envs/pkgs/bzip2-1.0.8-h4bc722e_7/lib:${LD_LIBRARY_PATH:-}
```

First audit all AD and risk-factor BWAS inputs:

```bash
Rscript scripts/00_check_inputs.R
```

This writes:

- `outputs/input_audit/ad_bwas_inventory.tsv`
- `outputs/input_audit/risk_factor_bwas_inventory.tsv`
- `outputs/input_audit/input_readiness.tsv`
- `outputs/input_audit/blockers.tsv`

Then build the manifests and pairwise design tables:

```bash
Rscript scripts/01_make_manifests.R
```

This writes:

- `manifests/ad_bwas_manifest.tsv`
- `manifests/risk_factor_bwas_manifest.tsv`
- `manifests/brainmapr_small_batch_design.tsv`
- `manifests/brainmapr_pairwise_design.tsv`

The current default design includes the 8 AD traits with confirmed sample
sizes and excludes AD traits whose sample sizes are not yet confirmed. It also
excludes obvious non-analysis UKB columns such as `eid`, `sex`, and pilot
duplicate files from the default risk-factor batch.

Before running brainMapR, create derived UKB inputs with `Probe` headers. For
the recommended small batch only:

```bash
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_small_batch_design.tsv \
  --force
```

For the full default batch:

```bash
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_pairwise_design.tsv \
  --force
```

Run a dry-run path check before submitting real jobs:

```bash
Rscript scripts/04_run_brainMapR_batch.R \
  --design manifests/brainmapr_small_batch_design.tsv \
  --dry-run
```

Submit the small batch through PBS:

```bash
qsub -v DESIGN=manifests/brainmapr_small_batch_design.tsv \
  scripts/run_brainMapR_batch.pbs
```

To run a single pair:

```bash
qsub -v DESIGN=manifests/brainmapr_small_batch_design.tsv,PAIR_ID=ADvsHC__hyperTension \
  scripts/run_brainMapR_batch.pbs
```

After jobs finish, collect outputs:

```bash
Rscript scripts/05_collect_brainMapR_outputs.R \
  --design manifests/brainmapr_small_batch_design.tsv \
  --output-dir outputs/batch/summary_small_batch
```

For full batch collection, use `manifests/brainmapr_pairwise_design.tsv` and
write to `outputs/batch/summary/`.

## Parallel full batch

After the small batch succeeds, run the full default grid as a PBS array rather
than one serial job. Each array task runs one row of
`manifests/brainmapr_pairwise_design.tsv`, so the scheduler can execute many
AD x risk-factor pairs in parallel.

Prepare all full-batch derived inputs first:

```bash
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_pairwise_design.tsv \
  --force
```

Confirm the number of jobs:

```bash
full_n=$(($(wc -l < manifests/brainmapr_pairwise_design.tsv) - 1))
echo "$full_n"
```

The current expected value is `544`.

Submit the full batch as an array:

```bash
qsub -J 1-${full_n} \
  -v DESIGN=manifests/brainmapr_pairwise_design.tsv \
  scripts/run_brainMapR_batch_array.pbs
```

If the cluster limits the number of simultaneous array tasks, use the local
PBS syntax for throttling, for example `-J 1-${full_n}%40` if supported. If
array throttling is not supported, submit smaller ranges manually:

```bash
qsub -J 1-100 -v DESIGN=manifests/brainmapr_pairwise_design.tsv scripts/run_brainMapR_batch_array.pbs
qsub -J 101-200 -v DESIGN=manifests/brainmapr_pairwise_design.tsv scripts/run_brainMapR_batch_array.pbs
```

Monitor:

```bash
qstat -u $USER
ls -lh logs/brainMapR_array.*.R.log | tail
```

After the array finishes, check pair statuses:

```bash
find outputs/batch/brainMapR_pairs -name run_status.tsv \
  -exec awk -F '\t' 'NR > 1 {print $4}' {} \; | sort | uniq -c
```

Any failed pair writes its own status file under
`outputs/batch/brainMapR_pairs/<pair_id>/run_status.tsv` and a compact failure
record under `outputs/batch/failures/`.

Collect full-batch outputs:

```bash
Rscript scripts/05_collect_brainMapR_outputs.R \
  --design manifests/brainmapr_pairwise_design.tsv \
  --output-dir outputs/batch/summary
```

## Main outputs

- Phenotype tables such as `df_*`
- BWAS result files under `Results/`
- AD summary statistics in `AlzDisease_LMM/`
