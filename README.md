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
  --design manifests/brainmapr_small_batch_design.tsv
```

For the full default batch:

```bash
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_pairwise_design.tsv
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

## Main outputs

- Phenotype tables such as `df_*`
- BWAS result files under `Results/`
- AD summary statistics in `AlzDisease_LMM/`
