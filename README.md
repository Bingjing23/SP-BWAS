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
  bwasSampleSize = c(NAD, "NMISS"),
  refPanel = c("AVERAGE"),
  outputPath = "outputs/pilot/"
)
```

For the first pilot run:

- AD-side input: `AlzDisease_LMM/BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat`
- Trait-side input: `SSTAT_BingJing/sstat_FS_All_moda_total_hyperTension.linear`
- AD sample size: replace `NAD` with the confirmed AD vs HC BWAS sample size.
- Trait sample size: use the `NMISS` column from the UKB trait file.
- Reference panel: `AVERAGE`.
- Output directory: `outputs/pilot/`.
- Log directory: `logs/`.

Raw input files should not be modified. UKB trait files use `Voxel` as the
spatial identifier, while brainMapR expects `Probe`. Create derived corrected
copies under `outputs/pilot/derived_inputs/` and rename only the header column
`Voxel` to `Probe`.

Server environment used for brainMapR:

```bash
module purge
module load R/4.5.0

export R_LIBS_USER=/mnt/backedup/home/bingjinZ/Rlibs/R-4.5.0
export RGL_USE_NULL=TRUE
export LD_LIBRARY_PATH=/software/conda-envs/pkgs/bzip2-1.0.8-h4bc722e_7/lib:$LD_LIBRARY_PATH
```

`brainMapR` installation:

```r
install.packages("devtools", lib = Sys.getenv("R_LIBS_USER"),
                 repos = "https://cloud.r-project.org")
devtools::install_github("baptisteCD/brainMapR", force = TRUE,
                         upgrade = "never")
library(brainMapR)
```

Before running the pilot, verify:

```bash
Rscript -e 'library(brainMapR); sessionInfo()'
Rscript -e 'args(brainMapR::sumR2_regression_bivariate)'
```

## Main outputs

- Phenotype tables such as `df_*`
- BWAS result files under `Results/`
- AD summary statistics in `AlzDisease_LMM/`
