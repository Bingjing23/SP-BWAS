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

## Main outputs

- Phenotype tables such as `df_*`
- BWAS result files under `Results/`
- AD summary statistics in `AlzDisease_LMM/`

