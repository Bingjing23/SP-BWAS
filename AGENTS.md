# AGENTS.md

## Project scope

This repository is an Alzheimer / BWAS / brainMapR downstream integration project.

The project should be interpreted using the existing repository structure:

```text
SP-BWAS/
├── AlzDisease_LMM/
├── SSTAT_BingJing/
└── README.md
```

Do not reorganise this repository unless explicitly instructed. The current folder layout is intentional and should be preserved.

This project is **not** an fMRI preprocessing project, **not** a deep learning embedding project, **not** a PRS project, and **not** a SNP-level GWAS project. Do not introduce or modify workflows related to fMRI alignment, foundation models, embeddings, disease prediction, PRS, or neural-network architectures unless explicitly requested.

The immediate goal is to determine whether the current files are sufficient to run downstream map-level integration using:

```r
brainMapR::sumR2_regression_bivariate()
```

## Scientific objective

The scientific goal is to compare brain-wide association maps for Alzheimer-related phenotypes against brain-wide association maps for Alzheimer-relevant risk factors.

The main question is:

> Which Alzheimer-relevant risk-factor BWAS maps share spatial signal with AD-related BWAS maps, including AD vs HC, MCI vs HC, conversion phenotypes, and cognitive/clinical measures?

The expected final output is a pairwise comparison matrix:

```text
rows    = Alzheimer-related BWAS outcomes from AlzDisease_LMM/
columns = risk-factor / trait BWAS maps from SSTAT_BingJing/
cells   = brainMapR / sumR2 outputs, such as shared signal, sum R2, correlation, SE, p-value, or CI
```

## Directory roles

### AlzDisease_LMM/

This directory contains Alzheimer-related BWAS summary statistics.

These files should be treated as the AD-related outcome maps. Examples include:

- `BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_AD_vs._MCI_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_MCI_vs._HC_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Conversion_at_1_year_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Conversion_at_2_years_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Conversion_at_3_years_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Conversion_at_4_years_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Conversion_at_5_years_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Alzheimer's_progression_-__3_stages_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Mini_Mental_State_Examination_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Clinical_Dementia_Rating_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Clinical_Dementia_Rating_-_sum_of_boxes_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_RAVLT_*.linear_random_oscaFormat`
- `BWAS_meta_Geriatric_Depression_Scale_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Functional_assessment_questionnaire_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Logical_memory_*.linear_random_oscaFormat`
- `BWAS_meta_Maternal_AD_QC_8SD.linear_random_oscaFormat`
- `BWAS_meta_Paternal_AD_QC_8SD.linear_random_oscaFormat`

Do not assume biological interpretation from filenames alone. Use filenames only to infer candidate phenotype labels for manifest construction.

### SSTAT_BingJing/

This directory contains UKB-derived BWAS summary statistics for Alzheimer-relevant risk factors and traits.

These files should be treated as candidate risk-factor / external trait maps. Examples include:

Vascular / metabolic:

- `sstat_FS_All_moda_total_hyperTension.linear`
- `sstat_FS_All_moda_total_htn_i10_preimg.linear`
- `sstat_FS_All_moda_total_T2D.linear`
- `sstat_FS_All_moda_total_t2d_e11_preimg.linear`
- `sstat_FS_All_moda_total_diabetes_doctor_dx.linear`
- `sstat_FS_All_moda_total_diabetes_preimg.linear`
- `sstat_FS_All_moda_total_hyperCholesterolemia.linear`
- `sstat_FS_All_moda_total_ldl_direct.linear`
- `sstat_FS_All_moda_total_bmi.linear`
- `sstat_FS_All_moda_total_stroke.linear`

Psychiatric:

- `sstat_FS_All_moda_total_majorDepression.linear`
- `sstat_FS_All_moda_total_mdd_broad_preimg.linear`
- `sstat_FS_All_moda_total_anxDis.linear`
- `sstat_FS_All_moda_total_anxiety_broad_preimg.linear`
- `sstat_FS_All_moda_total_PTSD.linear`
- `sstat_FS_All_moda_total_ptsd_broad_preimg.linear`
- `sstat_FS_All_moda_total_BD.linear`
- `sstat_FS_All_moda_total_bd_broad_preimg.linear`
- `sstat_FS_All_moda_total_SCZ.linear`
- `sstat_FS_All_moda_total_scz_broad_preimg.linear`
- `sstat_FS_All_moda_total_phobAnx.linear`

Lifestyle:

- `sstat_FS_All_moda_total_smoking_ever.linear`
- `sstat_FS_All_moda_total_pack_years.linear`
- `sstat_FS_All_moda_total_pack_years_adult_prop.linear`
- `sstat_FS_All_moda_total_sleep_duration.linear`
- `sstat_FS_All_moda_total_insomnia.linear`
- `sstat_FS_All_moda_total_daytime_napping.linear`
- `sstat_FS_All_moda_total_Alcohol_merge_phenotype.linear`
- `sstat_FS_All_moda_total_Alcohol_week_avg.linear`
- `sstat_FS_All_moda_total_alcohol_intake_frequency.linear`
- `sstat_FS_All_moda_total_ever_drinker.linear`
- `sstat_FS_All_moda_total_*_weekly_intake.linear`
- `sstat_FS_All_moda_total_pa_*.linear`

Social / socioeconomic:

- `sstat_FS_All_moda_total_age_completed_education.linear`
- `sstat_FS_All_moda_total_income_before_tax.linear`
- `sstat_FS_All_moda_total_imd_england.linear`
- `sstat_FS_All_moda_total_leisure_social_score.linear`
- `sstat_FS_All_moda_total_loneliness_isolation.linear`

Frailty / multimorbidity:

- `sstat_FS_All_moda_total_frailty_fi49.linear`
- `sstat_FS_All_moda_total_frailty_missing_count.linear`
- `sstat_FS_All_moda_total_multimorbidity.linear`
- `sstat_FS_All_moda_total_med_20003_n_distinct_i2.linear`

Hormonal / sex-specific:

- `sstat_FS_All_moda_total_Menopause.linear`
- `sstat_FS_All_moda_total_hrt_ever_used.linear`

Infection / inflammatory / other candidate risk factors:

- `sstat_FS_All_moda_total_infect_any_ever.linear`
- `sstat_FS_All_moda_total_infect_count_5y.linear`
- `sstat_FS_All_moda_total_infect_count_total.linear`
- `sstat_FS_All_moda_total_infect_recency_days.linear`
- `sstat_FS_All_moda_total_lrti_ever.linear`
- `sstat_FS_All_moda_total_pneumonia_ever.linear`
- `sstat_FS_All_moda_total_sepsis_ever.linear`
- `sstat_FS_All_moda_total_ssti_ever.linear`
- `sstat_FS_All_moda_total_uti_ever.linear`
- `sstat_FS_All_moda_total_periodontal_k05_preimg.linear`

The file `16_BingJing.Rmd` is the upstream phenotype construction / OSCA-running notebook. Do not rewrite it unless explicitly asked. It can be inspected to understand trait definitions.

## Do not change project structure

The user does not want to change the project structure.

Do not move existing files.

Do not create a new `data/`, `metadata/`, `scripts/`, or `output/` hierarchy unless explicitly instructed.

If new files are needed, prefer adding them at the project root or in minimally invasive folders such as:

```text
manifests/
scripts/
outputs/
logs/
```

Only create these if needed and after reporting the plan.

## Immediate audit task

Before making code changes, inspect and report:

1. Which AD-related BWAS files exist in `AlzDisease_LMM/`.
2. Which risk-factor BWAS files exist in `SSTAT_BingJing/`.
3. Whether each file is readable.
4. Whether each file has a `Probe` column.
5. If `Probe` is absent, whether it has a `Voxel` column that needs renaming to `Probe`.
6. Whether each file has `NMISS` or another sample-size column.
7. Whether effect-size and uncertainty columns exist, such as `BETA`, `SE`, `P`, or equivalent.
8. Whether files appear ready for `brainMapR`.
9. Whether reference panel files are present or expected externally.
10. Whether existing scripts already run `brainMapR::sumR2_regression_bivariate()`.

Do not edit files during the first audit unless explicitly instructed.

## Required audit tables

When auditing the repository, report these tables.

### AD-related BWAS inventory

```text
phenotype | file_path | file_name | format | first_columns | sample_size_source | ready_for_brainMapR | notes
```

### Risk-factor BWAS inventory

```text
trait | category | file_path | file_name | format | first_columns | sample_size_source | ready_for_brainMapR | notes
```

### Input readiness table

```text
file | has_Probe | has_Voxel | needs_Voxel_to_Probe | has_NMISS | has_beta | has_se | has_p | readable | ready | issue
```

### Blockers table

```text
blocker | affected_files | severity | proposed_fix
```

## Input-format rules

The most important input-format rule is:

```text
Voxel -> Probe
```

The brainMapR workflow expects a `Probe` identifier. If a file has `Voxel` instead of `Probe`, do not modify the raw file in place. Create a corrected derived file only when explicitly instructed.

Candidate BWAS summary statistics files should be checked for:

1. `Probe` column.
2. `Voxel` column if `Probe` is missing.
3. `NMISS` or another usable sample-size column.
4. Effect-size column, such as `BETA` or equivalent.
5. Standard error column, such as `SE` or equivalent.
6. P-value column, such as `P` or equivalent.
7. Tabular format readable by R.
8. Consistent identifiers across AD and risk-factor files.
9. No malformed headers.

## Reference panel assumptions

Assume the main analysis should use the `AVERAGE` reference panel unless code or documentation clearly says otherwise.

Assume `UKB` can be used as a sensitivity-analysis reference panel.

If reference panel files are not present in this repository, report them as missing or external dependencies. Do not invent paths.

## Minimal pilot analysis

The first pilot should be:

```text
AD-related file:
AlzDisease_LMM/BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat
```

paired with one strong AD risk-factor BWAS file from `SSTAT_BingJing/`, preferably one of:

```text
SSTAT_BingJing/sstat_FS_All_moda_total_hyperTension.linear
SSTAT_BingJing/sstat_FS_All_moda_total_T2D.linear
SSTAT_BingJing/sstat_FS_All_moda_total_hyperCholesterolemia.linear
SSTAT_BingJing/sstat_FS_All_moda_total_ldl_direct.linear
SSTAT_BingJing/sstat_FS_All_moda_total_majorDepression.linear
SSTAT_BingJing/sstat_FS_All_moda_total_smoking_ever.linear
```

The pilot should verify:

1. Both files are readable.
2. Required columns are present.
3. `Voxel -> Probe` is handled if required.
4. Sample size is correctly supplied using `NMISS` or a fixed N.
5. Reference panel can be loaded.
6. `brainMapR::sumR2_regression_bivariate()` runs successfully.
7. Output files are written to a clearly named output directory.

## Batch analysis design

After the pilot succeeds, batch analysis should be manifest-driven.

Use a design table conceptually like:

```text
ad_id | ad_file | ad_sample_size | trait_id | trait_file | trait_category | trait_sample_size | reference_panel | output_dir
```

Do not manually hard-code dozens of pairwise runs if a design table can be used.

The intended batch comparison is:

```text
AlzDisease_LMM/*.linear_random_oscaFormat
    x
SSTAT_BingJing/*.linear
```

but only after input readiness is verified.

## Script expectations

If scripts are added, prefer minimal and explicit R scripts.

Suggested script names, if needed:

```text
scripts/00_check_inputs.R
scripts/01_make_manifests.R
scripts/02_fix_sumstats_headers.R
scripts/03_run_brainMapR_pilot.R
scripts/04_run_brainMapR_batch.R
scripts/05_collect_brainMapR_outputs.R
```

Scripts should:

1. Preserve raw input files.
2. Fail early when required files or columns are missing.
3. Print all input and output paths.
4. Avoid hidden assumptions.
5. Avoid hard-coded absolute paths unless the user explicitly requests local execution.
6. Separate pilot output from batch output.
7. Save logs sufficient for debugging.
8. Record package versions where practical.

## Biological interpretation constraints

Do not over-interpret map-level results.

Use cautious terms:

- shared spatial signal
- map-level similarity
- brain-wide association pattern
- cross-map comparison
- AD-related BWAS map
- risk-factor-associated BWAS map

Avoid unsupported causal claims:

- do not say a risk factor causes AD through a specific brain region
- do not call these SNP-level GWAS results
- do not interpret file names as definitive phenotype definitions without checking code or documentation
- do not conflate this project with fMRI embedding or disease prediction work

## Reporting style

When reporting back to the user, be concise and actionable.

Always state:

1. What was found.
2. What is missing.
3. Whether the project is currently executable.
4. What exact pilot should be run first.
5. What must be fixed before batch analysis.
6. What files should be created next, if any.

Do not ask to reorganise the repository unless there is a technical reason.
