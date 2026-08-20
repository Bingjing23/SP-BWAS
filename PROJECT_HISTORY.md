# SP-BWAS Project History

Last updated: 2026-08-20

This file is the chronological master record for the SP-BWAS project. It
combines project motivation, data structure, upstream phenotype-code evidence,
brainMapR pipeline history, meeting-driven QC questions, Elise's reply, v2
summary-statistics comparison, and the current next steps.

Use this file for orientation. Use `PROJECT_HANDOFF.md` for the current
handoff and `PROJECT_TRACKER.md` for the live task tracker.

## 1. Project Inception and Scientific Aim

The project asks whether Alzheimer-related brain-wide association maps share
spatial signal with UK Biobank-derived risk-factor BWAS maps.

The target comparison is:

```text
rows    = Alzheimer-related BWAS maps
columns = UKB-derived risk-factor BWAS maps
cells   = brainMapR / sumR2 outputs such as rGM, SE, CI, p-value, and brain covariance
```

The analysis is map-level. It should not be described as a SNP-level GWAS, a
causal analysis, or a risk-factor ranking by causal importance.

Recommended language:

```text
shared spatial signal
map-level similarity
cross-map comparison
brain-wide association pattern
```

## 2. Main Data Inputs

The SP-BWAS repository contains:

```text
AlzDisease_LMM/
SSTAT_BingJing/
SSTAT_BingJing_v2/
scripts/
manifests/
outputs/
```

AD-side input:

```text
AlzDisease_LMM/*.linear_random_oscaFormat
```

UKB risk-factor input:

```text
SSTAT_BingJing/*.linear
SSTAT_BingJing_v2/*.linear
```

The 8 AD maps currently used in the clean analysis have confirmed fixed sample
sizes:

| AD map | Sample size |
| --- | ---: |
| `ADvsHC` | 3542 |
| `MCIvsHC` | 3976 |
| `Conversion1year` | 1257 |
| `Conversion2years` | 1199 |
| `Conversion3years` | 1031 |
| `Conversion4years` | 1285 |
| `Conversion5years` | 1197 |
| `MMSE` | 6981 |

UKB trait-side sample sizes are read from the `NMISS` column.

## 3. Upstream Phenotype Project

The upstream phenotype project is:

```text
/Users/junzhou/Desktop/Side Project/ukb-AD-traits
```

This project is the source of phenotype definitions used to generate many of
the UKB BWAS summary statistics. Its README describes a step-wise,
imaging-first pipeline:

1. Extract imaging participants.
2. Extract only required UKB fields listed in `fields.csv`.
3. Build derived traits using modular functions.
4. Run optional QC on the final trait table.

Important upstream files:

| File | Role |
| --- | --- |
| `fields.csv` | UKB field IDs and instances extracted from the raw UKB table |
| `fields_description.csv` | Field dictionary with descriptions |
| `ukbtraits/registry.py` | Exported trait list and builder order |
| `ukbtraits/traits.py` | Main trait construction functions |
| `ukbtraits/frailty.py` | Frailty-specific recoding |
| `ukbtraits/utils.py` | Missing-value helpers |
| `TRAITS_RULES.md` | Human-readable phenotype rules |
| `TRAITS_PSY.md` | Psychiatric phenotype rules |
| `TRAIT_COLUMNS.md` | Final trait-column dictionary |
| `outputs/qc_v2/qc_summary.tsv` | Upstream QC summary for 72 defined traits |

The script path invoked by `SSTAT_BingJing/16_BingJing.Rmd` matches this
upstream project structure:

```text
scripts.extract_imaging_eids
scripts.extract_fields_for_eids
scripts.build_final_table
outputs/final_traits.tsv
```

## 4. General Upstream Trait Rules

The upstream trait project uses these general rules:

- Most UKB negative codes such as `-1`, `-3`, and `-7` are treated as missing
  and converted to `NA`.
- Pre-imaging diagnosis logic uses imaging date as the index date.
- Events after imaging are set to `NA` rather than controls.
- Missing diagnosis dates do not automatically exclude a subject for ICD-based
  traits.
- Packaged-date placeholders such as `1900-01-01`, `1901-01-01`,
  `1902-02-02`, `1903-03-03`, `1909-09-09`, and `2037-07-07` are treated as
  invalid timing.
- Binary traits are encoded as `1`, `0`, and `NA`.

## 5. Key Phenotype Definitions Confirmed From Code

### Smoking

`smoking_ever` is a misleading label. It is not conventional ever-smoked.

Source logic:

```text
20116 instance 2: Smoking status
  0 = never
  1 = previous/former
  2 = current

20160 instance 2: Ever smoked
  0 = no
  1 = yes
```

Final coding:

```text
20116 == 2                    -> 1 current smoker
20116 == 1                    -> NA former smoker
20116 == 0 and 20160 == 0     -> 0 never smoker
else                          -> NA
```

Interpretation:

```text
current smoker vs never smoker
```

Do not label this as standard ever-smoked in slides, figures, or manuscript
text unless the phenotype is regenerated.

Other smoking traits:

```text
pack_years             = UKB field 20161, instance 2
pack_years_adult_prop  = UKB field 20162, instance 2
```

### Sleep

Sleep traits in `ukbtraits/registry.py`:

```text
sleep_duration = field 1160, instance 2
insomnia       = field 1200, prefer instance 2, fallback instance 0
daytime_napping = field 1220, instance 2
```

`insomnia` retains the UKB field 1200 coding direction:

```text
1 = Never/rarely
2 = Sometimes
3 = Usually
```

Higher values therefore mean more frequent insomnia. Negative missing codes are
converted to `NA`.

Interpretation:

- The unexpected negative/protective-looking insomnia direction is not explained
  by a simple coding reversal.
- Sleep duration is a scalar raw-hours phenotype unless a nonlinear or
  categorical upstream phenotype is generated.
- Because sleep duration may be U-shaped, do not interpret the current linear
  map as monotonic risk or protection without sensitivity analysis.

### Alcohol

Alcohol traits in the upstream code include:

```text
alcohol_intake_frequency
red_wine_weekly_intake
white_wine_weekly_intake
beer_cider_weekly_intake
spirits_weekly_intake
fortified_wine_weekly_intake
other_alcohol_weekly_intake
ever_drinker
```

General alcohol rules:

- Weekly intake is used if present.
- Monthly intake is converted to weekly by dividing by 4.
- If weekly and monthly/4 are both present, their average is used.
- Values greater than 20 are capped at 20 for weekly intake traits.
- Former drinkers are set to `NA`.
- `3731` is only collected when `1558 == 6`; if `1558 == 6` and `3731` is
  missing, alcohol traits are set to `NA`.

`ever_drinker`:

```text
yes if 1558 in 1..5 or any combined weekly field > 0
no if 1558 == 6, 3731 == 0, and all combined weekly fields are 0/NA
else NA
```

Elise's three additional alcohol variables:

```text
Alcohol_week_avg
Alcohol_intake_frequency
Alcohol_merge_phenotype
```

Elise clarified that:

- Former/heavy former drinkers were excluded.
- `Alcohol_week_avg` is mean weekly intake across beverage types, not total
  weekly intake.
- `Alcohol_merge_phenotype` combines frequency/week and average intake.

Interpretation:

- Do not describe alcohol results as protective without strong caveats.
- Former-drinker exclusion and current-drinker conditioning can create selection
  or collider-like interpretation problems.
- Composite alcohol measures need clear labels.

### Physical Activity

Physical activity traits use raw fields unless future code changes them:

```text
pa_moderate_days_10min = 884, instance 2
pa_vigorous_days_10min = 904, instance 2
pa_walk_days_10min     = 864, instance 2
pa_moderate_duration   = 894, instance 2
pa_light_time          = 104920, instance 2
pa_moderate_time       = 104910, instance 2
pa_vigorous_time       = 104900, instance 2
```

`pa_vigorous_time` code values are:

```text
0, 1, 13, 35, 57, 79, 912, 1200
```

This field drove earlier out-of-range concerns. Elise's v2 file appears to fix
the scale problem at the summary-statistics level, but it still needs a
post-v2 brainMapR rerun to verify that out-of-range rGM estimates disappear.

### Infection

Infection traits use ICD-10 codes from:

```text
41202 with dates in 41262
41270 with dates in 41280
```

Main infection definitions:

```text
infect_any_ever: A40-A41, J12-J18, J20-J22, N30, N39.0, L00-L08
sepsis_ever: A40-A41, R65.2
pneumonia_ever: J12-J18
lrti_ever: J20-J22
uti_ever: N30, N39.0
ssti_ever: L00-L08
infect_count_total: total matched infection codes
infect_count_5y: matched infection codes within 5 years before imaging
infect_recency_days: days from imaging to most recent pre-imaging infection
```

Interpretation:

- Pneumonia's old file had a clear row-count problem and UKB-panel failure.
- v2 fixes the pneumonia and sepsis file row counts.
- `infect_recency_days` remains byte-identical between old and v2 and still
  needs cautious interpretation because of missingness and scale.

### Psychiatric Traits

The upstream project defines strict and broad psychiatric traits.

Strict traits:

```text
scz_strict_preimg
bd_strict_preimg
anxiety_strict_preimg
ptsd_strict_preimg
mdd_strict_preimg
```

Broad traits:

```text
scz_broad_preimg
bd_broad_preimg
anxiety_broad_preimg
ptsd_broad_preimg
mdd_broad_preimg
```

The final summary-statistics folders include broad and several Rmd/manual
clinical psychiatric traits, but the strict psychiatric traits are present in
`qc_v2` and absent from both summary-statistics folders.

### Cardiometabolic and Vascular

Upstream definitions include:

```text
diabetes_doctor_dx = field 2443, prefer instance 2 then 0
t2d_e11_preimg = ICD-10 E11 pre-imaging
diabetes_preimg = doctor diagnosis or E11 pre-imaging
stroke_i64_preimg = packaged I64 pre-imaging flag
htn_i10_preimg = packaged I10 pre-imaging flag
```

The summary-statistics folders include `stroke`, but not `stroke_i64_preimg`.
This is an unresolved phenotype-choice issue.

### Hearing, Vision, Air Pollution

The upstream `qc_v2` phenotype set includes hearing, vision, and air-pollution
traits, but neither original nor v2 summary-statistics folders contain their
BWAS outputs. They are optional Lancet-aligned extensions, not current negative
results.

## 6. Initial Downstream Audit and Pilot

The initial repository audit established:

- AD-side files are in `AlzDisease_LMM/`.
- UKB risk-factor files are in `SSTAT_BingJing/`.
- UKB files generally use `Voxel`, while brainMapR expects `Probe`.
- Raw files must not be modified in place.
- Derived `Probe` files are created under `outputs/batch/derived_inputs/`.
- AD-side sample sizes are supplied as fixed numeric values.
- UKB-side sample sizes use `NMISS`.

The official pilot used:

```text
AD map:    BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat
UKB trait: sstat_FS_All_moda_total_hyperTension.linear
Panel:     AVERAGE
```

The pilot succeeded, establishing that the brainMapR route was executable.

## 7. Full Exploratory Batch

The exploratory batch used:

```text
8 AD maps x 66 UKB risk-factor maps = 528 jobs
```

Outcome:

```text
488 successful
40 failed
```

The 40 failures corresponded to 5 traits across all 8 AD maps:

```text
htn_i10_preimg
med_20003_n_distinct_i2
periodontal_k05_preimg
psy_any_strict_preimg
stroke
```

These failures were not random compute failures. The raw rows showed no usable
association statistics such as:

```text
b = 0
se = 0
p = NA
```

## 8. Clean AVERAGE and UKB Reruns

The clean analysis used:

```text
8 AD maps x 44 cleaned UKB risk-factor maps = 352 pairs
```

Clean AVERAGE:

```text
352 / 352 succeeded
```

Clean UKB sensitivity:

```text
344 / 352 succeeded
```

All 8 UKB failures corresponded to `pneumonia_ever`, with:

```text
object 'BWASsignif2' not found
```

The pre-v2 clean AVERAGE outputs had:

```text
126 Bonferroni/FWER-significant pairs
13 out-of-range rGM estimates
```

Out-of-range estimates were concentrated in:

```text
pa_vigorous_time
pneumonia_ever
infect_recency_days
```

## 9. Supervisor Meeting Interpretation Issues

The meeting raised these main issues:

1. `|rGM| > 1` can arise from unstable covariance, tiny variances, constrained
   variance but unconstrained covariance behavior, or pathological inputs.
2. Pneumonia had missing rows and UKB-panel failures.
3. Alcohol, smoking, insomnia, and sleep showed directions that were not
   straightforward to interpret.
4. Alcohol could be affected by former-drinker exclusion and selection into
   current drinkers.
5. Sleep duration may be U-shaped, so a linear map can be misleading.
6. Covariates were uncertain, especially whether BMI was included.
7. The project needs slides for Michelle only after QC flags are resolved or
   clearly documented.

## 10. Elise Reply and Current Interpretation

Elise replied that:

- Updated summary statistics with the right number of vertices were available.
- Problematic files such as `pneumonia_ever`, `stroke`, `pa_vigorous_time`, and
  `infect_recency_days` had v2 files or updates.
- Former/heavy former drinkers were excluded in alcohol phenotypes.
- `Alcohol_week_avg` is mean weekly intake across beverage types, not total
  weekly intake.
- `Alcohol_merge_phenotype` combines frequency/week and average intake.
- `smoking_ever` was part of the user's code.
- Sleep duration was not changed.
- The listed covariates did not include BMI.

Current reading:

```text
Elise likely generated v2 from the user's ukb-AD-traits code, but the exact
version should still be confirmed.
```

## 11. v2 Summary-Statistics Comparison

Local v2 folder:

```text
/Users/junzhou/Desktop/Side Project/SP-BWAS/SSTAT_BingJing_v2
```

Comparison with original:

```text
old SSTAT_BingJing:    71 trait files
new SSTAT_BingJing_v2: 71 trait files
common trait names:    71
v2-only traits:        0
old-only traits:       0
byte-identical files:  62
changed/fixed files:   9
```

Changed/fixed v2 files:

| Trait | Old state | v2 state |
| --- | --- | --- |
| `Menopause` | Header only | Complete 654003 rows |
| `hrt_ever_used` | Header only | Complete 654003 rows |
| `pneumonia_ever` | 19730 rows | Complete 654003 rows |
| `sepsis_ever` | 128581 rows | Complete 654003 rows |
| `diabetes_doctor_dx` | 255667 rows | Complete 654003 rows |
| `other_alcohol_weekly_intake` | 613514 rows | Complete 654003 rows |
| `stroke` | Complete rows but `b=0`, `se=0`, `p=NA` | Complete map with real statistics |
| `pa_vigorous_time` | Complete rows but unstable scale | Complete map with normal scale |
| `Alcohol_merge_phenotype` | Extra row-name/format issue | Clean `Probe` format |

Interpretation:

- v2 is a refreshed full folder, not a folder of new trait names.
- The 62 identical files can be treated as file-level already usable.
- The 9 changed files should be rerun through brainMapR before interpretation.

## 12. qc_v2 Traits Missing From Summary Stats

These upstream-defined traits are present in
`ukb-AD-traits/outputs/qc_v2/qc_summary.tsv` but absent from both original and
v2 summary-statistics folders:

```text
age_at_menopause
air_no2_2010
air_nox_2010
air_pm10_2010
air_pm25_10_2010
air_pm25_2010
anxiety_strict_preimg
bd_strict_preimg
hearing_aid_user
hearing_loss_current
mdd_strict_preimg
ptsd_strict_preimg
scz_strict_preimg
stroke_i64_preimg
vision_problems_current
```

These are not negative results. They are missing upstream BWAS summary
statistics unless Elise confirms intentional omission.

## 13. Current Next Steps

1. Make v2 the active input source.
   - Cleaner option: add `--risk-dir` support to audit/manifest scripts.
   - Conservative option: back up old files and replace the 9 changed files.
2. Rerun input audit against v2.
3. Regenerate manifests.
4. Force-regenerate derived `Probe` files.
5. Rerun clean AVERAGE brainMapR.
6. Rerun clean UKB brainMapR.
7. Recollect and replot outputs.
8. Check whether previous problems are resolved:
   - `pneumonia_ever` UKB failure
   - `pneumonia_ever` out-of-range rGM
   - `pa_vigorous_time` out-of-range rGM
   - `stroke` no-usable-statistics issue
   - `Menopause` and `hrt_ever_used` header-only issue
9. Ask Elise the remaining targeted questions:
   - Was v2 generated from the exact current `ukb-AD-traits` version?
   - Were the 15 `qc_v2`-present but summary-statistics-missing traits
     intentionally omitted?
   - Should `stroke_i64_preimg` be generated and used instead of or alongside
     current `stroke`?
   - Should strict psychiatric, hearing, vision, and air-pollution traits be
     generated for this project or treated as optional extensions?
10. Prepare Michelle-facing slides after v2 rerun QC.

