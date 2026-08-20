# SP-BWAS Project Tracker

Last updated: 2026-08-20

For the chronological project record and upstream phenotype-code evidence, use:

```text
PROJECT_HISTORY.md
```

Status labels:

```text
DONE       completed and documented
READY      ready to run or interpret after routine rerun
PENDING    waiting for rerun, collaborator confirmation, or decision
BLOCKED    cannot proceed without new upstream file or explicit decision
OPTIONAL   outside the current core analysis unless scope expands
```

## Workstream Tracker

| Workstream | Status | Current evidence | Next action |
| --- | --- | --- | --- |
| Original clean AVERAGE analysis | DONE | 352/352 clean AVERAGE pairs succeeded in previous run | Supersede with v2 rerun |
| Original clean UKB sensitivity | DONE | 344/352 succeeded; all 8 failures were `pneumonia_ever` | Supersede with v2 rerun |
| Elise v2 file receipt | DONE | `SSTAT_BingJing_v2/` exists locally and is git-ignored | Use as active input source |
| v2 vs old file comparison | DONE | 71 old, 71 v2, 71 common; 62 identical and 9 changed | Track changed files below |
| Upstream phenotype-code lookup | DONE | `ukb-AD-traits` inspected | Use confirmed coding in interpretation |
| `smoking_ever` coding | DONE | Current-vs-never smoker; former smoker is `NA` | Relabel in plots/slides/manuscript |
| `insomnia` coding | DONE | UKB field 1200, instance 2 then 0, higher = more insomnia | Interpret odd direction cautiously |
| `sleep_duration` functional form | PENDING | UKB field 1160 raw scalar; U-shaped concern remains | Consider nonlinear/categorical upstream BWAS |
| BMI covariate status | DONE | Elise listed covariates without BMI | Note results are not BMI-adjusted |
| v2 active input decision | PENDING | Scripts currently default to `SSTAT_BingJing/` | Choose copy/replace or add source-dir support |
| Post-v2 input audit | PENDING | Not yet run | Run audit against v2 |
| Post-v2 manifest generation | PENDING | Not yet run | Generate AVERAGE and UKB clean manifests |
| Post-v2 derived Probe inputs | PENDING | Old derived inputs may persist | Recreate with `--force` |
| Post-v2 AVERAGE rerun | PENDING | Not yet run | Submit clean AVERAGE jobs |
| Post-v2 UKB rerun | PENDING | Not yet run | Submit clean UKB jobs |
| Post-v2 result collection | PENDING | Not yet run | Run collector |
| Post-v2 figures | PENDING | Not yet run | Regenerate heatmaps, forest plot, sensitivity plot |
| Michelle slides | PENDING | Need v2 rerun status | Build after v2 QC is resolved |
| Elise follow-up | PENDING | Only targeted remaining questions needed | Ask about exact code version and missing qc_v2 traits |

## v2 Changed or Fixed Files

| Trait | Domain | Old issue | v2 evidence | Status | Next action |
| --- | --- | --- | --- | --- | --- |
| `Menopause` | Hormonal/sex-specific | Header-only file | 654003 rows, real statistics | READY | Include only if sex-specific interpretation is wanted |
| `hrt_ever_used` | Hormonal/sex-specific | Header-only file | 654003 rows, real statistics | READY | Include only if sex-specific interpretation is wanted |
| `pneumonia_ever` | Infection | 19730 rows; UKB-panel failures; out-of-range AVERAGE rGM | 654003 rows, real statistics | READY | Rerun AVERAGE/UKB and check failures/out-of-range |
| `sepsis_ever` | Infection | 128581 rows | 654003 rows, real statistics | READY | Rerun and compare effect size stability |
| `diabetes_doctor_dx` | Vascular/metabolic | 255667 rows | 654003 rows, real statistics | READY | Usually keep as sensitivity because `T2D` is main diabetes trait |
| `other_alcohol_weekly_intake` | Alcohol | 613514 rows | 654003 rows, real statistics | READY | Still beverage-specific; likely not main clean trait |
| `stroke` | Vascular/metabolic | 654003 rows but `b=0`, `se=0`, `p=NA` | 654003 rows, real statistics | READY | Rerun; decide whether to include main clean set |
| `pa_vigorous_time` | Physical activity | Large unstable beta scale and out-of-range rGM | 654003 rows, normal scale | READY | Rerun and verify out-of-range issue is resolved |
| `Alcohol_merge_phenotype` | Alcohol | Extra row-name/format issue | 654003 rows, clean `Probe` format | READY | Rerun; still interpret cautiously due composite definition |

## Byte-Identical Files

These 62 files are byte-identical between `SSTAT_BingJing/` and
`SSTAT_BingJing_v2/`. They can be treated as file-level already usable, but
coding and biological interpretation rules still apply.

```text
Alcohol_week_avg
BD
PTSD
SCZ
T2D
age_completed_education
alcohol_intake_frequency
anxDis
anxiety_broad_preimg
bd_broad_preimg
beer_cider_weekly_intake
bmi
daytime_napping
diabetes_preimg
eid
ever_drinker
fortified_wine_weekly_intake
frailty_fi49
frailty_missing_count
htn_i10_preimg
hyperCholesterolemia
hyperTension
imd_england
income_before_tax
infect_any_ever
infect_count_5y
infect_count_total
infect_recency_days
insomnia
ldl_direct
leisure_social_score
loneliness_isolation
lrti_ever
majorDepression
mdd_broad_preimg
med_20003_n_distinct_i2
multimorbidity
pa_light_time
pa_moderate_days_10min
pa_moderate_duration
pa_moderate_duration_pilot
pa_moderate_time
pa_vigorous_days_10min
pa_walk_days_10min
pack_years
pack_years_adult_prop
periodontal_k05_preimg
phobAnx
psy_any_broad_preimg
psy_any_strict_preimg
psy_clean_ctrl_mask
ptsd_broad_preimg
red_wine_weekly_intake
scz_broad_preimg
sex
sleep_duration
smoking_ever
spirits_weekly_intake
ssti_ever
t2d_e11_preimg
uti_ever
white_wine_weekly_intake
```

Special note: several byte-identical files such as `htn_i10_preimg`,
`med_20003_n_distinct_i2`, `periodontal_k05_preimg`, and
`psy_any_strict_preimg` remain old no-usable-statistics style outputs with
`b=0`, `se=0`, and `p=NA`. They should stay excluded unless regenerated or
explicitly needed for audit.

## qc_v2 Traits Missing From Summary Stats

These 15 traits are present in the upstream phenotype QC table but absent from
both original and v2 summary-statistics folders.

| Trait | Priority | Status | Next action |
| --- | --- | --- | --- |
| `age_at_menopause` | Medium | BLOCKED | Ask Elise whether to generate; not equivalent to `Menopause` |
| `stroke_i64_preimg` | High | BLOCKED | Ask whether this should replace/supplement current `stroke` |
| `anxiety_strict_preimg` | Medium | BLOCKED | Ask whether strict psychiatric traits were omitted intentionally |
| `bd_strict_preimg` | Medium | BLOCKED | Ask whether strict psychiatric traits were omitted intentionally |
| `mdd_strict_preimg` | Medium | BLOCKED | Ask whether strict psychiatric traits were omitted intentionally |
| `ptsd_strict_preimg` | Medium | BLOCKED | Ask whether strict psychiatric traits were omitted intentionally |
| `scz_strict_preimg` | Medium | BLOCKED | Ask whether strict psychiatric traits were omitted intentionally |
| `hearing_loss_current` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `hearing_aid_user` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `vision_problems_current` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `air_no2_2010` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `air_nox_2010` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `air_pm10_2010` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `air_pm25_2010` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |
| `air_pm25_10_2010` | Optional/Lancet | BLOCKED | Optional extension; needs summary stats |

## Interpretation Tracker

| Trait or domain | Status | Rule for current interpretation |
| --- | --- | --- |
| Alcohol composites | PENDING | Do not call protective; former-drinker exclusion and composite scale can create selection/interpretation issues |
| `Alcohol_week_avg` | READY with caution | Mean weekly intake across beverages, not total weekly intake |
| `Alcohol_merge_phenotype` | READY with caution | Composite of frequency and intake; rerun after v2 format fix |
| `smoking_ever` | DONE | Relabel as current-vs-never smoking |
| `pack_years` | READY | Safer smoking burden measure than `smoking_ever` label |
| `insomnia` | DONE | Higher values mean more insomnia; odd direction is not coding reversal |
| `sleep_duration` | PENDING | Linear raw-hours map only; consider nonlinear/categorical sensitivity |
| `pa_vigorous_time` | PENDING | v2 likely fixes scale; must rerun brainMapR to verify |
| `infect_recency_days` | PENDING | Still byte-identical and small-sample/high-scale concern remains |
| `pneumonia_ever` | PENDING | v2 fixes row count; rerun to verify UKB failure and rGM range |
| `stroke` | PENDING | v2 fixes no-usable stats; decide current `stroke` vs `stroke_i64_preimg` |
| Menopause/HRT | PENDING | v2 fixes files; include only with clear sex-specific analysis decision |
| BMI adjustment | DONE | Not included in Elise-listed covariates; note this in interpretation |

## Post-v2 Rerun Checklist

1. Create a pre-v2 backup or use source-dir-aware scripts.
2. Run input audit against v2.
3. Generate v2 manifests.
4. Force-regenerate derived `Probe` inputs.
5. Submit clean AVERAGE brainMapR batch.
6. Submit clean UKB brainMapR batch.
7. Collect AVERAGE and UKB outputs.
8. Plot heatmaps, top associations, and reference-panel sensitivity.
9. Check `out_of_range_rGM.tsv`.
10. Check UKB failed-pairs table.
11. Update tracker with resolved/unresolved traits.
12. Prepare Michelle-facing slides.

## Targeted Elise Follow-up Draft

```text
Hi Elise,

Thanks again for the updated files. I compared the v2 folder with the previous
summary-statistics folder. The v2 folder appears to contain the same 71 trait
names, with nine previously problematic files replaced or fixed.

Could you please confirm whether the v2 summary statistics were generated from
the current ukb-AD-traits code/version?

I also noticed that several traits present in the ukb-AD-traits qc_v2 summary
are not present in either summary-statistics folder: age_at_menopause,
stroke_i64_preimg, the strict psychiatric traits, hearing/vision traits, and
air-pollution traits. Were these intentionally omitted, or should any of them
be generated separately?

Best,
Bingjing
```
