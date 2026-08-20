# SP-BWAS Project Handoff

Last updated: 2026-08-20

## Current State

This handoff records the current project state after the supervisor meeting,
the targeted email to Elise, Elise's reply, inspection of the upstream
`ukb-AD-traits` phenotype code, and comparison of the original and v2 UKB
summary-statistics folders.

For the chronological project record from inception to the current state, use:

```text
PROJECT_HISTORY.md
```

The project remains a map-level BWAS integration analysis. Results should be
described as shared spatial signal, map-level similarity, or cross-map
comparison. Do not describe these results as causal evidence or SNP-level GWAS
findings.

## Main Local Paths

```text
SP-BWAS repo:
/Users/junzhou/Desktop/Side Project/SP-BWAS

Original UKB summary stats:
/Users/junzhou/Desktop/Side Project/SP-BWAS/SSTAT_BingJing

Elise v2 UKB summary stats:
/Users/junzhou/Desktop/Side Project/SP-BWAS/SSTAT_BingJing_v2

Upstream phenotype-code project:
/Users/junzhou/Desktop/Side Project/ukb-AD-traits

Upstream phenotype QC table:
/Users/junzhou/Desktop/Side Project/ukb-AD-traits/outputs/qc_v2/qc_summary.tsv
```

`SSTAT_BingJing_v2/` is ignored by git and should remain local data.

## Why Elise Was Contacted

The meeting focused on several possible interpretation and QC problems:

- Some genetic/map correlations exceeded the theoretical `[-1, 1]` range.
- Pneumonia and selected traits had incomplete or unusable summary-statistics
  files.
- Alcohol, smoking, and insomnia showed directions that looked protective or
  otherwise biologically unexpected.
- There was uncertainty about phenotype coding, former-drinker handling,
  covariates, and whether BMI was adjusted for.
- The analysis needed a clean story before sharing results with Michelle.

## Elise's Main Reply

Elise replied that:

- New summary statistics with the correct number of vertices were available on
  ownCloud.
- New/fixed files were available for `pneumonia_ever`, `pa_vigorous_time`,
  `infect_recency_days`, `stroke`, and related problematic traits.
- For some no-usable-output traits, the issue arose because the phenotype was
  part of the user's code and ended up with no patients/cases.
- Former/heavy former drinkers were excluded in the alcohol construction.
- `Alcohol_week_avg` is an average weekly-intake measure across beverages, not
  total weekly intake.
- `Alcohol_merge_phenotype` combines frequency/week and average intake.
- BMI was not included in the listed covariates.
- `smoking_ever` was part of the user's code.
- Sleep duration was not altered by Elise.

The practical reading is that Elise likely generated the v2 summary statistics
from the user's `ukb-AD-traits` phenotype code, but the exact code version or
commit should still be confirmed.

## v2 File Comparison

The original and v2 summary-statistics folders contain the same 71 trait names:

```text
old SSTAT_BingJing:    71 trait files
new SSTAT_BingJing_v2: 71 trait files
common trait names:    71
v2-only traits:        0
old-only traits:       0
byte-identical files:  62
changed/fixed files:   9
```

Changed or fixed v2 files:

| Trait | Old state | v2 state | Status |
| --- | --- | --- | --- |
| `Menopause` | Header only | Complete 654003-row map | File fixed |
| `hrt_ever_used` | Header only | Complete 654003-row map | File fixed |
| `pneumonia_ever` | 19730 rows | Complete 654003-row map | File fixed |
| `sepsis_ever` | 128581 rows | Complete 654003-row map | File fixed |
| `diabetes_doctor_dx` | 255667 rows | Complete 654003-row map | File fixed |
| `other_alcohol_weekly_intake` | 613514 rows | Complete 654003-row map | File fixed |
| `stroke` | Complete rows but unusable `b=0`, `se=0`, `p=NA` | Complete map with real statistics | File fixed |
| `pa_vigorous_time` | Complete rows but unstable effect scale | Complete map with normal scale | File fixed |
| `Alcohol_merge_phenotype` | Row-name/format issue | Clean `Probe` format | File fixed |

The 62 unchanged files can be treated as file-level usable, but unchanged files
are not automatically interpretation-ready.

## Upstream Phenotype-Code Findings

The upstream project contains the phenotype definitions:

```text
/Users/junzhou/Desktop/Side Project/ukb-AD-traits
```

Relevant source files:

```text
ukbtraits/registry.py
ukbtraits/traits.py
ukbtraits/utils.py
TRAITS_RULES.md
TRAIT_COLUMNS.md
fields.csv
fields_description.csv
```

Confirmed coding:

| Trait | Definition |
| --- | --- |
| `smoking_ever` | Misleading label. It is current smoker versus never smoker, not standard ever-smoked. `20116 == 2 -> 1`, `20116 == 1 -> NA`, `20116 == 0` and `20160 == 0 -> 0`, else `NA`. |
| `insomnia` | UKB field `1200`, prefer instance 2 then fallback to instance 0. Negative UKB missing codes are `NA`. Higher values mean more frequent insomnia. |
| `sleep_duration` | UKB field `1160`, instance 2. Treated as scalar raw hours unless a new upstream nonlinear/categorical phenotype is generated. |

Implication:

- Do not label `smoking_ever` as conventional ever smoking in slides or paper.
  Use current-vs-never smoking unless the phenotype is regenerated.
- Insomnia's unexpected direction is not explained by a simple coding reversal.
- Sleep duration should not be interpreted as a monotonic risk/protective
  variable without a nonlinear or categorical sensitivity analysis.

## qc_v2 Traits Missing From Both Summary-Statistics Folders

The upstream `qc_v2/qc_summary.tsv` contains 72 phenotype definitions. The
following 15 are present in `qc_v2` but absent from both original and v2
summary-statistics folders:

| Trait | Current decision |
| --- | --- |
| `age_at_menopause` | Ask Elise if this should be generated separately; not equivalent to `Menopause`. |
| `air_no2_2010` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `air_nox_2010` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `air_pm10_2010` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `air_pm25_10_2010` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `air_pm25_2010` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `anxiety_strict_preimg` | Ask whether strict psychiatric versions were intentionally omitted. |
| `bd_strict_preimg` | Ask whether strict psychiatric versions were intentionally omitted. |
| `mdd_strict_preimg` | Ask whether strict psychiatric versions were intentionally omitted. |
| `ptsd_strict_preimg` | Ask whether strict psychiatric versions were intentionally omitted. |
| `scz_strict_preimg` | Ask whether strict psychiatric versions were intentionally omitted. |
| `hearing_aid_user` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `hearing_loss_current` | Optional Lancet-relevant extension; currently no BWAS summary stats. |
| `stroke_i64_preimg` | Ask if this pre-imaging stroke phenotype should replace or supplement `stroke`. |
| `vision_problems_current` | Optional Lancet-relevant extension; currently no BWAS summary stats. |

These are not negative results and not previously passed summary stats. They
are missing upstream BWAS outputs unless Elise confirms otherwise.

## Remaining Questions for Elise

Ask only targeted remaining questions:

1. Were the v2 summary statistics generated from the exact current
   `ukb-AD-traits` code/version?
2. Were the 15 `qc_v2` traits missing from both summary-statistics folders
   intentionally omitted, or should any be generated separately?
3. Should `stroke_i64_preimg` be generated and used instead of, or alongside,
   the current `stroke` file?
4. Should strict psychiatric pre-imaging traits be generated, or should the
   current clinical/Rmd psychiatric traits remain the primary analysis set?
5. Should hearing, vision, and air-pollution traits be generated for a
   Lancet-aligned extension, or left out of the current paper scope?
6. Confirm that the listed covariates did not include BMI, and document that
   metabolic/vascular interpretations are not BMI-adjusted.

## Immediate Next Steps

1. Make v2 the active UKB summary-statistics input source.
   - Conservative option: back up `SSTAT_BingJing/`, then replace the 9 changed
     files with their `SSTAT_BingJing_v2/` versions.
   - Cleaner pipeline option: edit the manifest/audit scripts to accept a
     `--risk-dir` argument and run directly from `SSTAT_BingJing_v2/`.
2. Rerun input audit against v2 and save a separate post-v2 audit directory.
3. Regenerate manifests for AVERAGE and UKB reference panels.
4. Recreate derived `Probe` inputs with `--force`; otherwise old derived inputs
   could silently persist.
5. Rerun clean AVERAGE and UKB brainMapR batches.
6. Recollect and replot outputs.
7. Specifically check whether the old problems are resolved:
   - `pneumonia_ever` UKB-panel failures
   - `pneumonia_ever` out-of-range rGM
   - `pa_vigorous_time` out-of-range rGM
   - `stroke` no-usable-statistics issue
   - `Menopause` and `hrt_ever_used` header-only issue
8. Prepare Michelle-facing slides only after v2 QC and reruns show which traits
   are fixed, excluded, or still flagged.

## Suggested Wording for Current Status

```text
Elise's v2 folder appears to be a refreshed full set of the same 71 UKB summary
statistics. Most files are byte-identical to the previous version, while nine
previously problematic files were replaced or fixed. The upstream phenotype
code confirms that smoking_ever is actually current-vs-never smoking and that
insomnia retains the UKB field 1200 direction. The next step is to rerun the
input audit, derived-input generation, and AVERAGE/UKB brainMapR batches using
v2, then verify whether the previous out-of-range and failed traits are
resolved.
```
