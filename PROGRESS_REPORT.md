# AD-BWAS Project Progress Report

Last updated: 2026-05-08

## Project Aim

This project tests whether Alzheimer-related brain-wide association maps share
spatial signal with UK Biobank-derived risk-factor BWAS maps.

The core question is:

**Which risk-factor-associated brain patterns are most similar to
Alzheimer-related brain patterns?**

This is a brain-map-level comparison. It is not a SNP-level GWAS and should not
be interpreted causally.

## Current Analysis Scope

The current reviewed analysis uses:

```text
8 Alzheimer-related BWAS maps
x
44 cleaned UKB risk-factor BWAS maps
=
352 pairwise comparisons
```

The 8 Alzheimer-related maps have confirmed fixed sample sizes:

| Alzheimer-related map | Sample size |
| --- | ---: |
| AD vs HC | 3542 |
| MCI vs HC | 3976 |
| Conversion within 1 year | 1257 |
| Conversion within 2 years | 1199 |
| Conversion within 3 years | 1031 |
| Conversion within 4 years | 1285 |
| Conversion within 5 years | 1197 |
| MMSE | 6981 |

The cleaned UKB risk-factor set follows Baptiste's feedback:

- Specific alcohol beverage variables were removed.
- Global alcohol measures were retained.
- `T2D` was retained as the selected diabetes variable.
- More clinical psychiatric diagnoses were retained.
- Broad/pre-imaging psychiatric definitions were removed.
- Currently unusable or pending upstream BWAS traits were excluded from the
  clean rerun.

## Completed Work

### Package and Environment

- Resolved the `GFA::ldsc_rg()` dependency issue by using the GitHub version of
  `GFA` from `jean997/GFA`.
- Confirmed `brainMapR::sumR2_regression_bivariate()` is available.
- Confirmed the installed BrainMapR function includes `varConstrained`.
- Confirmed `varConstrained` is used inside the function body.
- Updated the batch wrapper to pass `varConstrained = TRUE` explicitly.
- Confirmed all clean AVERAGE rerun metadata files record
  `varConstrained = TRUE`.

### Input Preparation

- Audited AD BWAS files and UKB risk-factor BWAS files.
- Preserved all raw inputs.
- Created derived UKB inputs with `Voxel` renamed to `Probe`.
- Fixed UKB files where an extra row-name column shifted the data columns.
- Generated manifest-driven design tables for:
  - Full exploratory pairwise analysis.
  - Clean AVERAGE rerun.
  - Clean UKB reference-panel sensitivity rerun.

### Exploratory Full Batch

The first full exploratory batch used the broader default risk-factor set:

```text
8 AD maps x 66 UKB risk-factor maps = 528 jobs
```

Outcome:

```text
488 successful jobs
40 failed jobs
```

The 40 failures corresponded to 5 UKB traits across all 8 AD maps:

```text
htn_i10_preimg
med_20003_n_distinct_i2
periodontal_k05_preimg
psy_any_strict_preimg
stroke
```

These traits had no usable association statistics, such as `b = 0`, `se = 0`,
and `p = NA`, and were excluded from the clean rerun pending upstream checks.

### Clean Rerun

Clean AVERAGE reference-panel rerun:

```text
352 / 352 successful
```

Clean UKB reference-panel sensitivity rerun:

```text
344 / 352 successful
8 failed
```

The 8 UKB-panel failures all corresponded to:

```text
pneumonia_ever
```

This trait succeeded under the AVERAGE panel but failed under the UKB panel
with:

```text
object 'BWASsignif2' not found
```

This looks reference-panel or input-specific rather than a general workflow
failure.

### Figures and Tables

The cleaned results were collected into:

```text
outputs/batch/summary_clean_AVERAGE/
outputs/batch/summary_clean_UKB/
```

Updated outputs include:

- Vertical rGM heatmap with UKB traits on the Y axis and AD maps on the X axis.
- Vertical p-value heatmap.
- Top stable rGM forest plot.
- AVERAGE-vs-UKB reference-panel sensitivity scatter plot.
- Bonferroni/FWER-based significance markers.
- PNG and TIFF exports.
- Readable trait and AD labels.

Current clean AVERAGE plotting summary:

```text
AD maps: 8
Risk-factor maps: 44
Collected pairs: 352
Bonferroni/FWER-significant pairs: 126
Out-of-range rGM estimates: 13
Top stable associations written: 30
```

## Problems Identified and Current Interpretation

### Header-Only Traits

The following traits contained only headers and no usable data:

```text
Menopause
hrt_ever_used
```

These are pending upstream checks. Baptiste also suggested adding:

```text
3581 Age at menopause
```

### All-Zero or Unusable BWAS Traits

The following traits had unusable BWAS outputs in the exploratory full batch:

```text
htn_i10_preimg
med_20003_n_distinct_i2
periodontal_k05_preimg
psy_any_strict_preimg
stroke
```

Potential explanations include zero cases in the imaging subset, pre-imaging
phenotype derivation issues, or sex-covariate issues for sex-specific traits.
These need Elise's upstream confirmation.

### Remaining Out-of-Range rGM Estimates

Even after explicitly setting `varConstrained = TRUE`, the clean AVERAGE rerun
still has 13 estimates outside the theoretical `[-1, 1]` range.

Checks completed:

- BrainMapR includes `varConstrained`.
- The function body uses `varConstrained`.
- The batch wrapper passes `varConstrained = TRUE`.
- All 352 clean AVERAGE run metadata files record `varConstrained = TRUE`.
- Timestamps confirm the outputs are from the new rerun, not old results.

The 13 out-of-range estimates are concentrated in three traits:

| Trait | Out-of-range count | Current interpretation |
| --- | ---: | --- |
| `pa_vigorous_time` | 8 | Complete map, but very large beta/SE scale; likely unstable |
| `pneumonia_ever` | 4 | Incomplete map and UKB-panel failure; likely upstream/input issue |
| `infect_recency_days` | 1 | Complete map, but very large beta/SE scale; likely unstable |

Input checks:

- `pa_vigorous_time` has a complete map but a large beta range
  approximately `-370` to `368`.
- `infect_recency_days` has a complete map but an even larger beta range
  approximately `-1522` to `2030`.
- `pneumonia_ever` has only about `19.7k` rows rather than the expected
  approximately `654k`, and failed across all 8 AD maps under the UKB
  reference panel.

Current interpretation:

```text
The remaining |rGM| > 1 estimates are unlikely to be caused by missing
varConstrained or stale outputs. They appear concentrated in specific unstable
or incomplete trait inputs.
```

These traits should be flagged or excluded from the main interpretable figure
unless Baptiste recommends a BrainMapR-side fix.

## Files to Share with Collaborators

Recommended figure attachments:

```text
outputs/batch/summary_clean_AVERAGE/figures/figure_1_vertical_rGM_heatmap.png
outputs/batch/summary_clean_AVERAGE/figures/figure_2_vertical_pvalue_heatmap.png
outputs/batch/summary_clean_AVERAGE/figures/figure_3_top_rGM_forest.png
outputs/batch/summary_clean_AVERAGE/figures/figure_4_reference_panel_sensitivity.png
```

Recommended table attachments:

```text
outputs/batch/summary_clean_AVERAGE/top_associations_bonferroni.tsv
outputs/batch/summary_clean_AVERAGE/out_of_range_rGM.tsv
outputs/batch/summary_clean_AVERAGE/reference_panel_rGM_comparison.tsv
outputs/batch/ukb_clean_failed_pairs.tsv
```

## Next Steps

1. Send Baptiste the BrainMapR / `varConstrained` QC update.
2. Ask whether `pa_vigorous_time`, `infect_recency_days`, and
   `pneumonia_ever` should be excluded from the main interpretable figure or
   handled with another BrainMapR-side fix.
3. Follow up with Elise on:
   - `Menopause`
   - `hrt_ever_used`
   - `3581 Age at menopause`
   - `htn_i10_preimg`
   - `med_20003_n_distinct_i2`
   - `periodontal_k05_preimg`
   - `psy_any_strict_preimg`
   - `stroke`
   - `pneumonia_ever`
4. Prepare the main internal-update figure set using the clean AVERAGE outputs.
5. Use the UKB sensitivity results to report robustness of rGM estimates across
   reference panels.
6. Keep the manifest-driven pipeline ready so fixed or new upstream traits can
   be rerun and appended without reorganizing the project.
