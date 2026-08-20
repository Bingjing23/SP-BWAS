# Lancet 2024 Risk-Factor Comparison Handoff

This handoff summarizes the current window's work comparing Lancet 2024
dementia risk-factor priorities with the current SP-BWAS brainMapR results.

## Purpose

The question addressed here was:

> Which traits are described by the Lancet 2024 Commission as relatively
> important dementia risk factors, but appear weak in our current map-level
> brainMapR results?

The comparison is map-level only. It should be described as shared spatial
signal or map-level similarity, not as causal evidence.

## Source Files

- Lancet paper:
  `/Users/junzhou/Desktop/Side Project/PIIS0140673624012960.pdf`
- Current brainMapR summary directory:
  `/Users/junzhou/Desktop/Side Project/summary_clean_AVERAGE/`
- Main SP-BWAS working directory:
  `/Users/junzhou/Desktop/Side Project/SP-BWAS`

The reported brainMapR results use the `AVERAGE` reference panel. The main
metric used for this handoff is `rGM` from
`brainMapR::sumR2_regression_bivariate()`.

## PAF Definition Used Here

The PAF values in the figures and tables come from the Lancet 2024 Commission
paper, Table 1. They are not calculated from SP-BWAS data.

PAF means population attributable fraction. A simple single-risk-factor PAF is:

```text
PAF = Pe * (RR - 1) / (1 + Pe * (RR - 1))
```

where `Pe` is population exposure prevalence and `RR` is relative risk for
dementia.

The Lancet paper additionally accounts for clustering or overlap between risk
factors using communality estimates. The values used here are the paper's
rounded weighted PAF values, not raw unweighted PAF values and not a simple sum
of individual PAFs.

## Lancet 2024 Rounded Weighted PAF Ranking

| Rank | Risk factor | Rounded weighted PAF |
|---:|---|---:|
| 1 | Hearing loss | 7% |
| 2 | High LDL cholesterol | 7% |
| 3 | Less education | 5% |
| 4 | Social isolation | 5% |
| 5 | Depression | 3% |
| 6 | Traumatic brain injury | 3% |
| 7 | Air pollution | 3% |
| 8 | Physical inactivity | 2% |
| 9 | Smoking | 2% |
| 10 | Diabetes | 2% |
| 11 | Hypertension | 2% |
| 12 | Untreated vision loss | 2% |
| 13 | Obesity | 1% |
| 14 | Excessive alcohol consumption | 1% |

## Interpretation Rule for Current Results

For this review, weak current map-level signal was defined approximately as:

```text
abs(rGM) < 0.20
```

This is an effect-size rule. A result can have a small p-value but still be
called weak if the absolute `rGM` is small.

## Main Findings

### Lancet-Relevant Traits That Look Weak in Our Current Results

| Priority | Lancet risk factor | Current closest trait(s) | Current result summary |
|---:|---|---|---|
| 1 | Diabetes | `T2D` | AD vs HC `rGM = -0.033`, `p = 0.269`; max `abs(rGM)` across 8 AD-related maps is `0.067`. |
| 2 | Hypertension | `hyperTension` | AD vs HC `rGM = -0.045`, `p = 0.087`; max `abs(rGM)` across 8 AD-related maps is `0.133`. |
| 3 | Less education | `age_completed_education` | AD vs HC `rGM = -0.151`, `p = 3.1e-06`; max `abs(rGM)` is `0.176`. Statistically non-null but small in effect size. |
| 4 | Smoking burden | `pack_years`, `pack_years_adult_prop` | Dose-related smoking proxies are weak. `smoking_ever` is modest, with max `abs(rGM)` about `0.207`. |
| 5 | Excessive alcohol consumption | `Alcohol_week_avg`, `Alcohol_merge_phenotype`, `alcohol_intake_frequency` | Current alcohol proxies have weak map-level similarity, with max `abs(rGM)` about `0.166`. |
| 6 | Physical inactivity | Most `pa_*` proxies | Mixed evidence. Most proxies are weak-to-moderate. Avoid interpreting out-of-range physical-activity rGM values as biological signal. |

### Lancet-Relevant Traits With Clearer Current Map-Level Signal

| Lancet risk factor | Current closest trait(s) | Current result summary |
|---|---|---|
| High LDL cholesterol | `ldl_direct`, `hyperCholesterolemia` | `ldl_direct` shows clear signal, with AD vs HC `rGM = 0.344` and max `abs(rGM) = 0.444`. `hyperCholesterolemia` is weaker. |
| Social isolation | `loneliness_isolation`, `leisure_social_score` | `loneliness_isolation` shows clear signal, with AD vs HC `rGM = 0.608` and max `abs(rGM) = 0.700`. `leisure_social_score` is weaker. |
| Depression | `majorDepression` | Clear signal, with AD vs HC `rGM = 0.513` and max `abs(rGM) = 0.629`. |
| Obesity | `bmi` | Clear signal, with AD vs HC `rGM = 0.304` and max `abs(rGM) = 0.435`. |

### Lancet-Relevant Traits Not Directly Tested

No direct current BWAS trait was identified for:

- Hearing loss
- Traumatic brain injury
- Air pollution
- Untreated vision loss

These should be described as not tested, not as weak current results.

## Generated Shareable Tables

Full English comparison table:

```text
outputs/lancet_2024_risk_factor_comparison_for_supervisors.md
outputs/lancet_2024_risk_factor_comparison_for_supervisors.csv
```

Focused weak-signal table:

```text
outputs/lancet_2024_weak_signal_traits_for_supervisors.csv
```

## Generated Figures

Initial figure set:

```text
outputs/lancet_2024_supervisor_figures/
```

Preferred clean figure set:

```text
outputs/lancet_2024_supervisor_figures_clean/
```

Preferred files:

```text
outputs/lancet_2024_supervisor_figures_clean/figure_1_lancet_rank_lollipop.png
outputs/lancet_2024_supervisor_figures_clean/figure_1_lancet_rank_lollipop.svg
outputs/lancet_2024_supervisor_figures_clean/figure_2_lancet_paf_rgm_quadrant.png
outputs/lancet_2024_supervisor_figures_clean/figure_2_lancet_paf_rgm_quadrant.svg
```

Figure 1 is a row-wise Lancet ranking plot with PAF and max `abs(rGM)`.
Figure 2 is a quadrant-style scatter using numbered points and a right-side
key to avoid overlapping labels.

## Reproducible Plotting Script

The clean figures are generated by:

```text
scripts/07_plot_lancet_supervisor_figures.js
```

Run from the repository root:

```bash
node scripts/07_plot_lancet_supervisor_figures.js
```

The script writes SVG files and, if `rsvg-convert` is available, PNG files.

## Suggested Wording for Supervisors

The main mismatch between the Lancet 2024 risk-factor ranking and our current
brainMapR results is that several epidemiologically important dementia risk
factors, especially diabetes, hypertension, lower education, smoking burden,
and alcohol-related traits, show weak map-level similarity with the current
AD-related BWAS maps. In contrast, LDL cholesterol, loneliness/isolation,
depression, and BMI show clearer map-level signals. Hearing loss, traumatic
brain injury, air pollution, and untreated vision loss cannot be evaluated from
the current result set because no direct corresponding BWAS trait is available.

## Cautions

- Do not describe these results as causal.
- Do not describe them as SNP-level GWAS results.
- Do not interpret missing Lancet traits as negative results.
- Use cautious terms such as map-level similarity, shared spatial signal, and
  cross-map comparison.
- Physical-activity results require QC caution because one earlier proxy showed
  out-of-range rGM values and should not be treated as interpretable signal.
