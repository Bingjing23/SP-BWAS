# AD-BWAS Project Slide Content Source

Last updated: 2026-05-28

This document is a detailed content source for preparing slides about the
AD-BWAS / brainMapR downstream integration project. It is written as a
slide-planning reference rather than as a methods script or technical run log.

## 1. Short Project Summary

This project compares Alzheimer-related brain-wide association study (BWAS)
maps with UK Biobank-derived risk-factor BWAS maps.

The main question is:

**Which Alzheimer-relevant risk-factor brain maps show shared spatial signal
with Alzheimer-related brain maps?**

The analysis is performed at the level of brain-wide association maps. It is
not a SNP-level GWAS, and it should not be interpreted as causal evidence that
a risk factor directly causes Alzheimer disease through a specific brain
region.

The final scientific output is a pairwise comparison matrix:

```text
Rows    = Alzheimer-related BWAS maps
Columns = UKB risk-factor BWAS maps
Cells   = brainMapR outputs, such as rGM / map-level correlation,
          sumR2-style shared signal, standard errors, p-values, or confidence
          intervals where available
```

## 2. One-Sentence Version

We are using brainMapR to test whether brain-wide association patterns for
Alzheimer-related traits spatially overlap with brain-wide association patterns
for UK Biobank risk factors such as hypertension, diabetes, cholesterol,
depression, smoking, sleep, infection, frailty, and socioeconomic traits.

## 3. Non-Technical Explanation

Each BWAS file can be thought of as a brain map. Every location in the map has
a number that summarizes how strongly that brain location is associated with a
trait.

Examples:

- An AD vs HC map summarizes the brain-wide association pattern for Alzheimer
  dementia cases compared with cognitively normal controls.
- A hypertension map summarizes the brain-wide association pattern associated
  with hypertension.
- A depression map summarizes the brain-wide association pattern associated
  with major depression.

The project asks whether two maps have similar spatial patterns. For example:

```text
Does the AD vs HC brain map look spatially similar to the hypertension brain map?
Does the MCI vs HC brain map look spatially similar to the diabetes brain map?
Do early conversion-to-AD maps share spatial signal with lifestyle or
psychiatric risk-factor maps?
```

If two maps have similar spatial patterns, this may suggest shared
brain-wide association structure. It does not prove causality.

## 4. Scientific Motivation

Alzheimer disease is influenced by a broad range of risk factors, including
vascular, metabolic, psychiatric, lifestyle, social, inflammatory, and frailty
related traits.

Many of these risk factors are also associated with brain structure. Rather
than testing one region at a time, this project compares whole-brain
association patterns.

The scientific motivation is to identify which risk-factor-associated brain
patterns are most aligned with Alzheimer-related brain patterns. This can help
prioritize risk-factor domains for follow-up analyses and interpretation.

## 5. Key Terms for Slides

### BWAS

Brain-wide association study. In this project, a BWAS summary-statistic file
contains association statistics across brain locations for a phenotype or
trait.

### AD-related BWAS map

A brain-wide association map for an Alzheimer-related outcome, such as AD vs
healthy controls, MCI vs healthy controls, conversion to AD dementia, or MMSE.

### Risk-factor BWAS map

A brain-wide association map for a risk factor or external trait, derived from
UK Biobank data.

### Map-level similarity

Similarity between two brain-wide association maps across brain locations.

### Shared spatial signal

Evidence that two maps show related spatial association patterns.

### rGM

Gray-matter map-level correlation or related brainMapR output summarizing
cross-map similarity. Use cautious language because exact interpretation
depends on the brainMapR implementation and reference panel.

### sumR2

A brainMapR output intended to summarize shared signal / explained variance
style information between maps.

## 6. Main Scientific Question

The main question is:

```text
Which Alzheimer-relevant risk-factor BWAS maps share spatial signal with
Alzheimer-related BWAS maps?
```

More specifically:

- Do vascular and metabolic risk-factor maps resemble AD-related maps?
- Do psychiatric maps, such as depression or anxiety, resemble AD-related maps?
- Do lifestyle maps, such as smoking, sleep, or alcohol use, resemble AD-related
  maps?
- Do infection, frailty, or socioeconomic risk-factor maps show shared
  spatial signal with early AD conversion maps?
- Are the patterns stronger for dementia status, MCI, conversion phenotypes, or
  cognitive measures such as MMSE?

## 7. Data Overview

### Alzheimer-Related Inputs

Raw AD-related BWAS input directory:

```text
AlzDisease_LMM/
```

Initial inventory:

```text
24 Alzheimer-related BWAS summary-statistic files
```

These include:

- AD vs healthy controls
- MCI vs healthy controls
- AD vs MCI
- Conversion to AD dementia within 1 to 5 years
- Alzheimer progression phenotypes
- MMSE
- Clinical Dementia Rating
- Functional Assessment Questionnaire
- Geriatric Depression Scale
- Logical memory measures
- RAVLT measures
- Maternal / paternal AD
- Neuropsychiatric Inventory Score

The AD meta-analysis files have the required brain-map columns, but they do
not contain `NMISS` because sample size was not carried forward in the
meta-analysis summary statistics.

Therefore, AD sample sizes are supplied as fixed numeric values.

### Confirmed AD Maps Used in Current Clean Analysis

The current analysis uses the 8 Alzheimer-related maps with confirmed sample
sizes:

| AD map | Description | Sample size |
| --- | --- | ---: |
| AD vs HC | AD dementia cases vs healthy / cognitively normal controls | 3542 |
| MCI vs HC | Mild cognitive impairment vs healthy / cognitively normal controls | 3976 |
| Conversion within 1 year | Converted to AD dementia within 1 year from MRI vs non-converters followed for the same period | 1257 |
| Conversion within 2 years | Converted within 2 years vs matched non-converters | 1199 |
| Conversion within 3 years | Converted within 3 years vs matched non-converters | 1031 |
| Conversion within 4 years | Converted within 4 years vs matched non-converters | 1285 |
| Conversion within 5 years | Converted within 5 years vs matched non-converters | 1197 |
| MMSE | Mini-Mental State Examination score | 6981 |

These 8 maps represent a range from clinical AD dementia status to earlier
conversion phenotypes and cognitive performance.

### UKB Risk-Factor Inputs

Raw UKB-derived BWAS input directory:

```text
SSTAT_BingJing/
```

Initial inventory:

```text
71 UKB-derived BWAS summary-statistic files
```

These maps cover broad Alzheimer-relevant risk-factor domains:

#### Vascular and Metabolic

Examples:

- Hypertension
- Type 2 diabetes
- Doctor-diagnosed diabetes
- Hypercholesterolemia
- LDL direct
- BMI
- Stroke

#### Psychiatric

Examples:

- Major depression
- Anxiety
- PTSD
- Bipolar disorder
- Schizophrenia
- Phobic anxiety
- Broad psychiatric phenotype candidates

#### Lifestyle

Examples:

- Ever smoking
- Pack years
- Sleep duration
- Insomnia
- Daytime napping
- Alcohol measures
- Physical activity measures

#### Social and Socioeconomic

Examples:

- Age completed education
- Income before tax
- Index of multiple deprivation
- Leisure / social score
- Loneliness / isolation

#### Frailty and Multimorbidity

Examples:

- Frailty index
- Frailty missing count
- Multimorbidity
- Number of distinct medication codes

#### Infection and Inflammatory

Examples:

- Any infection ever
- Infection count over 5 years
- Lifetime infection count
- Infection recency
- Lower respiratory tract infection
- Pneumonia
- Sepsis
- Skin / soft tissue infection
- Urinary tract infection
- Periodontal disease

#### Hormonal / Sex-Specific

Examples:

- Menopause
- HRT ever used

Some of these upstream risk-factor files were later excluded from cleaned
analyses because they were header-only, unstable, duplicate / overly specific,
or pending upstream confirmation.

## 8. Input Harmonization

### AD Files

AD meta-analysis files are already in the expected spatial-identifier format:

```text
Probe b se p
```

They do not contain `NMISS`, so sample size is supplied manually using
confirmed fixed values.

### UKB Files

UKB risk-factor files generally use:

```text
Chr Voxel bp Gene Orientation b se p NMISS
```

brainMapR expects `Probe` as the spatial identifier. Therefore, UKB files are
converted into derived copies with:

```text
Voxel -> Probe
```

Raw files are preserved. Only derived copies are used for brainMapR.

This is important for reproducibility: the original BWAS summary statistics
remain unchanged, while all brainMapR-ready inputs are generated
deterministically.

## 9. Analysis Method

The project uses:

```text
brainMapR::sumR2_regression_bivariate()
```

Each analysis compares two BWAS maps:

```text
one AD-related map
x
one UKB risk-factor map
```

For each pair, brainMapR estimates map-level shared signal using a reference
panel.

Primary reference panel:

```text
AVERAGE
```

Sensitivity reference panel:

```text
UKB
```

The AVERAGE panel is used as the primary panel because it was recommended as
the main analysis panel. The UKB panel is used as a sensitivity analysis to
check whether the same broad patterns are stable under a different reference
panel.

The current official route uses:

```text
brainMapR_1.1.0.9000
GFA from jean997/GFA
```

This resolved the earlier issue where the CRAN version of GFA did not contain
the function required by brainMapR.

## 10. Analysis Design Evolution

### Stage 1: Pilot

The first official pilot compared:

```text
AD vs HC
x
UKB hypertension
```

Settings:

```text
AD sample size: 3542
UKB sample size: NMISS
Reference panel: AVERAGE
```

Purpose:

- Confirm brainMapR and GFA load correctly.
- Confirm the AD fixed sample size approach works.
- Confirm the UKB `NMISS` approach works.
- Confirm the derived `Probe` inputs work.
- Confirm that the output structure is usable.

Status:

```text
Pilot completed successfully.
```

### Stage 2: Small Batch

After the successful pilot, a controlled small batch was prepared:

```text
2 AD maps x 6 UKB risk-factor maps = 12 pairwise comparisons
```

AD maps:

- AD vs HC
- MCI vs HC

Risk-factor maps:

- Hypertension
- Type 2 diabetes
- Hypercholesterolemia
- LDL direct
- Major depression
- Ever smoking

Purpose:

- Move beyond a single pilot pair.
- Test representative vascular / metabolic, psychiatric, and lifestyle risk
  factors.
- Confirm that multiple pairwise outputs can be generated and collected.
- Avoid launching the full grid before confirming that the pipeline is stable.

Status:

```text
Small batch was submitted to PBS / cluster.
```

### Stage 3: Exploratory Full Batch

The first broader full-batch design used:

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

Current interpretation:

These failures were not random. They were concentrated in specific UKB traits
with unusable association statistics, such as `b = 0`, `se = 0`, and `p = NA`.
These traits were excluded from the clean rerun pending upstream checks.

### Stage 4: Clean AVERAGE Rerun

After review and cleanup, the clean AVERAGE reference-panel rerun used:

```text
8 AD maps x 44 cleaned UKB risk-factor maps = 352 comparisons
```

Outcome:

```text
352 / 352 successful
```

This is the main clean analysis result.

### Stage 5: Clean UKB Sensitivity Rerun

The same clean 8 x 44 design was rerun using the UKB reference panel.

Outcome:

```text
344 / 352 successful
8 failed
```

The 8 failures all corresponded to:

```text
pneumonia_ever
```

This trait succeeded under the AVERAGE panel but failed under the UKB panel
with an error involving:

```text
object 'BWASsignif2' not found
```

Current interpretation:

The pneumonia issue appears to be reference-panel or input-specific rather than
a general pipeline failure.

## 11. Clean Analysis Scope

The current reviewed clean analysis focuses on:

```text
8 AD maps
x
44 cleaned UKB risk-factor maps
=
352 pairwise comparisons
```

This clean set follows Baptiste's feedback:

- Remove specific alcohol beverage variables.
- Retain global alcohol measures.
- Retain `T2D` as the selected diabetes variable.
- Retain more clinical psychiatric diagnoses.
- Remove broad or pre-imaging psychiatric definitions.
- Exclude currently unusable or pending upstream BWAS traits.
- Exclude header-only hormonal / sex-specific files unless regenerated.

This clean design is more suitable for presentation than the broader
exploratory design because it removes traits with known upstream or stability
issues.

## 12. Current Results Summary

### Clean AVERAGE Run

```text
AD maps: 8
Risk-factor maps: 44
Collected pairs: 352
Successful pairs: 352
```

### Clean UKB Sensitivity Run

```text
AD maps: 8
Risk-factor maps: 44
Expected pairs: 352
Successful pairs: 344
Failed pairs: 8
Failed trait: pneumonia_ever
```

### Clean AVERAGE Plotting Summary

Current clean AVERAGE plotting summary:

```text
AD maps: 8
Risk-factor maps: 44
Collected pairs: 352
Bonferroni / FWER-significant pairs: 126
Out-of-range rGM estimates: 13
Top stable associations written: 30
```

These values should be presented with caution. The number of significant pairs
is useful for reporting progress, but biological interpretation should focus on
patterns that remain stable, interpretable, and not driven by known unstable
traits.

## 13. Output Products

The clean results were collected into:

```text
outputs/batch/summary_clean_AVERAGE/
outputs/batch/summary_clean_UKB/
```

Main figure outputs:

```text
Figure 1: Vertical rGM heatmap
Figure 2: Vertical p-value heatmap
Figure 3: Top stable rGM forest plot
Figure 4: AVERAGE-vs-UKB reference-panel sensitivity scatter plot
```

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

## 14. Known Data and Interpretation Issues

### Header-Only Traits

The following UKB traits contained headers but no usable data rows:

```text
Menopause
hrt_ever_used
```

These are pending upstream checks.

Baptiste suggested adding:

```text
3581 Age at menopause
```

This may provide a usable menopause-related phenotype if upstream BWAS can be
generated.

### All-Zero or Unusable Traits

The following traits had unusable BWAS outputs in the exploratory full batch:

```text
htn_i10_preimg
med_20003_n_distinct_i2
periodontal_k05_preimg
psy_any_strict_preimg
stroke
```

Potential explanations:

- Too few or zero cases in the imaging subset.
- Phenotype derivation issue.
- Pre-imaging phenotype definition issue.
- Sex-covariate or sample-filtering issue for sex-specific traits.
- Upstream BWAS output generated but not statistically usable.

These traits should not be included in the main interpretable figure unless
they are regenerated or confirmed upstream.

### Pneumonia Under UKB Reference Panel

`pneumonia_ever` succeeded under the AVERAGE reference panel but failed under
the UKB reference panel.

Current interpretation:

This is likely input-specific or reference-panel-specific, not a general
failure of the workflow.

### Out-of-Range rGM Estimates

The clean AVERAGE rerun still contains 13 rGM estimates outside the theoretical
`[-1, 1]` range, even after explicitly setting `varConstrained = TRUE`.

The out-of-range estimates are concentrated in:

| Trait | Out-of-range count | Current interpretation |
| --- | ---: | --- |
| `pa_vigorous_time` | 8 | Complete map, but very large beta / SE scale; likely unstable |
| `pneumonia_ever` | 4 | Incomplete map and UKB-panel failure; likely upstream / input issue |
| `infect_recency_days` | 1 | Complete map, but very large beta / SE scale; likely unstable |

Current interpretation:

```text
The remaining |rGM| > 1 estimates are unlikely to be caused by missing
varConstrained or stale outputs. They appear concentrated in specific unstable
or incomplete trait inputs.
```

These traits should be flagged or excluded from the main interpretable figure
unless Baptiste recommends a BrainMapR-side fix.

## 15. How to Interpret Results

Recommended interpretation language:

- shared spatial signal
- map-level similarity
- gray-matter map-level correlation
- brain-wide association pattern overlap
- AD-related BWAS map
- risk-factor-associated BWAS map
- cross-map comparison

Avoid:

- "This risk factor causes AD."
- "This brain region mediates the causal effect."
- "These are SNP-level GWAS results."
- "This proves a mechanistic pathway."
- "This trait is definitively biologically linked to AD."

Correct framing:

```text
The analysis identifies risk-factor-associated brain-wide patterns that show
similarity to AD-related brain-wide patterns. These results can guide
prioritization and follow-up analyses, but they do not establish causality.
```

## 16. Suggested Slide Deck Structure

### Slide 1: Project Title

Suggested title:

```text
Map-Level Integration of Alzheimer-Related and Risk-Factor BWAS Maps
```

Subtitle:

```text
Using brainMapR to identify shared spatial signal between AD-related brain maps
and UKB-derived risk-factor brain maps
```

Key message:

```text
We compare brain-wide association maps, not SNP-level GWAS results.
```

### Slide 2: Background and Motivation

Content:

- Alzheimer disease is associated with multiple risk-factor domains.
- Many risk factors also show associations with brain structure.
- Whole-brain map-level comparison can reveal whether AD-related and
  risk-factor-related brain patterns are spatially similar.

Take-home message:

```text
The goal is to prioritize risk-factor domains whose brain-wide patterns overlap
with AD-related patterns.
```

### Slide 3: Data Sources

Content:

- AD BWAS maps from `AlzDisease_LMM/`.
- UKB risk-factor BWAS maps from `SSTAT_BingJing/`.
- 24 available AD-related files.
- 71 available UKB risk-factor files.
- Current clean analysis: 8 AD maps x 44 cleaned UKB maps.

Potential visual:

```text
AD BWAS maps  --->  brainMapR  <---  UKB risk-factor BWAS maps
```

### Slide 4: AD Maps Included

Content:

Table of 8 AD maps:

- AD vs HC
- MCI vs HC
- Conversion 1 year
- Conversion 2 years
- Conversion 3 years
- Conversion 4 years
- Conversion 5 years
- MMSE

Take-home message:

```text
The current analysis focuses on AD maps with confirmed sample sizes.
```

### Slide 5: Risk-Factor Domains

Content:

Broad risk-factor categories:

- Vascular / metabolic
- Psychiatric
- Lifestyle
- Social / socioeconomic
- Frailty / multimorbidity
- Infection / inflammatory
- Hormonal / sex-specific, where usable

Take-home message:

```text
The risk-factor set spans multiple AD-relevant biological and social domains.
```

### Slide 6: Method Overview

Content:

- Pairwise comparison of maps.
- One AD BWAS map x one risk-factor BWAS map.
- brainMapR estimates shared spatial signal / rGM / sumR2-style metrics.
- AVERAGE reference panel is primary.
- UKB reference panel is sensitivity.

Take-home message:

```text
Each cell in the final matrix is one AD map x one risk-factor map comparison.
```

### Slide 7: Workflow Progress

Content:

Timeline:

```text
Input audit -> pilot -> small batch -> exploratory full batch -> clean rerun
```

Status:

- Pilot completed.
- Exploratory full batch completed with identified failures.
- Clean AVERAGE rerun completed: 352 / 352.
- Clean UKB sensitivity rerun mostly completed: 344 / 352.

Take-home message:

```text
The core workflow is now operational, and the clean AVERAGE analysis completed.
```

### Slide 8: Clean Analysis Scope

Content:

```text
8 AD maps x 44 cleaned UKB risk-factor maps = 352 comparisons
```

Why cleaned:

- Excludes unusable upstream BWAS traits.
- Removes header-only files.
- Removes unstable or redundant trait definitions.
- Follows Baptiste's review comments.

Take-home message:

```text
The clean analysis is the primary interpretable result set.
```

### Slide 9: Main Output: rGM Heatmap

Content:

- Show vertical rGM heatmap.
- UKB traits on Y axis.
- AD maps on X axis.
- Highlight clusters or domains with stronger similarity.

Suggested figure:

```text
figure_1_vertical_rGM_heatmap.png
```

Take-home message:

```text
The heatmap summarizes which risk-factor maps are most similar to each AD map.
```

### Slide 10: Significance / P-Value Heatmap

Content:

- Show p-value heatmap.
- Mark Bonferroni / FWER-significant pairs.
- Discuss statistical evidence separately from effect direction / magnitude.

Suggested figure:

```text
figure_2_vertical_pvalue_heatmap.png
```

Take-home message:

```text
Significance helps prioritize pairs, but interpretation should also consider
stability and data quality.
```

### Slide 11: Top Stable Associations

Content:

- Show top stable rGM forest plot.
- Emphasize stable, interpretable associations.
- Avoid over-interpreting unstable or out-of-range estimates.

Suggested figure:

```text
figure_3_top_rGM_forest.png
```

Take-home message:

```text
Top associations should be selected from stable estimates, not only from
nominal significance.
```

### Slide 12: Reference Panel Sensitivity

Content:

- Compare AVERAGE vs UKB panel estimates.
- Show scatter plot.
- Identify pairs that are stable across reference panels.
- Flag pneumonia_ever issue under UKB panel.

Suggested figure:

```text
figure_4_reference_panel_sensitivity.png
```

Take-home message:

```text
Reference-panel sensitivity helps distinguish robust signals from panel- or
input-specific behavior.
```

### Slide 13: Data Quality Issues

Content:

- Header-only traits: Menopause, HRT ever used.
- Unusable traits from exploratory batch.
- Pneumonia UKB-panel failure.
- Out-of-range rGM concentrated in a small number of unstable traits.

Take-home message:

```text
Most workflow issues are trait-specific input quality issues rather than global
pipeline failures.
```

### Slide 14: Interpretation Boundaries

Content:

- Results indicate map-level similarity.
- They do not prove causality.
- They do not identify SNP-level mechanisms.
- They can guide follow-up analyses and trait prioritization.

Take-home message:

```text
The appropriate interpretation is shared spatial signal, not causal mechanism.
```

### Slide 15: Next Steps

Content:

- Review clean AVERAGE heatmap and top associations.
- Decide whether to exclude / flag unstable traits.
- Decide how to handle pneumonia_ever in UKB sensitivity analysis.
- Ask Elise / upstream phenotype team to confirm or regenerate problematic
  BWAS traits.
- Add or regenerate age at menopause if needed.
- Prepare final result matrices and figures for discussion.

Take-home message:

```text
The analysis is ready for scientific review, with a small set of upstream input
issues to resolve.
```

## 17. Suggested Spoken Narrative

Suggested narrative for presenting:

```text
This project compares Alzheimer-related BWAS maps with UKB-derived
risk-factor BWAS maps. The analysis is at the brain-map level. We are not
performing SNP-level GWAS or causal inference. Instead, we ask whether the
spatial pattern of association across the brain is similar between AD-related
traits and risk-factor traits.

We started with 24 AD-related BWAS files and 71 UKB risk-factor BWAS files.
For the current clean analysis, we use 8 AD maps with confirmed sample sizes
and 44 cleaned UKB risk-factor maps. The analysis produces a pairwise matrix of
8 by 44 comparisons.

The workflow has been validated through a pilot and an exploratory full batch.
After excluding unstable or unusable upstream traits, the clean AVERAGE
reference-panel rerun completed all 352 comparisons successfully. We also ran a
UKB reference-panel sensitivity analysis, where 344 of 352 comparisons
completed successfully; the remaining failures were all for pneumonia_ever.

The main outputs are rGM and p-value heatmaps, a top-association forest plot,
and a reference-panel sensitivity plot. The goal is to identify which
risk-factor-associated brain patterns show the strongest shared spatial signal
with AD-related maps, while avoiding causal over-interpretation.
```

## 18. Suggested Slide-Level Take-Home Messages

Use one take-home message per slide:

1. The project compares whole-brain association maps, not SNP-level results.
2. AD-related maps are compared against UKB risk-factor maps pairwise.
3. Current clean scope is 8 AD maps x 44 cleaned risk-factor maps.
4. AD sample sizes are fixed because meta-analysis files lack `NMISS`.
5. UKB sample sizes are read from `NMISS`.
6. AVERAGE is the main reference panel; UKB is sensitivity.
7. The clean AVERAGE run completed 352 / 352 comparisons.
8. The UKB sensitivity run completed 344 / 352 comparisons.
9. Remaining failures and unstable estimates are concentrated in specific
   traits.
10. Results should be interpreted as shared spatial signal, not causality.

## 19. Questions to Ask Baptiste / Collaborators

Useful discussion questions:

1. Should `pneumonia_ever` be excluded from the reference-panel sensitivity
   figure or kept with a warning?
2. Should `pa_vigorous_time`, `infect_recency_days`, and `pneumonia_ever` be
   excluded from the main interpretable rGM heatmap because of out-of-range or
   unstable estimates?
3. Should the clean figure focus on all 44 traits or only a smaller prioritized
   core set?
4. Should UKB reference-panel sensitivity be included in the main presentation
   or kept as supplementary QC?
5. Should Elise regenerate or check the upstream BWAS files for:
   - Menopause
   - HRT ever used
   - Age at menopause
   - Stroke
   - Periodontal disease
   - Pre-imaging hypertension
   - Medication count
6. How should out-of-range rGM estimates be reported?
7. Should results be ranked primarily by rGM magnitude, p-value, stability
   across reference panels, or a combined criterion?

## 20. Recommended Wording for Conclusions

Strong but cautious conclusion:

```text
The workflow now enables systematic map-level comparison between AD-related
BWAS maps and UKB-derived risk-factor BWAS maps. The clean AVERAGE analysis
completed all 352 planned comparisons, producing matrix-style results that can
be used to prioritize risk-factor domains showing shared spatial signal with
AD-related brain patterns.
```

If presenting preliminary results:

```text
Preliminary clean results suggest that multiple risk-factor domains can be
compared systematically against AD-related maps, but final biological
interpretation should focus on stable estimates and should exclude or flag
traits with known upstream data-quality issues.
```

Avoid saying:

```text
These risk factors cause AD.
```

Use instead:

```text
These risk-factor BWAS maps show map-level similarity or shared spatial signal
with AD-related BWAS maps.
```

