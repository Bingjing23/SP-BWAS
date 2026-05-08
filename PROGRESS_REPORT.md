# AD-BWAS Project Progress Report

## Project Aim

This project is investigating whether Alzheimer-related brain-wide association
maps share spatial signal with brain-wide association maps for
Alzheimer-relevant risk factors.

The core question is:

**Which risk-factor-associated brain patterns are most similar to
Alzheimer-related brain patterns?**

This is a brain-map-level comparison, not a SNP-level GWAS and not a causal
analysis.

## Data Prepared So Far

Two main groups of BWAS maps are being used.

### Alzheimer-Related BWAS Maps

There are 24 Alzheimer-related BWAS summary-statistic files available. These
include AD vs healthy controls, MCI vs healthy controls, conversion-to-AD
phenotypes, cognitive measures, and clinical measures.

At this stage, 8 Alzheimer-related maps have confirmed sample sizes and are
ready for the default analysis:

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

Other Alzheimer-related maps are available but are being held back until the
correct sample sizes are confirmed.

### Risk-Factor BWAS Maps

There are 71 UK Biobank-derived risk-factor BWAS maps available. These cover
several broad categories:

- Vascular and metabolic traits, such as hypertension, type 2 diabetes,
  cholesterol, LDL, BMI, and stroke
- Psychiatric traits, such as depression, anxiety, PTSD, bipolar disorder, and
  schizophrenia
- Lifestyle traits, such as smoking, alcohol use, sleep, and physical activity
- Social and socioeconomic traits, such as education, income, deprivation,
  loneliness, and social activity
- Frailty and multimorbidity traits
- Infection and inflammatory traits
- Hormonal and sex-specific traits

After excluding technical or non-analysis variables, 68 risk-factor maps are
included in the default analysis design.

## Pilot Completed

A first official pilot analysis has been completed successfully.

The pilot compared:

```text
AD vs HC map
x
UKB hypertension map
```

This confirmed that the main analysis strategy works for one Alzheimer map and
one risk-factor map.

The pilot also confirmed that the available data can be prepared in the format
required for the brain-map comparison, and that the chosen analysis workflow
can run successfully.

## Current Small-Batch Analysis

After the successful pilot, the project moved to a controlled small batch.

The small batch includes:

**Alzheimer-related maps**

- AD vs HC
- MCI vs HC

**Risk-factor maps**

- Hypertension
- Type 2 diabetes
- Hypercholesterolemia
- LDL direct
- Major depression
- Ever smoking

This gives:

```text
2 Alzheimer maps x 6 risk-factor maps = 12 pairwise comparisons
```

The small batch has already been submitted to the computing cluster.

The purpose of this small batch is to confirm that the workflow remains stable
beyond a single pilot comparison before scaling up to the full analysis. These
results will be used mainly for workflow validation and early signal checking,
not for final biological interpretation.

## Planned Full Analysis

If the small batch completes successfully, the next planned step is the full
default analysis:

```text
8 Alzheimer-related maps
x
68 UKB risk-factor maps
=
544 pairwise comparisons
```

The final output will be a matrix-style summary where:

```text
Rows = Alzheimer-related BWAS maps
Columns = UKB risk-factor BWAS maps
Cells = map-level similarity / shared spatial signal metrics
```

This will allow us to identify which risk-factor brain-wide association
patterns are most similar to Alzheimer-related brain-wide association patterns.

The full analysis grid can still be adjusted if Baptiste prefers a smaller or
different prioritized risk-factor set.

## Current Status

Current status:

```text
Data inventory: completed
Input readiness check: completed
Official pilot: completed successfully
Small batch: submitted to PBS / cluster
Full batch: pending small-batch validation
```

## Next Steps

1. Wait for the small batch to complete.
2. Check whether all 12 pairwise comparisons finished successfully.
3. Summarize the small-batch outputs.
4. If the small batch is successful, proceed to the full 544-comparison
   analysis.
5. After the full batch, organize the results into matrix-style summaries.
6. Interpret results cautiously as shared spatial signal or map-level
   similarity, not causal effects.
