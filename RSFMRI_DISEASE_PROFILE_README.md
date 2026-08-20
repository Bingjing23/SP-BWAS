# rsfMRI Disease Profile README

This document is the human-facing handoff for a separate rsfMRI disease-count
task that lives alongside the SP-BWAS downstream integration project.

## Why This Exists

The main SP-BWAS project compares AD-related BWAS maps against UKB risk-factor
BWAS maps with brainMapR. The rsfMRI disease profile task is different: it
counts disease-positive participants inside processed NeuroSTORM/UKB rsfMRI
cohorts, so the user can decide which diseases have enough usable subjects for
separate embedding or regression analyses.

This task should use Python because the upstream rsfMRI subject-list and
embedding workflows were built with Python-style file handling. R remains
appropriate for the brainMapR/BWAS workflow only.

## Cohorts

The cohort manifest is:

```text
manifests/rsfmri_cohorts.tsv
```

It defines two cohorts:

```text
primary:
  gwas_core_57485
  /working/lab_puyag/bingjinZ/UKBB/outputs/gwas_inputs_neurostorm_59053/max_subjects_core_covariates/subject_ids.txt

supplementary:
  embedding_59057
  /working/lab_puyag/bingjinZ/UKBB/outputs/neurostorm_embeddings_20227/neurostorm_embeddings_20227_8batch_59057.eids.txt
```

Use `gwas_core_57485` as the default main-analysis cohort because it is the
subject set already used for GWAS and phenotypic mapping and has core
covariates. Use `embedding_59057` as a supplementary coverage cohort because it
represents all subjects with formal eight-batch embeddings.

Do not merge diseases on raw UKB imaging case IDs such as
`1013356_20227_2_0`. Disease phenotype tables use participant-level `eid`. The
Python script normalizes case-like IDs to bare leading EID when reading usable
lists, but the preferred input is still one bare EID per line.

## Disease Groups

The disease definition table is:

```text
manifests/rsfmri_disease_definitions.tsv
```

Current default groups:

```text
neurodivergent:
  ASD, ADHD

neurodegenerative:
  AD, PD

mental_disorder:
  MDD, BD, SCZ, PTSD, phobAnx, anxDis
```

The user corrected the neurodegenerative group: it should include `AD`, not
`ADHD`.

The definition table currently maps each disease to candidate phenotype column
names. For `ASD`, `ADHD`, `AD`, and `PD`, the table assumes a confirmed wide
phenotype table will contain matching columns. For mental-disorder traits, the
columns match traits already seen in `SSTAT_BingJing/16_BingJing.Rmd` and
`SSTAT_BingJing/*.linear`.

## Instance and Disease Timing

UKB field `20227` has imaging instances for the original imaging visit and
first repeat imaging visit. The user processed all available 20227 instances
and then selected the first available instance per subject.

For disease counting, this selected first usable scan can be treated as the
index scan. If disease dates are available, the preferred disease status should
be anchored to that scan date:

```text
pre_scan_case      = first disease date <= selected scan date
post_scan_incident = first disease date > selected scan date
control            = no disease record
missing/uncertain  = missing disease/date information
```

The current Python script counts disease status from existing 0/1 phenotype
columns. If scan-date-aware counts are needed, prepare a phenotype table with
separate columns such as `AD_pre_scan`, `AD_post_scan_incident`, and update
`manifests/rsfmri_disease_definitions.tsv` accordingly.

## How To Run

Run from the project root on the system that can access `/working/lab_puyag/...`:

```bash
python3 scripts/08_summarize_rsfmri_disease_counts.py \
  --cohort-manifest manifests/rsfmri_cohorts.tsv \
  --phenotype-table /path/to/phenotype_table.tsv
```

The phenotype table must be a wide TSV or CSV with one participant EID column
and one or more disease columns. Recognized EID column names include:

```text
eid, EID, f.eid, IID, ID, participant_id, subject_id, V1
```

## Outputs

Outputs are written under:

```text
outputs/rsfmri_disease_profile/
```

Key files:

```text
rsfmri_disease_counts.tsv
rsfmri_group_counts.tsv
rsfmri_disease_missing_definitions.tsv
subjects_with_disease_in_gwas_core.tsv
subjects_with_disease_in_embedding_cohort.tsv
subjects_with_disease_all_cohorts.tsv
disease_positive_embedding_not_gwas_core.tsv
rsfmri_disease_summary_report.md
```

Use `rsfmri_disease_counts.tsv` for disease-level case/control/missing counts.
Use `rsfmri_group_counts.tsv` for group-level any-case counts. Use
`disease_positive_embedding_not_gwas_core.tsv` to identify disease-positive
subjects who have embeddings but are not in the GWAS/core covariate cohort.

## Current Blocker

The cohort lists live on the HPC path and are not visible from the current
desktop workspace. The actual run also needs a phenotype table path with the
confirmed disease definitions.

## Conversation Sync: 2026-08-17

The current window established these handoff decisions:

- Use Python for the rsfMRI disease-profile task. R should remain limited to
  the main SP-BWAS brainMapR workflow.
- The user previously downloaded and processed all available UKB field-20227
  instances, then selected the first available processed instance per subject.
- That first selected instance is an image-processing index scan, not by itself
  a disease-valid case/control definition.
- For disease analyses, disease status should ideally be anchored to the
  selected scan date. If diagnosis dates are available, separate pre-scan cases
  from post-scan incident diagnoses.
- Use the 57,485 GWAS/core covariate cohort as the default main-analysis
  cohort.
- Also report the 59,057 full embedding cohort as supplementary coverage.
- Compare disease-positive subjects present in the embedding cohort but absent
  from the GWAS/core covariate cohort.
- Do not use raw case IDs such as `1013356_20227_2_0` for disease merges. Use
  bare participant EID.
- The disease grouping correction is: neurodegenerative = `AD` and `PD`, not
  `PD` and `ADHD`.

The current implementation is:

```text
manifests/rsfmri_cohorts.tsv
manifests/rsfmri_disease_definitions.tsv
scripts/08_summarize_rsfmri_disease_counts.py
```

The current missing input is still:

```text
/path/to/phenotype_table.tsv
```

That phenotype table should contain `eid` and disease columns such as `ASD`,
`ADHD`, `AD`, `PD`, `majorDepression`, `BD`, `SCZ`, `PTSD`, `phobAnx`, and
`anxDis`, or scan-date-aware equivalents such as `AD_pre_scan`.
