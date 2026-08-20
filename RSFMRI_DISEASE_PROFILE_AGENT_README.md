# rsfMRI Disease Profile Agent Handoff

This is the AI-agent-facing handoff for the rsfMRI disease-count task. It
should be read with `AGENTS.md`, but this document narrows scope for the
separate NeuroSTORM/UKB rsfMRI cohort-count work.

## Scope Boundary

This task is not the brainMapR BWAS map-level integration. It is a separate
subject-level disease-count and cohort-comparison task used by another project.

Use Python for this task. The user explicitly noted that upstream rsfMRI and
embedding subject filtering was done with Python. Do not add R scripts for this
rsfMRI disease-count module unless the user asks.

Do not reorganize the repository. Keep new files in:

```text
manifests/
scripts/
outputs/rsfmri_disease_profile/
```

## Files Added For This Task

```text
manifests/rsfmri_cohorts.tsv
manifests/rsfmri_disease_definitions.tsv
scripts/08_summarize_rsfmri_disease_counts.py
RSFMRI_DISEASE_PROFILE_README.md
RSFMRI_DISEASE_PROFILE_AGENT_README.md
```

## Cohort Rules

Primary cohort:

```text
cohort_id: gwas_core_57485
path: /working/lab_puyag/bingjinZ/UKBB/outputs/gwas_inputs_neurostorm_59053/max_subjects_core_covariates/subject_ids.txt
role: default main analysis
```

Supplementary cohort:

```text
cohort_id: embedding_59057
path: /working/lab_puyag/bingjinZ/UKBB/outputs/neurostorm_embeddings_20227/neurostorm_embeddings_20227_8batch_59057.eids.txt
role: embedding coverage analysis
```

Never use raw 20227 case IDs as disease merge keys. Convert entries such as
`1013356_20227_2_0` to bare `1013356`. The Python script already strips
leading digits for usable-list IDs.

## Disease Groups

Use the current disease definition manifest unless the user provides a better
phenotype table or disease-date-aware definitions:

```text
manifests/rsfmri_disease_definitions.tsv
```

Current grouping:

```text
neurodivergent: ASD, ADHD
neurodegenerative: AD, PD
mental_disorder: MDD, BD, SCZ, PTSD, phobAnx, anxDis
```

Important correction from the user: neurodegenerative includes `AD`, not
`ADHD`.

## Required Input Still Missing

The script cannot produce real counts until the user supplies a phenotype table
or a path to one. It must be a wide table with EID plus disease columns, or an
already scan-date-anchored disease table.

If the user supplies only UKB raw `.tab` plus disease definitions, build a
separate phenotype table first; do not silently infer disease definitions from
filenames.

## Preferred Disease Timing

The user processed all field-20227 instances and selected the first usable
instance per subject. For disease analysis this should be treated as the index
scan only if disease status is anchored correctly.

If first diagnosis dates are available, prefer:

```text
pre_scan_case      = first disease date <= selected scan date
post_scan_incident = first disease date > selected scan date
control            = no disease record
missing/uncertain  = missing disease/date information
```

The current script counts existing 0/1 disease columns. It does not calculate
pre-scan/post-scan status unless those columns already exist in the phenotype
table.

## Command

Run on the HPC or another environment that can access `/working/lab_puyag/...`:

```bash
python3 scripts/08_summarize_rsfmri_disease_counts.py \
  --cohort-manifest manifests/rsfmri_cohorts.tsv \
  --phenotype-table /path/to/phenotype_table.tsv
```

Expected outputs:

```text
outputs/rsfmri_disease_profile/rsfmri_disease_counts.tsv
outputs/rsfmri_disease_profile/rsfmri_group_counts.tsv
outputs/rsfmri_disease_profile/subjects_with_disease_in_gwas_core.tsv
outputs/rsfmri_disease_profile/subjects_with_disease_in_embedding_cohort.tsv
outputs/rsfmri_disease_profile/disease_positive_embedding_not_gwas_core.tsv
```

## Existing Project Summary

The rest of the repository has already done or documented:

1. UKB `.tab` phenotype extraction and trait construction in
   `SSTAT_BingJing/16_BingJing.Rmd`.
2. OSCA BWAS runs producing `SSTAT_BingJing/sstat_FS_All_moda_total_*.linear`.
3. AD-related BWAS summary-stat ingestion from `AlzDisease_LMM/`.
4. Input audit and manifest generation for brainMapR.
5. Derived UKB BWAS header correction from `Voxel` to `Probe`.
6. brainMapR pilot using ADvsHC x hypertension.
7. Manifest-driven AD x UKB risk-factor batch analysis.
8. Collection of pairwise brainMapR outputs into long tables and matrices.
9. Figure generation and supervisor-facing summaries.
10. New rsfMRI disease-profile framework described in this document.

## Current Window Sync: 2026-08-17

If a future agent resumes from this file, preserve the following conclusions:

1. The user asked why R appeared in the rsfMRI disease-count task. The correct
   resolution is that R belongs to the brainMapR workflow, while the rsfMRI
   disease-count task should use Python.
2. The R implementation for `scripts/08_summarize_rsfmri_disease_counts.*` was
   replaced with `scripts/08_summarize_rsfmri_disease_counts.py`.
3. The Python script was tested with a minimal synthetic cohort and phenotype
   table under `/private/tmp`; it produced disease counts and correctly
   normalized case-like IDs such as `1001_20227_2_0` to bare EID `1001`.
4. The user's two cohort paths are HPC paths and are not readable from the
   current desktop workspace. Do not conclude they are invalid; run on the HPC
   or another mounted environment.
5. The task is not complete until a phenotype table path is supplied. The
   phenotype table must contain disease columns, or precomputed scan-date-aware
   disease columns if disease timing is required.
6. If the user later supplies a raw UKB `.tab` table instead of a phenotype
   table, build an explicit phenotype table first and document every disease
   definition. Do not infer disease definitions from file names alone.
7. When reporting results, provide both per-disease counts and group-level
   any-case counts for the primary `gwas_core_57485` cohort and supplementary
   `embedding_59057` cohort.
8. Also report `disease_positive_embedding_not_gwas_core.tsv` so the user can
   see which disease-positive subjects have embeddings but are not in the
   GWAS/core covariate cohort.

When modifying documentation, keep this module's scope separate from the
brainMapR BWAS workflow. Do not move existing directories or raw data.
