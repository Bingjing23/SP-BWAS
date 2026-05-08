# SP-BWAS

This repository contains the code and analysis notes for an AD-related brain-wide association study (BWAS) workflow.

## What is here

- `SSTAT_BingJing/16_BingJing.Rmd`: R Markdown script for building UK Biobank-derived phenotypes and running BWAS/BLUP-style analyses with `OSCA`.
- `AlzDisease_LMM/`: AD BWAS result files and related outputs.

## Notes

- Data files are not included in version control.
- The scripts assume the original UK Biobank and REML/LDSC files already exist on the same server paths used in the Rmd.
- `brainMapR` is used downstream for summary-statistic based analyses.

## Quick start

1. Open [`SSTAT_BingJing/16_BingJing.Rmd`](SSTAT_BingJing/16_BingJing.Rmd).
2. Update the file paths if needed for your environment.
3. Run the phenotype construction and BWAS steps on the compute environment that has `OSCA` and the required UK Biobank inputs.

## Downstream brainMapR analysis

The intended downstream analysis compares Alzheimer BWAS summary statistics in
`AlzDisease_LMM/` against UKB-derived risk-factor BWAS summary statistics in
`SSTAT_BingJing/` using:

```r
brainMapR::sumR2_regression_bivariate(
  inputPath = c("AlzDisease_LMM/", "outputs/pilot/derived_inputs/"),
  bwasFile = c(
    "BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat",
    "sstat_FS_All_moda_total_hyperTension.Probe.linear"
  ),
  bwasSampleSize = c(3542, "NMISS"),
  refPanel = c("AVERAGE"),
  outputPath = "outputs/pilot/brainMapR_ADvsHC_hyperTension/"
)
```

The official pilot using the inputs below has run successfully:

- AD-side input: `AlzDisease_LMM/BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat`
- Trait-side input: `outputs/pilot/derived_inputs/sstat_FS_All_moda_total_hyperTension.Probe.linear`
- AD sample size: fixed numeric value `3542`.
- Trait sample size: use the `NMISS` column from the UKB trait file.
- Reference panel: `AVERAGE`.
- Output directory: `outputs/pilot/brainMapR_ADvsHC_hyperTension/`.
- Log directory: `logs/`.

Raw input files should not be modified. UKB trait files use `Voxel` as the
spatial identifier, while brainMapR expects `Probe`. Create derived corrected
copies under `outputs/pilot/derived_inputs/` and rename only the header column
`Voxel` to `Probe`.

AD meta-analysis files in `AlzDisease_LMM/` do not carry `NMISS` because the
sample size is not preserved in the meta-analysis summary statistics. Supply
fixed numeric sample sizes for the AD-side input, and keep using `"NMISS"` for
UKB BWAS files that contain an `NMISS` column. The currently confirmed AD
sample sizes are:

| AD trait | Sample size |
| --- | ---: |
| ADvsHC | 3542 |
| MCIvsHC | 3976 |
| Conversion1year | 1257 |
| Conversion2years | 1199 |
| Conversion3years | 1031 |
| Conversion4years | 1285 |
| Conversion5years | 1197 |
| MMSE | 6981 |

Server environment used for brainMapR:

```bash
module purge
module load R/4.5.0

export R_LIBS_USER=/mnt/backedup/home/bingjinZ/Rlibs/R-4.5.0
export RGL_USE_NULL=TRUE
export LD_LIBRARY_PATH=/software/conda-envs/pkgs/bzip2-1.0.8-h4bc722e_7/lib:${LD_LIBRARY_PATH:-}
```

Use the narrow `bzip2` library path above for `rgl` on the QIMR HPC. Do not use
the broader `/software/conda-envs/lib` path because it can cause
`magick`/ImageMagick dynamic-library conflicts.

Official package route:

```r
install.packages("devtools", lib = Sys.getenv("R_LIBS_USER"),
                 repos = "https://cloud.r-project.org")
devtools::install_github("jean997/GFA", force = TRUE)
devtools::install_github("baptisteCD/brainMapR", force = TRUE)
library(brainMapR)
library(GFA)
```

The validated package versions were:

- `brainMapR_1.1.0.9000`
- `GFA_1.0.0.0449`

Do not use the CRAN `GFA_1.0.5`, the older RDM
`brainMapR_0.8.0.9000.tar.gz`, or a local replacement for `GFA::ldsc_rg()` for
the official analysis.

Before running pilot or batch jobs, verify:

```bash
Rscript -e 'library(brainMapR); library(GFA); stopifnot("sumR2_regression_bivariate" %in% ls("package:brainMapR")); stopifnot("ldsc_rg" %in% ls("package:GFA")); stopifnot("snp_ldsc" %in% ls("package:GFA")); sessionInfo()'
```

Submit long brainMapR jobs through PBS on a compute node rather than running
long `Rscript` jobs on the login node.

## Analysis plan

The downstream analysis proceeds from a validated pilot to the manifest-driven
full batch. The pilot confirms the environment, package versions, input format,
sample-size handling, and reference-panel loading. The full batch then applies
the same workflow to the selected AD x UKB risk-factor grid.

The scientific comparison is map-level rather than SNP-level: each job compares
one AD-related BWAS map from `AlzDisease_LMM/` with one UKB risk-factor BWAS map
from `SSTAT_BingJing/` using
`brainMapR::sumR2_regression_bivariate()`. The target output is a pairwise
comparison table or matrix whose rows are AD-related BWAS maps and whose
columns are UKB risk-factor BWAS maps.

Current design decisions:

- Use the official `brainMapR_1.1.0.9000` plus `jean997/GFA` route only.
- Use `AVERAGE` as the primary reference panel for pilot and full batch.
  `UKB` can be added later as a sensitivity reference panel after the workflow
  is stable.
- Run `brainMapR::sumR2_regression_bivariate()` with
  `varConstrained = TRUE` explicitly set. This is the current BrainMapR default,
  but the batch wrapper passes it explicitly so variance-constrained rGM
  estimates are reproducible.
- Preserve all raw input files.
- Use fixed numeric sample sizes for AD meta-analysis BWAS files. The currently
  included AD sample sizes were confirmed by Baptiste because these
  meta-analysis summary statistics do not carry `NMISS`.
- Use the `NMISS` column for UKB risk-factor BWAS files.
- Convert UKB `Voxel` headers to `Probe` only in derived copies under
  `outputs/batch/derived_inputs/`.
- Some UKB files were written with an extra row-name column before `Chr`; the
  derived-input script removes that extra column when detected. Header-only
  UKB files are excluded from regenerated batch manifests.
- Include only AD traits with confirmed sample sizes in default batch designs.
- Exclude obvious non-analysis UKB variables such as `eid`, `sex`, pilot
  duplicate files, and header-only files from the default risk-factor batch.

The final submitted full-batch design was:

```text
8 AD maps with confirmed sample sizes x 66 UKB risk-factor maps = 528 jobs
```

Biological interpretation should focus on the collected full-batch outputs, not
on the pilot alone. Use cautious language such as shared spatial signal,
map-level similarity, and gray-matter correlation; do not interpret these
map-level BWAS comparisons as SNP-level GWAS or causal results.

### Known UKB input quirks

Three alcohol-related UKB BWAS files were observed to have 9 header fields but
10 fields in data rows because `write.table()` preserved row names:

```text
SSTAT_BingJing/sstat_FS_All_moda_total_Alcohol_merge_phenotype.linear
SSTAT_BingJing/sstat_FS_All_moda_total_Alcohol_week_avg.linear
SSTAT_BingJing/sstat_FS_All_moda_total_alcohol_intake_frequency.linear
```

Their raw rows start with an extra row-name field before the actual `Chr`
value. If this field is not removed, columns shift and `b` is read as `*`,
causing `brainMapR` to fail with an error such as
`non-numeric argument to binary operator`. Do not edit these raw files in place.
`scripts/02_fix_sumstats_headers.R` detects this pattern in derived copies and
removes the extra first field while also renaming `Voxel` to `Probe`.

Two hormonal/sex-specific UKB files were observed to contain headers but no
data rows and should be excluded unless regenerated:

```text
SSTAT_BingJing/sstat_FS_All_moda_total_Menopause.linear
SSTAT_BingJing/sstat_FS_All_moda_total_hrt_ever_used.linear
```

Regenerating manifests after these checks should exclude header-only files from
the default full batch.

## Batch workflow

The batch workflow is manifest-driven. Do not hand-code pairwise comparisons
one by one.

From the project root on the HPC:

```bash
module purge
module load R/4.5.0

export R_LIBS_USER=/mnt/backedup/home/bingjinZ/Rlibs/R-4.5.0
export RGL_USE_NULL=TRUE
export LD_LIBRARY_PATH=/software/conda-envs/pkgs/bzip2-1.0.8-h4bc722e_7/lib:${LD_LIBRARY_PATH:-}
```

First audit all AD and risk-factor BWAS inputs:

```bash
Rscript scripts/00_check_inputs.R
```

This writes:

- `outputs/input_audit/ad_bwas_inventory.tsv`
- `outputs/input_audit/risk_factor_bwas_inventory.tsv`
- `outputs/input_audit/input_readiness.tsv`
- `outputs/input_audit/blockers.tsv`

Then build the manifests and pairwise design tables:

```bash
Rscript scripts/01_make_manifests.R
```

This writes:

- `manifests/ad_bwas_manifest.tsv`
- `manifests/risk_factor_bwas_manifest.tsv`
- `manifests/brainmapr_pairwise_design.tsv`
- `manifests/brainmapr_clean_average_design.tsv`
- `manifests/brainmapr_clean_ukb_design.tsv`

The current default design includes the 8 AD traits with confirmed sample
sizes and excludes AD traits whose sample sizes are not yet confirmed. It also
excludes obvious non-analysis UKB columns such as `eid`, `sex`, pilot duplicate
files, and header-only files from the default risk-factor batch.

The clean design follows Baptiste's review comments. It removes specific
alcohol beverage variables, keeps global alcohol measures, keeps `T2D` as the
selected diabetes variable, keeps more clinical psychiatric diagnoses rather
than broad/pre-imaging definitions, and excludes currently unusable upstream
BWAS traits until they are regenerated or confirmed.

Before running brainMapR, create derived UKB inputs with `Probe` headers for
the clean AVERAGE and UKB designs:

```bash
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_clean_average_design.tsv \
  --force
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_clean_ukb_design.tsv \
  --force
```

Run a dry-run path check before submitting real jobs:

```bash
Rscript scripts/04_run_brainMapR_batch.R \
  --design manifests/brainmapr_clean_average_design.tsv \
  --dry-run
Rscript scripts/04_run_brainMapR_batch.R \
  --design manifests/brainmapr_clean_ukb_design.tsv \
  --dry-run
```

To run a single pair:

```bash
qsub -v DESIGN=manifests/brainmapr_clean_average_design.tsv,PAIR_ID=ADvsHC__hyperTension \
  scripts/run_brainMapR_batch.pbs
```

## Parallel clean batch

Run the clean grid as a PBS array rather than one serial job. Each array
task runs one row of
`manifests/brainmapr_clean_average_design.tsv` or
`manifests/brainmapr_clean_ukb_design.tsv`, so the scheduler can execute many
AD x risk-factor pairs in parallel.

Prepare all full-batch derived inputs first:

```bash
Rscript scripts/00_check_inputs.R
Rscript scripts/01_make_manifests.R
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_clean_average_design.tsv \
  --force
Rscript scripts/02_fix_sumstats_headers.R \
  --design manifests/brainmapr_clean_ukb_design.tsv \
  --force
```

Verify that the known row-name alcohol files were corrected in derived copies:

```bash
for f in \
  outputs/batch/derived_inputs/sstat_FS_All_moda_total_Alcohol_merge_phenotype.Probe.linear \
  outputs/batch/derived_inputs/sstat_FS_All_moda_total_Alcohol_week_avg.Probe.linear \
  outputs/batch/derived_inputs/sstat_FS_All_moda_total_alcohol_intake_frequency.Probe.linear
do
  echo "==== $f ===="
  awk 'NR == 1 {print "header_NF=" NF, $0}
       NR == 2 {print "row2_NF=" NF, $0}' "$f"
done
```

Each file should report `header_NF=9` and `row2_NF=9`, and row 2 should start
with the true `Chr` value rather than a row name.

Confirm the number of jobs:

```bash
average_n=$(($(wc -l < manifests/brainmapr_clean_average_design.tsv) - 1))
ukb_n=$(($(wc -l < manifests/brainmapr_clean_ukb_design.tsv) - 1))
echo "AVERAGE jobs: $average_n"
echo "UKB jobs: $ukb_n"
```

Submit the clean AVERAGE batch as an array:

```bash
qsub -J 1-${average_n} \
  -v DESIGN=manifests/brainmapr_clean_average_design.tsv \
  scripts/run_brainMapR_batch_array.pbs
```

Submit the clean UKB reference-panel sensitivity batch as a second array:

```bash
qsub -J 1-${ukb_n} \
  -v DESIGN=manifests/brainmapr_clean_ukb_design.tsv \
  scripts/run_brainMapR_batch_array.pbs
```

If the cluster limits the number of simultaneous array tasks, use the local
PBS syntax for throttling, for example `-J 1-${full_n}%40` if supported. If
array throttling is not supported, submit smaller ranges manually:

```bash
qsub -J 1-100 -v DESIGN=manifests/brainmapr_pairwise_design.tsv scripts/run_brainMapR_batch_array.pbs
qsub -J 101-200 -v DESIGN=manifests/brainmapr_pairwise_design.tsv scripts/run_brainMapR_batch_array.pbs
```

Monitor:

```bash
qstat -u $USER
ls -lh logs/brainMapR_array.*.R.log | tail
```

After the array finishes, check pair statuses:

```bash
find outputs/batch/brainMapR_pairs_clean_AVERAGE -name run_status.tsv \
  -exec awk -F '\t' 'NR > 1 {print $4}' {} \; | sort | uniq -c
find outputs/batch/brainMapR_pairs_clean_UKB -name run_status.tsv \
  -exec awk -F '\t' 'NR > 1 {print $4}' {} \; | sort | uniq -c
```

Any failed pair writes its own status file under
`outputs/batch/brainMapR_pairs/<pair_id>/run_status.tsv` and a compact failure
record under `outputs/batch/failures/`.

Collect full-batch outputs:

```bash
Rscript scripts/05_collect_brainMapR_outputs.R \
  --design manifests/brainmapr_clean_average_design.tsv \
  --output-dir outputs/batch/summary_clean_AVERAGE
Rscript scripts/05_collect_brainMapR_outputs.R \
  --design manifests/brainmapr_clean_ukb_design.tsv \
  --output-dir outputs/batch/summary_clean_UKB
```

Generate PNG/TIFF figures with readable labels, vertical UKB-trait layout,
Bonferroni/FWER markers, and optional AVERAGE-vs-UKB sensitivity scatter:

```bash
Rscript scripts/06_plot_brainMapR_summary.R \
  --summary-dir outputs/batch/summary_clean_AVERAGE \
  --design manifests/brainmapr_clean_average_design.tsv \
  --compare-summary-dir outputs/batch/summary_clean_UKB
```

## Full batch run summary

The full manifest-driven batch was run after successful pilot validation.

Input distribution after manifest generation:

```text
AD BWAS files audited: 24
AD maps included in default batch: 8
AD maps excluded because sample size is not yet confirmed: 16

UKB risk-factor BWAS files audited: 71
UKB risk-factor maps included in default batch: 66
UKB files excluded as non-analysis or unusable inputs: 5

Full submitted grid: 8 AD maps x 66 UKB risk-factor maps = 528 jobs
```

The 8 included AD maps are the traits with Baptiste-confirmed fixed sample
sizes:

```text
ADvsHC, MCIvsHC, Conversion1year, Conversion2years, Conversion3years,
Conversion4years, Conversion5years, MMSE
```

Two UKB risk-factor files were excluded because they contained a header but no
data rows:

```text
SSTAT_BingJing/sstat_FS_All_moda_total_Menopause.linear
SSTAT_BingJing/sstat_FS_All_moda_total_hrt_ever_used.linear
```

Three alcohol-related UKB files had an extra row-name column in the raw files.
These were retained after derived-input correction because they contain usable
statistics:

```text
SSTAT_BingJing/sstat_FS_All_moda_total_Alcohol_merge_phenotype.linear
SSTAT_BingJing/sstat_FS_All_moda_total_Alcohol_week_avg.linear
SSTAT_BingJing/sstat_FS_All_moda_total_alcohol_intake_frequency.linear
```

Full array outcome:

```text
Submitted jobs: 528
Successful jobs: 488
Failed jobs: 40
```

The 40 failures were not random. They correspond exactly to 5 UKB traits across
all 8 AD maps:

```text
htn_i10_preimg
med_20003_n_distinct_i2
periodontal_k05_preimg
psy_any_strict_preimg
stroke
```

Each failed with the same `brainMapR` internal error:

```text
object 'BWASsignif2' not found
```

Manual QC showed that these five derived trait files contain rows but no usable
summary statistics for regression:

```text
rows = 654002
usable = 0
p < 0.05 = 0
p < 5e-8 = 0
```

Inspection of the raw rows showed `b = 0`, `se = 0`, and `p = NA`, so these are
treated as unusable upstream BWAS outputs in the current batch rather than
random compute failures. They should be excluded from interpretation unless
regenerated upstream.

The successful 488 jobs were collected with
`scripts/05_collect_brainMapR_outputs.R` into `outputs/batch/summary/`.

Collected summary outputs:

```text
outputs/batch/summary/brainmapr_output_files.tsv
outputs/batch/summary/brainmapr_pairwise_results_long.tsv
outputs/batch/summary/matrix_rGM.tsv
outputs/batch/summary/matrix_pvalue_rGM.tsv
outputs/batch/summary/matrix_rGM_CI_lb.tsv
outputs/batch/summary/matrix_rGM_CI_ub.tsv
outputs/batch/summary/matrix_rGM_se.tsv
outputs/batch/summary/matrix_braincov.tsv
outputs/batch/summary/matrix_m2_1.tsv
outputs/batch/summary/matrix_m2_2.tsv
outputs/batch/summary/failed_pairs.tsv
outputs/batch/summary/failed_trait_error_summary.txt
outputs/batch/summary/excluded_unusable_traits.txt
```

For the reviewed clean rerun, generate PNG/TIFF figures and helper tables from
the collected matrices with:

```bash
Rscript scripts/06_plot_brainMapR_summary.R \
  --summary-dir outputs/batch/summary_clean_AVERAGE \
  --design manifests/brainmapr_clean_average_design.tsv \
  --compare-summary-dir outputs/batch/summary_clean_UKB
```

This writes:

```text
outputs/batch/summary_clean_AVERAGE/figures/figure_1_vertical_rGM_heatmap.png
outputs/batch/summary_clean_AVERAGE/figures/figure_1_vertical_rGM_heatmap.tiff
outputs/batch/summary_clean_AVERAGE/figures/figure_2_vertical_pvalue_heatmap.png
outputs/batch/summary_clean_AVERAGE/figures/figure_2_vertical_pvalue_heatmap.tiff
outputs/batch/summary_clean_AVERAGE/figures/figure_3_top_rGM_forest.png
outputs/batch/summary_clean_AVERAGE/figures/figure_3_top_rGM_forest.tiff
outputs/batch/summary_clean_AVERAGE/figures/figure_4_reference_panel_sensitivity.png
outputs/batch/summary_clean_AVERAGE/figures/figure_4_reference_panel_sensitivity.tiff
outputs/batch/summary_clean_AVERAGE/top_associations_bonferroni.tsv
outputs/batch/summary_clean_AVERAGE/out_of_range_rGM.tsv
outputs/batch/summary_clean_AVERAGE/reference_panel_rGM_comparison.tsv
outputs/batch/summary_clean_AVERAGE/figure_generation_report.txt
```

The earlier exploratory full-grid plotting summary was:

```text
AD maps: 8
Risk-factor maps with successful results: 61
Collected successful pairs: 488
FDR-significant pairs: 299
Out-of-range rGM estimates: 15
Top stable associations written: 50
```

The rGM heatmap clips colors to `[-1, 1]` for readability. A small number of
rGM estimates exceed the theoretical correlation range, for example values
greater than 1. These are flagged in `out_of_range_rGM.tsv` and marked with `!`
in the rGM heatmap. Treat them as unstable or boundary estimates rather than
interpretable biological effects.

## Main outputs

- Phenotype tables such as `df_*`
- BWAS result files under `Results/`
- AD summary statistics in `AlzDisease_LMM/`
