#!/usr/bin/env Rscript

# Build manifest tables and pairwise brainMapR design files.
#
# Purpose:
#   Create reproducible design tables for comparing AD BWAS maps against
#   UKB-derived risk-factor BWAS maps with brainMapR.
#
# Inputs:
#   AlzDisease_LMM/*.linear_random_oscaFormat
#   SSTAT_BingJing/*.linear
#
# Outputs:
#   manifests/ad_bwas_manifest.tsv
#   manifests/risk_factor_bwas_manifest.tsv
#   manifests/brainmapr_pairwise_design.tsv
#   manifests/brainmapr_small_batch_design.tsv
#
# How to run:
#   Rscript scripts/01_make_manifests.R

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  opts <- list(
    project_root = ".",
    output_dir = "manifests",
    derived_input_dir = "outputs/batch/derived_inputs",
    pair_output_root = "outputs/batch/brainMapR_pairs",
    reference_panel = "AVERAGE"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--project-root") {
      i <- i + 1
      opts$project_root <- args[[i]]
    } else if (key == "--output-dir") {
      i <- i + 1
      opts$output_dir <- args[[i]]
    } else if (key == "--derived-input-dir") {
      i <- i + 1
      opts$derived_input_dir <- args[[i]]
    } else if (key == "--pair-output-root") {
      i <- i + 1
      opts$pair_output_root <- args[[i]]
    } else if (key == "--reference-panel") {
      i <- i + 1
      opts$reference_panel <- args[[i]]
    } else if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/01_make_manifests.R [--project-root .]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 1
  }

  opts
}

normalize_path <- function(path) {
  gsub("\\\\", "/", path)
}

relative_path <- function(path, root) {
  path <- normalize_path(normalizePath(path, mustWork = FALSE))
  root <- normalize_path(normalizePath(root, mustWork = FALSE))
  sub(paste0("^", root, "/?"), "", path)
}

safe_read_header <- function(path) {
  header <- readLines(path, n = 1, warn = FALSE)
  scan(text = header, what = character(), quiet = TRUE)
}

has_data_rows <- function(path) {
  length(readLines(path, n = 2, warn = FALSE)) >= 2
}

infer_ad_id <- function(file_name) {
  known <- c(
    "BWAS_meta_AD_vs._HC_QC_8SD.linear_random_oscaFormat" = "ADvsHC",
    "BWAS_meta_MCI_vs._HC_QC_8SD.linear_random_oscaFormat" = "MCIvsHC",
    "BWAS_meta_Conversion_at_1_year_QC_8SD.linear_random_oscaFormat" = "Conversion1year",
    "BWAS_meta_Conversion_at_2_years_QC_8SD.linear_random_oscaFormat" = "Conversion2years",
    "BWAS_meta_Conversion_at_3_years_QC_8SD.linear_random_oscaFormat" = "Conversion3years",
    "BWAS_meta_Conversion_at_4_years_QC_8SD.linear_random_oscaFormat" = "Conversion4years",
    "BWAS_meta_Conversion_at_5_years_QC_8SD.linear_random_oscaFormat" = "Conversion5years",
    "BWAS_meta_Mini_Mental_State_Examination_QC_8SD.linear_random_oscaFormat" = "MMSE"
  )
  if (file_name %in% names(known)) {
    return(unname(known[[file_name]]))
  }

  id <- file_name
  id <- sub("^BWAS_meta_", "", id)
  id <- sub("_QC_8SD\\.linear_random_oscaFormat$", "", id)
  id <- gsub("%", "percent", id)
  id <- gsub("vs\\._", "vs", id)
  id <- gsub("_-_", "_", id)
  id <- gsub("_+", "_", id)
  id <- gsub("[^A-Za-z0-9]+", "_", id)
  gsub("^_|_$", "", id)
}

confirmed_ad_sample_sizes <- function() {
  c(
    ADvsHC = 3542,
    MCIvsHC = 3976,
    Conversion1year = 1257,
    Conversion2years = 1199,
    Conversion3years = 1031,
    Conversion4years = 1285,
    Conversion5years = 1197,
    MMSE = 6981
  )
}

infer_trait_id <- function(file_name) {
  id <- sub("^sstat_FS_All_moda_total_", "", file_name)
  id <- sub("\\.linear$", "", id)
  gsub("[^A-Za-z0-9]+", "_", id)
}

infer_category <- function(trait_id) {
  id <- tolower(trait_id)
  if (id %in% c("eid", "sex") || grepl("_pilot$", id)) {
    return("exclude_candidate")
  }
  if (grepl("hypert|htn|t2d|diabetes|cholesterol|ldl|bmi|stroke", id)) {
    return("vascular_metabolic")
  }
  if (grepl("depression|mdd|anx|ptsd|bd|scz|phob|psy_", id)) {
    return("psychiatric")
  }
  if (grepl("smoking|pack_year|sleep|insomnia|napping|alcohol|drinker|wine|beer|spirits|pa_", id)) {
    return("lifestyle")
  }
  if (grepl("education|income|imd|social|loneliness|isolation", id)) {
    return("social_socioeconomic")
  }
  if (grepl("frailty|multimorbidity|med_20003", id)) {
    return("frailty_multimorbidity")
  }
  if (grepl("menopause|hrt", id)) {
    return("hormonal_sex_specific")
  }
  if (grepl("infect|lrti|pneumonia|sepsis|ssti|uti|periodontal", id)) {
    return("infection_inflammatory")
  }
  "other_candidate"
}

probe_file_name <- function(file_name) {
  sub("\\.linear$", ".Probe.linear", file_name)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(opts$project_root, mustWork = TRUE)
manifest_dir <- file.path(root, opts$output_dir)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

ad_files <- sort(Sys.glob(file.path(root, "AlzDisease_LMM", "*.linear_random_oscaFormat")))
risk_files <- sort(Sys.glob(file.path(root, "SSTAT_BingJing", "*.linear")))
sample_sizes <- confirmed_ad_sample_sizes()

ad_manifest <- do.call(rbind, lapply(ad_files, function(path) {
  file_name <- basename(path)
  ad_id <- infer_ad_id(file_name)
  columns <- safe_read_header(path)
  sample_size <- if (ad_id %in% names(sample_sizes)) {
    unname(sample_sizes[[ad_id]])
  } else {
    NA_real_
  }
  has_sample_size <- !is.null(sample_size) && !is.na(sample_size)
  has_data <- has_data_rows(path)
  data.frame(
    ad_id = ad_id,
    phenotype = ad_id,
    ad_file = relative_path(path, root),
    file_name = file_name,
    sample_size = if (has_sample_size) sample_size else NA_real_,
    sample_size_source = if (has_sample_size) "fixed_numeric_confirmed" else "missing_unconfirmed",
    include = has_sample_size && has_data && "Probe" %in% columns,
    notes = if (!has_data) {
      "no data rows"
    } else if (has_sample_size) {
      "confirmed by Baptiste"
    } else {
      "sample size not yet confirmed"
    },
    stringsAsFactors = FALSE
  )
}))

risk_manifest <- do.call(rbind, lapply(risk_files, function(path) {
  file_name <- basename(path)
  trait_id <- infer_trait_id(file_name)
  columns <- safe_read_header(path)
  has_probe <- "Probe" %in% columns
  has_voxel <- "Voxel" %in% columns
  derived_file <- if (!has_probe && has_voxel) {
    file.path(opts$derived_input_dir, probe_file_name(file_name))
  } else {
    relative_path(path, root)
  }
  category <- infer_category(trait_id)
  has_data <- has_data_rows(path)
  include <- category != "exclude_candidate" &&
    has_data &&
    "NMISS" %in% columns &&
    (has_probe || has_voxel)
  data.frame(
    trait_id = trait_id,
    trait = trait_id,
    category = category,
    raw_file = relative_path(path, root),
    file_name = file_name,
    analysis_file = derived_file,
    sample_size_source = if ("NMISS" %in% columns) "NMISS" else "missing",
    has_Probe = has_probe,
    has_Voxel = has_voxel,
    needs_Voxel_to_Probe = !has_probe && has_voxel,
    include = include,
    notes = if (!has_data) {
      "no data rows"
    } else if (include) {
      "candidate risk-factor BWAS"
    } else {
      "excluded from default batch"
    },
    stringsAsFactors = FALSE
  )
}))

included_ad <- ad_manifest[ad_manifest$include, , drop = FALSE]
included_risk <- risk_manifest[risk_manifest$include, , drop = FALSE]

make_design <- function(ad_df, risk_df) {
  if (nrow(ad_df) == 0 || nrow(risk_df) == 0) {
    return(data.frame())
  }

  rows <- vector("list", nrow(ad_df) * nrow(risk_df))
  idx <- 1
  for (i in seq_len(nrow(ad_df))) {
    for (j in seq_len(nrow(risk_df))) {
      ad <- ad_df[i, ]
      risk <- risk_df[j, ]
      pair_id <- paste(ad$ad_id, risk$trait_id, sep = "__")
      rows[[idx]] <- data.frame(
        pair_id = pair_id,
        ad_id = ad$ad_id,
        ad_file = ad$ad_file,
        ad_input_path = dirname(ad$ad_file),
        ad_bwas_file = basename(ad$ad_file),
        ad_sample_size = ad$sample_size,
        trait_id = risk$trait_id,
        trait_file = risk$analysis_file,
        trait_category = risk$category,
        trait_input_path = dirname(risk$analysis_file),
        trait_bwas_file = basename(risk$analysis_file),
        trait_sample_size = risk$sample_size_source,
        reference_panel = opts$reference_panel,
        output_dir = file.path(opts$pair_output_root, pair_id),
        include = TRUE,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  do.call(rbind, rows)
}

full_design <- make_design(included_ad, included_risk)

small_ad_ids <- c("ADvsHC", "MCIvsHC")
small_trait_ids <- c(
  "hyperTension",
  "T2D",
  "hyperCholesterolemia",
  "ldl_direct",
  "majorDepression",
  "smoking_ever"
)
small_design <- make_design(
  included_ad[included_ad$ad_id %in% small_ad_ids, , drop = FALSE],
  included_risk[included_risk$trait_id %in% small_trait_ids, , drop = FALSE]
)

write.table(ad_manifest, file.path(manifest_dir, "ad_bwas_manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(risk_manifest, file.path(manifest_dir, "risk_factor_bwas_manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_design, file.path(manifest_dir, "brainmapr_pairwise_design.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(small_design, file.path(manifest_dir, "brainmapr_small_batch_design.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Manifest generation complete.\n")
cat("Included AD maps:", nrow(included_ad), "\n")
cat("Included risk-factor maps:", nrow(included_risk), "\n")
cat("Full pairwise jobs:", nrow(full_design), "\n")
cat("Small-batch jobs:", nrow(small_design), "\n")
cat("Output directory:", relative_path(manifest_dir, root), "\n")
