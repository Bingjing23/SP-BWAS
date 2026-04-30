#!/usr/bin/env Rscript

# Audit AD and UKB BWAS summary-statistic inputs for brainMapR.
#
# Purpose:
#   Inspect AlzDisease_LMM/*.linear_random_oscaFormat and
#   SSTAT_BingJing/*.linear before batch-running
#   brainMapR::sumR2_regression_bivariate().
#
# Inputs:
#   AlzDisease_LMM/
#   SSTAT_BingJing/
#
# Outputs:
#   outputs/input_audit/ad_bwas_inventory.tsv
#   outputs/input_audit/risk_factor_bwas_inventory.tsv
#   outputs/input_audit/input_readiness.tsv
#   outputs/input_audit/blockers.tsv
#
# How to run:
#   Rscript scripts/00_check_inputs.R

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  opts <- list(
    project_root = ".",
    output_dir = "outputs/input_audit"
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
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/00_check_inputs.R [--project-root .] [--output-dir outputs/input_audit]\n"
      )
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
  if (!file.exists(path)) {
    return(list(readable = FALSE, columns = character(), issue = "file_missing"))
  }
  if (file.info(path)$size == 0) {
    return(list(readable = FALSE, columns = character(), issue = "file_empty"))
  }

  header <- tryCatch(
    readLines(path, n = 1, warn = FALSE),
    error = function(e) NA_character_
  )
  if (length(header) == 0 || is.na(header) || !nzchar(header)) {
    return(list(readable = FALSE, columns = character(), issue = "header_unreadable"))
  }

  columns <- tryCatch(
    scan(text = header, what = character(), quiet = TRUE),
    error = function(e) character()
  )
  if (length(columns) == 0) {
    return(list(readable = FALSE, columns = character(), issue = "header_parse_failed"))
  }

  list(readable = TRUE, columns = columns, issue = "")
}

has_any <- function(columns, candidates) {
  any(tolower(columns) %in% tolower(candidates))
}

first_existing <- function(columns, candidates) {
  idx <- match(tolower(candidates), tolower(columns), nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) {
    return("")
  }
  columns[[idx[[1]]]]
}

infer_ad_id <- function(file_name) {
  id <- file_name
  id <- sub("^BWAS_meta_", "", id)
  id <- sub("_QC_8SD\\.linear_random_oscaFormat$", "", id)
  id <- gsub("%", "percent", id)
  id <- gsub("vs\\._", "vs", id)
  id <- gsub("_-_", "_", id)
  id <- gsub("_+", "_", id)
  id <- gsub("[^A-Za-z0-9]+", "_", id)
  id <- gsub("^_|_$", "", id)
  id
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

readiness_row <- function(path, root, group) {
  header <- safe_read_header(path)
  columns <- header$columns

  has_probe <- "Probe" %in% columns
  has_voxel <- "Voxel" %in% columns
  has_nmiss <- "NMISS" %in% columns
  has_beta <- has_any(columns, c("b", "beta", "BETA"))
  has_se <- has_any(columns, c("se", "SE"))
  has_p <- has_any(columns, c("p", "P", "p_value", "PVAL"))

  needs_voxel_to_probe <- !has_probe && has_voxel
  issue <- header$issue
  if (header$readable) {
    missing <- character()
    if (!has_probe && !has_voxel) missing <- c(missing, "missing Probe/Voxel")
    if (!has_beta) missing <- c(missing, "missing beta")
    if (!has_se) missing <- c(missing, "missing se")
    if (!has_p) missing <- c(missing, "missing p")
    if (group == "risk" && !has_nmiss) missing <- c(missing, "missing NMISS")
    if (length(missing) > 0) {
      issue <- paste(missing, collapse = "; ")
    } else if (needs_voxel_to_probe) {
      issue <- "needs derived Probe copy"
    }
  }

  ready <- header$readable &&
    has_probe &&
    has_beta &&
    has_se &&
    has_p &&
    (group != "risk" || has_nmiss)

  data.frame(
    file = relative_path(path, root),
    has_Probe = has_probe,
    has_Voxel = has_voxel,
    needs_Voxel_to_Probe = needs_voxel_to_probe,
    has_NMISS = has_nmiss,
    has_beta = has_beta,
    has_se = has_se,
    has_p = has_p,
    readable = header$readable,
    ready = ready,
    issue = issue,
    stringsAsFactors = FALSE
  )
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(opts$project_root, mustWork = TRUE)
out_dir <- file.path(root, opts$output_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ad_files <- sort(Sys.glob(file.path(root, "AlzDisease_LMM", "*.linear_random_oscaFormat")))
risk_files <- sort(Sys.glob(file.path(root, "SSTAT_BingJing", "*.linear")))

if (length(ad_files) == 0) {
  stop("No AD BWAS files found under AlzDisease_LMM/")
}
if (length(risk_files) == 0) {
  stop("No risk-factor BWAS files found under SSTAT_BingJing/")
}

readiness <- do.call(
  rbind,
  c(
    lapply(ad_files, readiness_row, root = root, group = "ad"),
    lapply(risk_files, readiness_row, root = root, group = "risk")
  )
)

ad_inventory <- do.call(rbind, lapply(ad_files, function(path) {
  header <- safe_read_header(path)
  columns <- header$columns
  file_name <- basename(path)
  sample_source <- if ("NMISS" %in% columns) "NMISS" else "fixed_numeric_required"
  ready <- header$readable &&
    "Probe" %in% columns &&
    has_any(columns, c("b", "beta", "BETA")) &&
    has_any(columns, c("se", "SE")) &&
    has_any(columns, c("p", "P", "p_value", "PVAL"))
  data.frame(
    phenotype = infer_ad_id(file_name),
    file_path = relative_path(path, root),
    file_name = file_name,
    format = "linear_random_oscaFormat",
    first_columns = paste(head(columns, 12), collapse = ","),
    sample_size_source = sample_source,
    ready_for_brainMapR = ready,
    notes = if (ready) "AD meta-analysis requires fixed numeric sample size" else header$issue,
    stringsAsFactors = FALSE
  )
}))

risk_inventory <- do.call(rbind, lapply(risk_files, function(path) {
  header <- safe_read_header(path)
  columns <- header$columns
  file_name <- basename(path)
  trait <- infer_trait_id(file_name)
  has_probe <- "Probe" %in% columns
  has_voxel <- "Voxel" %in% columns
  ready <- header$readable &&
    has_probe &&
    "NMISS" %in% columns &&
    has_any(columns, c("b", "beta", "BETA")) &&
    has_any(columns, c("se", "SE")) &&
    has_any(columns, c("p", "P", "p_value", "PVAL"))
  notes <- if (ready) {
    "ready"
  } else if (!has_probe && has_voxel) {
    "requires derived Probe copy from Voxel header"
  } else {
    header$issue
  }
  data.frame(
    trait = trait,
    category = infer_category(trait),
    file_path = relative_path(path, root),
    file_name = file_name,
    format = "linear",
    first_columns = paste(head(columns, 12), collapse = ","),
    sample_size_source = if ("NMISS" %in% columns) "NMISS" else "missing",
    ready_for_brainMapR = ready,
    notes = notes,
    stringsAsFactors = FALSE
  )
}))

blockers <- list()
missing_or_unreadable <- readiness$file[!readiness$readable]
if (length(missing_or_unreadable) > 0) {
  blockers <- c(blockers, list(data.frame(
    blocker = "Unreadable or empty input files",
    affected_files = paste(missing_or_unreadable, collapse = ";"),
    severity = "high",
    proposed_fix = "Restore or replace the affected input files before batch analysis",
    stringsAsFactors = FALSE
  )))
}

missing_required <- readiness$file[readiness$readable & !readiness$has_Probe & !readiness$has_Voxel]
if (length(missing_required) > 0) {
  blockers <- c(blockers, list(data.frame(
    blocker = "Missing spatial identifier column",
    affected_files = paste(missing_required, collapse = ";"),
    severity = "high",
    proposed_fix = "Confirm the correct spatial identifier column; brainMapR expects Probe",
    stringsAsFactors = FALSE
  )))
}

needs_probe <- readiness$file[readiness$needs_Voxel_to_Probe]
if (length(needs_probe) > 0) {
  blockers <- c(blockers, list(data.frame(
    blocker = "Voxel column must be renamed to Probe in derived copies",
    affected_files = paste(needs_probe, collapse = ";"),
    severity = "medium",
    proposed_fix = "Run scripts/02_fix_sumstats_headers.R to create outputs/batch/derived_inputs/*.Probe.linear",
    stringsAsFactors = FALSE
  )))
}

missing_nmiss <- readiness$file[grepl("^SSTAT_BingJing/", readiness$file) & !readiness$has_NMISS]
if (length(missing_nmiss) > 0) {
  blockers <- c(blockers, list(data.frame(
    blocker = "Missing NMISS in UKB risk-factor BWAS",
    affected_files = paste(missing_nmiss, collapse = ";"),
    severity = "high",
    proposed_fix = "Confirm a fixed sample size or regenerate the UKB BWAS with NMISS",
    stringsAsFactors = FALSE
  )))
}

if (length(blockers) == 0) {
  blockers_df <- data.frame(
    blocker = "No high-severity input-format blockers detected",
    affected_files = "",
    severity = "none",
    proposed_fix = "Proceed to manifest generation and derived Probe copies",
    stringsAsFactors = FALSE
  )
} else {
  blockers_df <- do.call(rbind, blockers)
}

write.table(ad_inventory, file.path(out_dir, "ad_bwas_inventory.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(risk_inventory, file.path(out_dir, "risk_factor_bwas_inventory.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(readiness, file.path(out_dir, "input_readiness.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(blockers_df, file.path(out_dir, "blockers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Input audit complete.\n")
cat("AD files:", length(ad_files), "\n")
cat("Risk-factor files:", length(risk_files), "\n")
cat("Output directory:", relative_path(out_dir, root), "\n")
