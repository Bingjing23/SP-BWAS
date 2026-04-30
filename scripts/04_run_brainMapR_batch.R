#!/usr/bin/env Rscript

# Run manifest-driven brainMapR bivariate SumR2 analyses.
#
# Purpose:
#   Execute one or more rows from a pairwise design table using the official
#   brainMapR + GFA route.
#
# Inputs:
#   manifests/brainmapr_small_batch_design.tsv or
#   manifests/brainmapr_pairwise_design.tsv
#
# Outputs:
#   outputs/batch/brainMapR_pairs/<pair_id>/
#
# How to run:
#   Rscript scripts/04_run_brainMapR_batch.R \
#     --design manifests/brainmapr_small_batch_design.tsv \
#     --pair-id ADvsHC__hyperTension

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  opts <- list(
    project_root = ".",
    design = "manifests/brainmapr_small_batch_design.tsv",
    pair_id = "",
    max_pairs = NA_integer_,
    dry_run = FALSE,
    continue_on_error = FALSE
  )

  for (arg in args) {
    if (arg == "--dry-run") opts$dry_run <- TRUE
    if (arg == "--continue-on-error") opts$continue_on_error <- TRUE
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--project-root") {
      i <- i + 1
      opts$project_root <- args[[i]]
    } else if (key == "--design") {
      i <- i + 1
      opts$design <- args[[i]]
    } else if (key == "--pair-id") {
      i <- i + 1
      opts$pair_id <- args[[i]]
    } else if (key == "--max-pairs") {
      i <- i + 1
      opts$max_pairs <- as.integer(args[[i]])
    } else if (key %in% c("--dry-run", "--continue-on-error")) {
      # already handled
    } else if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/04_run_brainMapR_batch.R [--design <tsv>] [--pair-id <id>] [--dry-run]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 1
  }

  opts
}

check_package_functions <- function() {
  suppressPackageStartupMessages({
    library(brainMapR)
    library(GFA)
  })

  if (!"sumR2_regression_bivariate" %in% ls("package:brainMapR")) {
    stop("brainMapR::sumR2_regression_bivariate() is not available")
  }
  if (!"ldsc_rg" %in% ls("package:GFA")) {
    stop("GFA::ldsc_rg() is not available")
  }
  if (!"snp_ldsc" %in% ls("package:GFA")) {
    stop("GFA::snp_ldsc() is not available")
  }

  cat("brainMapR:", as.character(packageVersion("brainMapR")), "\n")
  cat("GFA:", as.character(packageVersion("GFA")), "\n")
}

write_status <- function(output_dir, row, status, error = "") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status_df <- data.frame(
    pair_id = row$pair_id,
    ad_id = row$ad_id,
    trait_id = row$trait_id,
    status = status,
    error = error,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  write.table(status_df, file.path(output_dir, "run_status.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

validate_pair_inputs <- function(root, row) {
  ad_path <- file.path(root, row$ad_input_path, row$ad_bwas_file)
  trait_path <- file.path(root, row$trait_input_path, row$trait_bwas_file)
  if (!file.exists(ad_path)) {
    stop("AD input does not exist: ", ad_path)
  }
  if (!file.exists(trait_path)) {
    stop("Trait input does not exist: ", trait_path)
  }
  invisible(TRUE)
}

run_pair <- function(root, row, dry_run = FALSE) {
  output_dir <- file.path(root, row$output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  validate_pair_inputs(root, row)

  metadata <- data.frame(
    key = c(
      "pair_id", "ad_id", "trait_id", "ad_file", "trait_file",
      "ad_sample_size", "trait_sample_size", "reference_panel"
    ),
    value = c(
      row$pair_id, row$ad_id, row$trait_id, row$ad_file, row$trait_file,
      row$ad_sample_size, row$trait_sample_size, row$reference_panel
    ),
    stringsAsFactors = FALSE
  )
  write.table(metadata, file.path(output_dir, "run_metadata.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat("Pair:", row$pair_id, "\n")
  cat("AD:", row$ad_file, "\n")
  cat("Trait:", row$trait_file, "\n")
  cat("Output:", row$output_dir, "\n")

  if (dry_run) {
    write_status(output_dir, row, "dry_run")
    return(invisible(NULL))
  }

  result <- brainMapR::sumR2_regression_bivariate(
    inputPath = c(
      paste0(row$ad_input_path, "/"),
      paste0(row$trait_input_path, "/")
    ),
    bwasFile = c(row$ad_bwas_file, row$trait_bwas_file),
    bwasSampleSize = c(as.character(row$ad_sample_size), row$trait_sample_size),
    refPanel = c(row$reference_panel),
    outputPath = row$output_dir
  )

  saveRDS(result, file.path(output_dir, "sumR2_regression_bivariate_result.rds"))
  write_status(output_dir, row, "success")
  print(result)
  invisible(result)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(opts$project_root, mustWork = TRUE)
setwd(root)

design_path <- file.path(root, opts$design)
if (!file.exists(design_path)) {
  stop("Design table does not exist: ", design_path)
}

design <- read.delim(design_path, check.names = FALSE)
if (!"pair_id" %in% names(design)) {
  stop("Design table is missing pair_id column: ", design_path)
}
if (nzchar(opts$pair_id)) {
  design <- design[design$pair_id == opts$pair_id, , drop = FALSE]
}
if (!is.na(opts$max_pairs)) {
  design <- head(design, opts$max_pairs)
}
if (nrow(design) == 0) {
  stop("No design rows selected")
}

if (!opts$dry_run) {
  check_package_functions()
}

failures <- list()
for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  err <- tryCatch(
    {
      run_pair(root, row, dry_run = opts$dry_run)
      NULL
    },
    error = function(e) e
  )

  if (!is.null(err)) {
    output_dir <- file.path(root, row$output_dir)
    write_status(output_dir, row, "failed", conditionMessage(err))
    failures[[length(failures) + 1]] <- data.frame(
      pair_id = row$pair_id,
      error = conditionMessage(err),
      stringsAsFactors = FALSE
    )
    message("ERROR in ", row$pair_id, ": ", conditionMessage(err))
    if (!opts$continue_on_error) {
      stop(err)
    }
  }
}

if (length(failures) > 0) {
  failure_df <- do.call(rbind, failures)
  dir.create(file.path(root, "outputs/batch"), recursive = TRUE, showWarnings = FALSE)
  write.table(failure_df, file.path(root, "outputs/batch/brainmapr_failures.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  stop("One or more pairs failed. See outputs/batch/brainmapr_failures.tsv")
}

cat("Selected brainMapR jobs completed successfully.\n")
