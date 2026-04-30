#!/usr/bin/env Rscript

# Create derived brainMapR-ready UKB BWAS inputs by renaming Voxel to Probe.
#
# Purpose:
#   Preserve raw SSTAT_BingJing/*.linear files while creating derived copies
#   with a Probe column for brainMapR.
#
# Inputs:
#   manifests/risk_factor_bwas_manifest.tsv
#
# Outputs:
#   outputs/batch/derived_inputs/*.Probe.linear
#
# How to run:
#   Rscript scripts/02_fix_sumstats_headers.R
#   Rscript scripts/02_fix_sumstats_headers.R --design manifests/brainmapr_small_batch_design.tsv

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  opts <- list(
    project_root = ".",
    manifest = "manifests/risk_factor_bwas_manifest.tsv",
    design = "",
    force = FALSE
  )

  for (arg in args) {
    if (arg == "--force") {
      opts$force <- TRUE
    }
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--project-root") {
      i <- i + 1
      opts$project_root <- args[[i]]
    } else if (key == "--manifest") {
      i <- i + 1
      opts$manifest <- args[[i]]
    } else if (key == "--design") {
      i <- i + 1
      opts$design <- args[[i]]
    } else if (key == "--force") {
      # already handled
    } else if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/02_fix_sumstats_headers.R [--design <pairwise-design.tsv>] [--force]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 1
  }

  opts
}

stream_copy_with_probe_header <- function(input_file, output_file, force = FALSE) {
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  if (file.exists(output_file) && !force) {
    return("skipped_existing")
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  in_con <- file(input_file, open = "r")
  on.exit(close(in_con), add = TRUE)
  out_con <- file(output_file, open = "w")
  on.exit(close(out_con), add = TRUE)

  header <- readLines(in_con, n = 1, warn = FALSE)
  if (length(header) != 1 || !nzchar(header)) {
    stop("Could not read header from: ", input_file)
  }

  columns <- scan(text = header, what = character(), quiet = TRUE)
  if (!"Voxel" %in% columns) {
    stop("Input file does not contain Voxel column: ", input_file)
  }
  if ("Probe" %in% columns) {
    stop("Input already contains Probe column: ", input_file)
  }

  first_data_line <- readLines(in_con, n = 1, warn = FALSE)
  if (length(first_data_line) == 0) {
    stop("Input file has a header but no data rows: ", input_file)
  }

  count_fields <- function(line) {
    length(scan(text = line, what = character(), quiet = TRUE))
  }
  data_field_count <- count_fields(first_data_line)
  expected_field_count <- length(columns)
  has_extra_row_name <- data_field_count == expected_field_count + 1
  if (!has_extra_row_name && data_field_count != expected_field_count) {
    stop(
      "Unexpected data-column count in ", input_file,
      ": header has ", expected_field_count,
      " fields but first data row has ", data_field_count
    )
  }

  drop_first_field <- function(lines) {
    sub("^\\s*(\"[^\"]*\"|\\S+)\\s+", "", lines)
  }

  columns[columns == "Voxel"] <- "Probe"
  # Keep the derived header whitespace-delimited like the original OSCA files.
  # Mixing a tab-delimited header with space-delimited data can make vroom fail
  # to guess the delimiter inside brainMapR/GFA.
  writeLines(paste(columns, collapse = " "), out_con)
  writeLines(
    if (has_extra_row_name) drop_first_field(first_data_line) else first_data_line,
    out_con
  )

  repeat {
    chunk <- readLines(in_con, n = 100000, warn = FALSE)
    if (length(chunk) == 0) {
      break
    }
    if (has_extra_row_name) {
      chunk <- drop_first_field(chunk)
    }
    writeLines(chunk, out_con)
  }

  "created"
}

as_true <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(opts$project_root, mustWork = TRUE)
manifest_path <- file.path(root, opts$manifest)

if (!file.exists(manifest_path)) {
  stop("Manifest does not exist. Run scripts/01_make_manifests.R first: ", manifest_path)
}

manifest <- read.delim(manifest_path, check.names = FALSE)
needed <- manifest[as_true(manifest$needs_Voxel_to_Probe), , drop = FALSE]

if (nzchar(opts$design)) {
  design_path <- file.path(root, opts$design)
  if (!file.exists(design_path)) {
    stop("Design table does not exist: ", design_path)
  }
  design <- read.delim(design_path, check.names = FALSE)
  needed <- needed[needed$analysis_file %in% unique(design$trait_file), , drop = FALSE]
}

if (nrow(needed) == 0) {
  cat("No files require Voxel to Probe conversion.\n")
  quit(save = "no", status = 0)
}

results <- lapply(seq_len(nrow(needed)), function(i) {
  row <- needed[i, ]
  input_file <- file.path(root, row$raw_file)
  output_file <- file.path(root, row$analysis_file)
  status <- stream_copy_with_probe_header(input_file, output_file, force = opts$force)
  data.frame(
    trait_id = row$trait_id,
    input_file = row$raw_file,
    output_file = row$analysis_file,
    status = status,
    stringsAsFactors = FALSE
  )
})

results <- do.call(rbind, results)
log_dir <- file.path(root, "outputs/batch")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
write.table(results, file.path(log_dir, "derived_input_creation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Derived input creation complete.\n")
cat("Created:", sum(results$status == "created"), "\n")
cat("Skipped existing:", sum(results$status == "skipped_existing"), "\n")
cat("Log:", "outputs/batch/derived_input_creation.tsv", "\n")
