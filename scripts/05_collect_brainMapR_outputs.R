#!/usr/bin/env Rscript

# Collect manifest-driven brainMapR pair outputs into summary tables.
#
# Purpose:
#   Index all per-pair output files and extract scalar tabular results where
#   possible into a long table and metric matrices.
#
# Inputs:
#   manifests/brainmapr_pairwise_design.tsv or a clean/sensitivity design
#
# Outputs:
#   outputs/batch/summary/brainmapr_output_files.tsv
#   outputs/batch/summary/brainmapr_pairwise_results_long.tsv
#   outputs/batch/summary/matrix_*.tsv
#
# How to run:
#   Rscript scripts/05_collect_brainMapR_outputs.R \
#     --design manifests/brainmapr_clean_average_design.tsv \
#     --output-dir outputs/batch/summary_clean_AVERAGE

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  opts <- list(
    project_root = ".",
    design = "manifests/brainmapr_pairwise_design.tsv",
    output_dir = "outputs/batch/summary"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--project-root") {
      i <- i + 1
      opts$project_root <- args[[i]]
    } else if (key == "--design") {
      i <- i + 1
      opts$design <- args[[i]]
    } else if (key == "--output-dir") {
      i <- i + 1
      opts$output_dir <- args[[i]]
    } else if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/05_collect_brainMapR_outputs.R [--design <tsv>]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 1
  }

  opts
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

try_read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("rsq", "tsv", "csv", "txt")) {
    return(NULL)
  }
  if (file.info(path)$size == 0) {
    return(NULL)
  }

  sep <- if (ext == "csv") "," else ""
  tryCatch(
    read.table(path, header = TRUE, sep = sep, check.names = FALSE,
               stringsAsFactors = FALSE, comment.char = "", quote = "\""),
    error = function(e) NULL
  )
}

extract_scalar_metrics <- function(path) {
  dat <- try_read_table(path)
  if (is.null(dat) || nrow(dat) == 0) {
    return(data.frame())
  }

  first_row <- dat[1, , drop = FALSE]
  rows <- lapply(names(first_row), function(col) {
    value <- first_row[[col]][1]
    if (length(value) == 0 || is.na(value)) {
      return(NULL)
    }
    data.frame(
      metric = col,
      value = as.character(value),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

write_metric_matrices <- function(results, out_dir) {
  if (nrow(results) == 0) {
    return(invisible(NULL))
  }

  metrics <- unique(results$metric)
  for (metric in metrics) {
    dat <- results[results$metric == metric, , drop = FALSE]
    numeric_value <- suppressWarnings(as.numeric(dat$value))
    if (all(is.na(numeric_value))) {
      next
    }
    dat$value_numeric <- numeric_value
    dat <- dat[!is.na(dat$value_numeric), , drop = FALSE]
    if (nrow(dat) == 0) {
      next
    }

    ad_ids <- sort(unique(dat$ad_id))
    trait_ids <- sort(unique(dat$trait_id))
    mat <- matrix(NA_real_, nrow = length(ad_ids), ncol = length(trait_ids),
                  dimnames = list(ad_ids, trait_ids))
    for (i in seq_len(nrow(dat))) {
      mat[dat$ad_id[[i]], dat$trait_id[[i]]] <- dat$value_numeric[[i]]
    }

    matrix_df <- data.frame(ad_id = rownames(mat), mat, check.names = FALSE)
    out_file <- file.path(out_dir, paste0("matrix_", sanitize_name(metric), ".tsv"))
    write.table(matrix_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(opts$project_root, mustWork = TRUE)
design_path <- file.path(root, opts$design)
out_dir <- file.path(root, opts$output_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(design_path)) {
  stop("Design table does not exist: ", design_path)
}

design <- read.delim(design_path, check.names = FALSE)
file_rows <- list()
metric_rows <- list()

for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  output_dir <- file.path(root, row$output_dir)
  files <- if (dir.exists(output_dir)) {
    list.files(output_dir, recursive = TRUE, full.names = TRUE)
  } else {
    character()
  }

  if (length(files) == 0) {
    file_rows[[length(file_rows) + 1]] <- data.frame(
      pair_id = row$pair_id,
      ad_id = row$ad_id,
      trait_id = row$trait_id,
      output_dir = row$output_dir,
      file_path = "",
      size_bytes = NA_real_,
      modified_time = "",
      stringsAsFactors = FALSE
    )
    next
  }

  for (path in files) {
    info <- file.info(path)
    rel_path <- sub(paste0("^", root, "/?"), "", gsub("\\\\", "/", path))
    file_rows[[length(file_rows) + 1]] <- data.frame(
      pair_id = row$pair_id,
      ad_id = row$ad_id,
      trait_id = row$trait_id,
      output_dir = row$output_dir,
      file_path = rel_path,
      size_bytes = info$size,
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )

    metrics <- extract_scalar_metrics(path)
    if (nrow(metrics) > 0) {
      metrics$pair_id <- row$pair_id
      metrics$ad_id <- row$ad_id
      metrics$trait_id <- row$trait_id
      metrics$source_file <- rel_path
      metric_rows[[length(metric_rows) + 1]] <- metrics[
        c("pair_id", "ad_id", "trait_id", "source_file", "metric", "value")
      ]
    }
  }
}

file_index <- do.call(rbind, file_rows)
metrics_long <- if (length(metric_rows) > 0) {
  do.call(rbind, metric_rows)
} else {
  data.frame(
    pair_id = character(),
    ad_id = character(),
    trait_id = character(),
    source_file = character(),
    metric = character(),
    value = character(),
    stringsAsFactors = FALSE
  )
}

write.table(file_index, file.path(out_dir, "brainmapr_output_files.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metrics_long, file.path(out_dir, "brainmapr_pairwise_results_long.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write_metric_matrices(metrics_long, out_dir)

cat("Collection complete.\n")
cat("Indexed files:", nrow(file_index), "\n")
cat("Scalar metric rows:", nrow(metrics_long), "\n")
cat("Output directory:", opts$output_dir, "\n")
