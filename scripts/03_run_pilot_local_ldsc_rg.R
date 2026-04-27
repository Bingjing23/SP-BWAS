#!/usr/bin/env Rscript

# Technical smoke-test implementation for SP-BWAS.
#
# This script runs the UKB Alzheimer 40k x hyperTension pilot using brainMapR
# input/reference handling and a local non-official replacement for
# GFA::ldsc_rg(), which is unavailable in CRAN GFA_1.0.5.
#
# This is for technical pipeline validation only. Do not interpret the output
# as final scientific evidence until validated against the original
# GFA::ldsc_rg() implementation expected by brainMapR.

options(stringsAsFactors = FALSE)

shim_warning <- paste(
  "This output was generated with a local non-official ldsc_rg shim because",
  "GFA::ldsc_rg() is unavailable in CRAN GFA_1.0.5. This is for technical",
  "pipeline validation only and must not be used as final scientific evidence."
)

message("WARNING: ", shim_warning)

suppressPackageStartupMessages({
  library(brainMapR)
  library(vroom)
  library(plyr)
  library(bigsnpr)
})

params <- list(
  inputPath = c("outputs/pilot/derived_inputs/", "outputs/pilot/derived_inputs/"),
  bwasFile = c(
    "BWAS_Alzheimer_40k_vtxName.Probe.linear",
    "sstat_FS_All_moda_total_hyperTension.Probe.linear"
  ),
  bwasSampleSize = c("NMISS", "NMISS"),
  refPanel = c("AVERAGE"),
  nblock = 200,
  chi2Threshold = 80,
  varConstrained = TRUE,
  outputPath = "outputs/pilot/AD_UKB40k__hyperTension_local_ldsc_rg"
)

required_columns <- c("Probe", "b", "se", "p", "NMISS")

as_numeric_vector <- function(x) {
  as.numeric(unlist(x, use.names = FALSE))
}

read_sumstats <- function(path) {
  if (!file.exists(path)) {
    stop("Input file is missing: ", path)
  }
  if (file.info(path)$size <= 0) {
    stop("Input file is empty: ", path)
  }
  dat <- vroom::vroom(path, show_col_types = FALSE, progress = FALSE)
  missing_columns <- setdiff(required_columns, names(dat))
  if (length(missing_columns) > 0) {
    stop(
      "Input file is missing required columns: ", path, "\n",
      "Missing: ", paste(missing_columns, collapse = ", "), "\n",
      "Observed columns: ", paste(names(dat), collapse = ", ")
    )
  }
  dat
}

get_brainmapr_function <- function(name) {
  if (exists(name, envir = asNamespace("brainMapR"), inherits = FALSE)) {
    get(name, envir = asNamespace("brainMapR"), inherits = FALSE)
  } else if (exists(name, mode = "function")) {
    get(name, mode = "function")
  } else {
    stop("Could not find required brainMapR helper function: ", name)
  }
}

formatBWAScortical <- get_brainmapr_function("formatBWAScortical")
formatBWASsubcortical <- get_brainmapr_function("formatBWASsubcortical")

resolve_sample_size <- function(dat, sample_size_spec, suffix) {
  numeric_value <- suppressWarnings(as.numeric(sample_size_spec))
  if (!is.na(numeric_value)) {
    return(numeric_value)
  }

  candidates <- unique(c(sample_size_spec, paste0(sample_size_spec, suffix)))
  matched <- candidates[candidates %in% names(dat)]
  if (length(matched) == 0) {
    stop(
      "Could not resolve sample size column. Requested: ", sample_size_spec,
      "; suffix: ", suffix, "; available columns include: ",
      paste(head(names(dat), 40), collapse = ", ")
    )
  }

  value <- as_numeric_vector(dat[[matched[1]]])[1]
  if (!is.finite(value) || value <= 0) {
    stop("Invalid sample size from column ", matched[1], ": ", value)
  }
  value
}

extract_snp_ldsc <- function(ldsc_result) {
  values <- as_numeric_vector(ldsc_result)
  names(values) <- names(unlist(ldsc_result, use.names = TRUE))

  if (length(values) < 4) {
    stop("bigsnpr::snp_ldsc() returned fewer than four numeric values.")
  }

  name_lower <- tolower(names(values))
  int_idx <- which(name_lower %in% c("int", "intercept"))[1]
  int_se_idx <- which(name_lower %in% c("int_se", "intercept_se", "intercept.se"))[1]
  h2_idx <- which(name_lower %in% c("h2", "h2_obs", "m2", "coef", "slope"))[1]
  h2_se_idx <- which(name_lower %in% c("h2_se", "h2_obs_se", "m2_se", "coef_se", "slope_se"))[1]

  if (any(is.na(c(int_idx, int_se_idx, h2_idx, h2_se_idx)))) {
    # bigsnpr::snp_ldsc() commonly returns intercept, intercept_se, h2, h2_se
    # as its first four numeric values. Use that stable order if names differ.
    int_idx <- 1
    int_se_idx <- 2
    h2_idx <- 3
    h2_se_idx <- 4
  }

  c(
    int = unname(values[int_idx]),
    int_se = unname(values[int_se_idx]),
    m2 = unname(values[h2_idx]),
    m2_se = unname(values[h2_se_idx])
  )
}

weighted_covariance_estimate <- function(ld_score, ld_size, z1, z2,
                                         sample_size_1, sample_size_2,
                                         chi2_thr2 = 80) {
  chi2_1 <- z1^2
  chi2_2 <- z2^2
  keep <- is.finite(ld_score) &
    is.finite(z1) & is.finite(z2) &
    is.finite(chi2_1) & is.finite(chi2_2) &
    chi2_1 < chi2_thr2 & chi2_2 < chi2_thr2

  if (sum(keep) < 100) {
    stop("Too few finite rows for cross-trait covariance regression.")
  }

  ld_score <- ld_score[keep]
  z1 <- z1[keep]
  z2 <- z2[keep]

  x <- ld_score / ld_size * sqrt(sample_size_1 * sample_size_2)
  y <- z1 * z2
  weights <- 1 / pmax(ld_score, 1)

  fit <- stats::lm.wfit(x = cbind(1, x), y = y, w = weights)
  coefficients <- as.numeric(fit$coefficients)
  c(int = coefficients[1], braincov = coefficients[2])
}

block_ids <- function(n, blocks) {
  block_count <- min(blocks, n)
  ids <- ceiling(seq_len(n) / n * block_count)
  ids[ids < 1] <- 1
  ids[ids > block_count] <- block_count
  ids
}

jackknife_from_leave_one <- function(full_estimate, leave_one_estimates) {
  h_blocks <- nrow(leave_one_estimates)
  pseudovalues <- h_blocks *
    matrix(rep(full_estimate, each = h_blocks), nrow = h_blocks) -
    (h_blocks - 1) * leave_one_estimates
  colnames(pseudovalues) <- names(full_estimate)
  estimate_j <- colMeans(pseudovalues, na.rm = TRUE)
  se_j <- sqrt(colMeans((pseudovalues - rep(estimate_j, each = h_blocks))^2 / (h_blocks - 1)))
  list(estimate = estimate_j, se = se_j)
}

local_ldsc_rg <- function(ld_score, ld_size, z1, z2, sample_size_1, sample_size_2,
                          blocks = 200, chi2_thr2 = 80) {
  ld_score <- as_numeric_vector(ld_score)
  z1 <- as_numeric_vector(z1)
  z2 <- as_numeric_vector(z2)
  sample_size_1 <- as.numeric(sample_size_1)[1]
  sample_size_2 <- as.numeric(sample_size_2)[1]

  keep <- is.finite(ld_score) & is.finite(z1) & is.finite(z2)
  ld_score <- ld_score[keep]
  z1 <- z1[keep]
  z2 <- z2[keep]

  if (length(ld_score) < 100) {
    stop("Too few finite rows for local_ldsc_rg().")
  }

  trait1_ldsc <- extract_snp_ldsc(
    bigsnpr::snp_ldsc(
      ld_score = ld_score,
      ld_size = ld_size,
      chi2 = z1^2,
      sample_size = sample_size_1,
      blocks = blocks,
      chi2_thr1 = chi2_thr2,
      chi2_thr2 = chi2_thr2
    )
  )

  trait2_ldsc <- extract_snp_ldsc(
    bigsnpr::snp_ldsc(
      ld_score = ld_score,
      ld_size = ld_size,
      chi2 = z2^2,
      sample_size = sample_size_2,
      blocks = blocks,
      chi2_thr1 = chi2_thr2,
      chi2_thr2 = chi2_thr2
    )
  )

  full_cov <- weighted_covariance_estimate(
    ld_score = ld_score,
    ld_size = ld_size,
    z1 = z1,
    z2 = z2,
    sample_size_1 = sample_size_1,
    sample_size_2 = sample_size_2,
    chi2_thr2 = chi2_thr2
  )

  m2_1 <- unname(trait1_ldsc["m2"])
  m2_2 <- unname(trait2_ldsc["m2"])
  braincov <- unname(full_cov["braincov"])

  if (!is.finite(m2_1) || !is.finite(m2_2) || m2_1 <= 0 || m2_2 <= 0) {
    warning("m2_1 <= 0 or m2_2 <= 0; rGM and rGM_se set to NA.")
    rGM <- NA_real_
  } else {
    rGM <- braincov / sqrt(m2_1 * m2_2)
  }

  ids <- block_ids(length(ld_score), blocks)
  h_blocks <- length(unique(ids))
  if (h_blocks < 2) {
    stop("At least two jackknife blocks are required.")
  }

  full_jackknife_estimate <- c(int = unname(full_cov["int"]), braincov = braincov, rGM = rGM)
  leave_one <- matrix(
    NA_real_,
    nrow = h_blocks,
    ncol = length(full_jackknife_estimate),
    dimnames = list(NULL, names(full_jackknife_estimate))
  )

  for (block in seq_len(h_blocks)) {
    use <- ids != block
    cov_leave <- weighted_covariance_estimate(
      ld_score = ld_score[use],
      ld_size = sum(use),
      z1 = z1[use],
      z2 = z2[use],
      sample_size_1 = sample_size_1,
      sample_size_2 = sample_size_2,
      chi2_thr2 = chi2_thr2
    )

    braincov_leave <- unname(cov_leave["braincov"])
    rGM_leave <- if (is.finite(m2_1) && is.finite(m2_2) && m2_1 > 0 && m2_2 > 0) {
      braincov_leave / sqrt(m2_1 * m2_2)
    } else {
      NA_real_
    }
    leave_one[block, ] <- c(
      int = unname(cov_leave["int"]),
      braincov = braincov_leave,
      rGM = rGM_leave
    )
  }

  jack <- jackknife_from_leave_one(full_jackknife_estimate, leave_one)

  out <- c(
    int = jack$estimate["int"],
    int_se = jack$se["int"],
    braincov = jack$estimate["braincov"],
    braincov_se = jack$se["braincov"],
    t1_int = unname(trait1_ldsc["int"]),
    t2_int = unname(trait2_ldsc["int"]),
    m2_1 = m2_1,
    m2_1_se = unname(trait1_ldsc["m2_se"]),
    m2_2 = m2_2,
    m2_2_se = unname(trait2_ldsc["m2_se"]),
    rGM = jack$estimate["rGM"],
    rGM_se = jack$se["rGM"]
  )

  names(out) <- c(
    "int", "int_se",
    "braincov", "braincov_se",
    "t1_int", "t2_int",
    "m2_1", "m2_1_se",
    "m2_2", "m2_2_se",
    "rGM", "rGM_se"
  )

  if (!is.numeric(out) || length(out) != 12 || any(is.na(names(out)))) {
    stop("local_ldsc_rg() did not return the required named numeric vector of length 12.")
  }
  out
}

result_row <- function(panel, ldres) {
  pvalue_1 <- 2 * stats::pnorm(-abs(ldres["m2_1"] / ldres["m2_1_se"]))
  pvalue_2 <- 2 * stats::pnorm(-abs(ldres["m2_2"] / ldres["m2_2_se"]))
  pvalue_rGM <- 2 * stats::pnorm(-abs(ldres["rGM"] / ldres["rGM_se"]))

  data.frame(
    sumR2RefPanel = panel,
    int = ldres["int"],
    int_se = ldres["int_se"],
    braincov = ldres["braincov"],
    braincov_se = ldres["braincov_se"],
    t1_int = ldres["t1_int"],
    t2_int = ldres["t2_int"],
    m2_1 = ldres["m2_1"],
    m2_1_se = ldres["m2_1_se"],
    pvalue_1 = pvalue_1,
    m2_1_CI_lb = ldres["m2_1"] - 1.96 * ldres["m2_1_se"],
    m2_1_CI_ub = ldres["m2_1"] + 1.96 * ldres["m2_1_se"],
    m2_2 = ldres["m2_2"],
    m2_2_se = ldres["m2_2_se"],
    pvalue_2 = pvalue_2,
    m2_2_CI_lb = ldres["m2_2"] - 1.96 * ldres["m2_2_se"],
    m2_2_CI_ub = ldres["m2_2"] + 1.96 * ldres["m2_2_se"],
    rGM = ldres["rGM"],
    rGM_se = ldres["rGM_se"],
    pvalue_rGM = pvalue_rGM,
    rGM_CI_lb = ldres["rGM"] - 1.96 * ldres["rGM_se"],
    rGM_CI_ub = ldres["rGM"] + 1.96 * ldres["rGM_se"],
    check.names = FALSE
  )
}

dir.create(params$outputPath, recursive = TRUE, showWarnings = FALSE)
writeLines(
  shim_warning,
  file.path(params$outputPath, "README_NONOFFICIAL_SHIM.txt")
)

message("Reading input BWAS files.")
bwas_path_1 <- paste0(params$inputPath[1], params$bwasFile[1])
bwas_path_2 <- paste0(params$inputPath[2], params$bwasFile[2])
bwas_1 <- read_sumstats(bwas_path_1)
bwas_2 <- read_sumstats(bwas_path_2)

message("Formatting first BWAS with brainMapR helper functions.")
bwas_annotated <- NULL
for (mod in c("thickness", "area", "thick", "LogJacs")) {
  for (hemi in c("lh", "rh")) {
    formatted <- if (mod %in% c("area", "thickness")) {
      formatBWAScortical(BWASsumstat = bwas_1, hemi = hemi, mod = mod)
    } else {
      formatBWASsubcortical(BWASsumstat = bwas_1, hemi = hemi, mod = mod)
    }
    formatted$CHI2 <- (formatted$b / formatted$se)^2
    bwas_annotated <- plyr::rbind.fill(bwas_annotated, formatted)
  }
}

for (panel in params$refPanel) {
  sumr2_col <- paste0("SumR2_", panel)
  if (!sumr2_col %in% names(bwas_annotated)) {
    stop("Required reference-panel column missing after formatting: ", sumr2_col)
  }
}

message("Merging BWAS files by Probe.")
bwas_merged <- merge(
  bwas_annotated,
  bwas_2,
  by = "Probe",
  suffixes = c("_1", "_2")
)

required_merged <- c("b_1", "se_1", "p_1", "b_2", "se_2", "p_2")
missing_merged <- setdiff(required_merged, names(bwas_merged))
if (length(missing_merged) > 0) {
  stop("Merged data is missing required columns: ", paste(missing_merged, collapse = ", "))
}

bwas_merged$z_1 <- as_numeric_vector(bwas_merged$b_1) / as_numeric_vector(bwas_merged$se_1)
bwas_merged$z_2 <- as_numeric_vector(bwas_merged$b_2) / as_numeric_vector(bwas_merged$se_2)

missing_p <- is.na(bwas_merged$p_1) | is.na(bwas_merged$p_2)
if (any(missing_p)) {
  bwas_merged <- bwas_merged[!missing_p, ]
}

sample_size_1 <- resolve_sample_size(
  bwas_merged,
  sample_size_spec = params$bwasSampleSize[1],
  suffix = "_1"
)
sample_size_2 <- resolve_sample_size(
  bwas_merged,
  sample_size_spec = params$bwasSampleSize[2],
  suffix = "_2"
)

message("Resolved sample size 1: ", sample_size_1)
message("Resolved sample size 2: ", sample_size_2)

all_results <- NULL
for (panel in params$refPanel) {
  sumr2_col <- paste0("SumR2_", panel)
  panel_data <- bwas_merged[!is.na(bwas_merged[[sumr2_col]]), ]
  if (nrow(panel_data) < 100) {
    stop("Too few rows with non-missing ", sumr2_col, ".")
  }

  message(
    "Estimating rGM using local_ldsc_rg. Using ",
    nrow(panel_data), " vertices and SumR2 from reference panel: ", panel
  )

  ldres <- local_ldsc_rg(
    ld_score = panel_data[[sumr2_col]],
    ld_size = nrow(panel_data),
    z1 = panel_data$z_1,
    z2 = panel_data$z_2,
    sample_size_1 = sample_size_1,
    sample_size_2 = sample_size_2,
    blocks = params$nblock,
    chi2_thr2 = params$chi2Threshold
  )

  if (isTRUE(params$varConstrained)) {
    if (is.finite(ldres["m2_1"]) && ldres["m2_1"] > 1) {
      ldres["m2_1"] <- 1
      ldres["rGM"] <- ldres["braincov"] / sqrt(ldres["m2_1"] * ldres["m2_2"])
      message("Variance of trait 1 constrained to 1.")
    }
    if (is.finite(ldres["m2_2"]) && ldres["m2_2"] > 1) {
      ldres["m2_2"] <- 1
      ldres["rGM"] <- ldres["braincov"] / sqrt(ldres["m2_1"] * ldres["m2_2"])
      message("Variance of trait 2 constrained to 1.")
    }
    if (is.finite(ldres["m2_1"]) && ldres["m2_1"] < 0) {
      ldres["m2_1"] <- 0.00001
      ldres["rGM"] <- ldres["braincov"] / sqrt(ldres["m2_1"] * ldres["m2_2"])
      message("Variance of trait 1 constrained to 0.")
    }
    if (is.finite(ldres["m2_2"]) && ldres["m2_2"] < 0) {
      ldres["m2_2"] <- 0.00001
      ldres["rGM"] <- ldres["braincov"] / sqrt(ldres["m2_1"] * ldres["m2_2"])
      message("Variance of trait 2 constrained to 0.")
    }
  }

  all_results <- rbind(all_results, result_row(panel, ldres))
}

output_file <- file.path(
  params$outputPath,
  "SumR2_regression_bivariate_local_ldsc_rg.rsq"
)
write.table(all_results, output_file, col.names = TRUE, row.names = FALSE, quote = FALSE)

message("Wrote output: ", output_file)
message("WARNING: ", shim_warning)
print(all_results)
sessionInfo()
