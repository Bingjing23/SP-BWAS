#!/usr/bin/env Rscript

# Create publication-review figures from collected brainMapR outputs.
#
# Inputs:
#   outputs/batch/summary*/matrix_rGM.tsv
#   outputs/batch/summary*/matrix_pvalue_rGM.tsv
#   outputs/batch/summary*/matrix_rGM_CI_lb.tsv
#   outputs/batch/summary*/matrix_rGM_CI_ub.tsv
#   optional pairwise design with readable labels
#
# Outputs:
#   PNG and TIFF figures
#   top_associations_bonferroni.tsv
#   out_of_range_rGM.tsv
#   figure_generation_report.txt
#
# How to run:
#   Rscript scripts/06_plot_brainMapR_summary.R \
#     --summary-dir outputs/batch/summary_clean_AVERAGE \
#     --design manifests/brainmapr_clean_average_design.tsv
#   Rscript scripts/06_plot_brainMapR_summary.R \
#     --summary-dir outputs/batch/summary_clean_AVERAGE \
#     --design manifests/brainmapr_clean_average_design.tsv \
#     --compare-summary-dir outputs/batch/summary_clean_UKB

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  opts <- list(
    summary_dir = "outputs/batch/summary",
    design = "",
    compare_summary_dir = "",
    out_dir = ""
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--summary-dir") {
      i <- i + 1
      opts$summary_dir <- args[[i]]
    } else if (key == "--design") {
      i <- i + 1
      opts$design <- args[[i]]
    } else if (key == "--compare-summary-dir") {
      i <- i + 1
      opts$compare_summary_dir <- args[[i]]
    } else if (key == "--out-dir") {
      i <- i + 1
      opts$out_dir <- args[[i]]
    } else if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/06_plot_brainMapR_summary.R --summary-dir <dir> [--design <tsv>] [--compare-summary-dir <dir>]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 1
  }

  if (!nzchar(opts$out_dir)) {
    opts$out_dir <- file.path(opts$summary_dir, "figures")
  }
  opts
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg)
  }
}

read_matrix_long <- function(path, value_name) {
  if (!file.exists(path)) {
    stop("Missing matrix file: ", path)
  }
  mat <- read.delim(path, check.names = FALSE)
  if (!"ad_id" %in% names(mat)) {
    stop("Matrix is missing ad_id column: ", path)
  }
  value_cols <- setdiff(names(mat), "ad_id")
  rows <- vector("list", length(value_cols))
  for (i in seq_along(value_cols)) {
    trait_id <- value_cols[[i]]
    rows[[i]] <- data.frame(
      ad_id = mat$ad_id,
      trait_id = trait_id,
      value = suppressWarnings(as.numeric(mat[[trait_id]])),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  names(out)[names(out) == "value"] <- value_name
  out
}

ad_label_default <- function(ad_id) {
  labels <- c(
    ADvsHC = "AD dementia vs HC",
    MCIvsHC = "MCI vs HC",
    Conversion1year = "Conversion <=1y",
    Conversion2years = "Conversion <=2y",
    Conversion3years = "Conversion <=3y",
    Conversion4years = "Conversion <=4y",
    Conversion5years = "Conversion <=5y",
    MMSE = "MMSE"
  )
  ifelse(ad_id %in% names(labels), labels[ad_id], ad_id)
}

category_label <- function(category) {
  labels <- c(
    vascular_metabolic = "Vascular/metabolic",
    psychiatric = "Psychiatric",
    lifestyle = "Lifestyle",
    social_socioeconomic = "Social/SES",
    frailty_multimorbidity = "Frailty/multimorbidity",
    infection_inflammatory = "Infection/inflammatory",
    hormonal_sex_specific = "Hormonal/sex-specific",
    other_candidate = "Other",
    other = "Other"
  )
  ifelse(category %in% names(labels), labels[category], category)
}

infer_category <- function(trait_id) {
  id <- tolower(trait_id)
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
  "other"
}

read_label_map <- function(design_path) {
  if (!nzchar(design_path) || !file.exists(design_path)) {
    return(NULL)
  }
  design <- read.delim(design_path, check.names = FALSE)
  trait_cols <- intersect(c("trait_id", "trait_label", "trait_category"), names(design))
  ad_cols <- intersect(c("ad_id", "ad_label"), names(design))
  list(
    traits = unique(design[trait_cols]),
    ads = unique(design[ad_cols])
  )
}

load_summary <- function(summary_dir, design_path = "") {
  r <- read_matrix_long(file.path(summary_dir, "matrix_rGM.tsv"), "rGM")
  p <- read_matrix_long(file.path(summary_dir, "matrix_pvalue_rGM.tsv"), "pvalue")
  lb <- read_matrix_long(file.path(summary_dir, "matrix_rGM_CI_lb.tsv"), "ci_lb")
  ub <- read_matrix_long(file.path(summary_dir, "matrix_rGM_CI_ub.tsv"), "ci_ub")

  records <- merge(r, p, by = c("ad_id", "trait_id"), all = TRUE)
  records <- merge(records, lb, by = c("ad_id", "trait_id"), all = TRUE)
  records <- merge(records, ub, by = c("ad_id", "trait_id"), all = TRUE)

  labels <- read_label_map(design_path)
  records$ad_label <- ad_label_default(records$ad_id)
  records$trait_label <- gsub("_", " ", records$trait_id)
  records$trait_category <- vapply(records$trait_id, infer_category, character(1))

  if (!is.null(labels)) {
    if (all(c("ad_id", "ad_label") %in% names(labels$ads))) {
      records <- merge(records, labels$ads, by = "ad_id", all.x = TRUE, suffixes = c("", ".design"))
      records$ad_label <- ifelse(!is.na(records$ad_label.design), records$ad_label.design, records$ad_label)
      records$ad_label.design <- NULL
    }
    if (all(c("trait_id", "trait_label") %in% names(labels$traits))) {
      records <- merge(records, labels$traits, by = "trait_id", all.x = TRUE, suffixes = c("", ".design"))
      records$trait_label <- ifelse(!is.na(records$trait_label.design), records$trait_label.design, records$trait_label)
      records$trait_label.design <- NULL
      if ("trait_category.design" %in% names(records)) {
        records$trait_category <- ifelse(!is.na(records$trait_category.design), records$trait_category.design, records$trait_category)
        records$trait_category.design <- NULL
      }
    }
  }

  n_tests <- sum(is.finite(records$pvalue))
  records$bonferroni_p <- pmin(records$pvalue * n_tests, 1)
  records$fwer_significant <- is.finite(records$bonferroni_p) & records$bonferroni_p < 0.05
  records$out_of_range <- is.finite(records$rGM) & abs(records$rGM) > 1
  records$rGM_clipped <- pmax(pmin(records$rGM, 1), -1)
  records$neg_log10_p <- pmin(-log10(records$pvalue), 20)
  records$category_label <- category_label(records$trait_category)
  records
}

order_records <- function(records) {
  ad_order <- c(
    "ADvsHC", "MCIvsHC", "Conversion1year", "Conversion2years",
    "Conversion3years", "Conversion4years", "Conversion5years", "MMSE"
  )
  ad_ids <- c(ad_order[ad_order %in% records$ad_id], setdiff(unique(records$ad_id), ad_order))
  ad_labels <- unique(records[match(ad_ids, records$ad_id), c("ad_id", "ad_label")])

  category_order <- c(
    "Vascular/metabolic", "Psychiatric", "Lifestyle", "Social/SES",
    "Frailty/multimorbidity", "Infection/inflammatory",
    "Hormonal/sex-specific", "Other"
  )
  trait_stats <- aggregate(
    abs(records$rGM),
    by = list(
      trait_id = records$trait_id,
      trait_label = records$trait_label,
      category_label = records$category_label
    ),
    FUN = function(x) max(x, na.rm = TRUE)
  )
  names(trait_stats)[names(trait_stats) == "x"] <- "max_abs_rGM"
  trait_stats$category_label <- factor(trait_stats$category_label, levels = category_order)
  trait_stats <- trait_stats[order(trait_stats$category_label, -trait_stats$max_abs_rGM, trait_stats$trait_label), ]

  records$ad_label <- factor(records$ad_label, levels = ad_labels$ad_label)
  records$trait_label <- factor(records$trait_label, levels = rev(trait_stats$trait_label))
  records$category_label <- factor(records$category_label, levels = category_order)
  records
}

save_plot <- function(plot, out_dir, stem, width, height) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  tiff_file <- file.path(out_dir, paste0(stem, ".tiff"))
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(tiff_file, plot, width = width, height = height, units = "in",
                  dpi = 300, device = "tiff", compression = "lzw")
}

write_table <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

make_plots <- function(records, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  records <- order_records(records)

  n_traits <- length(unique(records$trait_label))
  height <- max(8, 2.4 + 0.24 * n_traits)
  width <- 9.5

  heatmap <- ggplot2::ggplot(records, ggplot2::aes(x = ad_label, y = trait_label, fill = rGM_clipped)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_point(
      data = records[records$fwer_significant, ],
      ggplot2::aes(x = ad_label, y = trait_label),
      inherit.aes = FALSE,
      size = 1.5,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = records[records$out_of_range, ],
      ggplot2::aes(x = ad_label, y = trait_label, label = "!"),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    ) +
    ggplot2::facet_grid(category_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#2c7bb6", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1), name = "rGM"
    ) +
    ggplot2::labs(
      title = "AD x risk-factor grey-matter correlations",
      subtitle = "Black dot: Bonferroni/FWER < 0.05. Exclamation mark: |rGM| > 1, flagged as unstable.",
      x = "AD-related BWAS map",
      y = "UKB risk-factor BWAS map"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 8),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  save_plot(heatmap, out_dir, "figure_1_vertical_rGM_heatmap", width, height)

  pheatmap <- ggplot2::ggplot(records, ggplot2::aes(x = ad_label, y = trait_label, fill = neg_log10_p)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_point(
      data = records[records$fwer_significant, ],
      ggplot2::aes(x = ad_label, y = trait_label),
      inherit.aes = FALSE,
      size = 1.5,
      color = "black"
    ) +
    ggplot2::facet_grid(category_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient(low = "#fff5f0", high = "#67000d", limits = c(0, 20), name = "-log10(p)") +
    ggplot2::labs(
      title = "AD x risk-factor rGM statistical evidence",
      subtitle = "Color shows -log10(p), clipped at 20. Black dot: Bonferroni/FWER < 0.05.",
      x = "AD-related BWAS map",
      y = "UKB risk-factor BWAS map"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 8),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  save_plot(pheatmap, out_dir, "figure_2_vertical_pvalue_heatmap", width, height)

  top <- records[
    is.finite(records$rGM) &
      is.finite(records$ci_lb) &
      is.finite(records$ci_ub) &
      records$fwer_significant &
      !records$out_of_range,
  ]
  top <- top[order(top$bonferroni_p, -abs(top$rGM)), ]
  top <- head(top, 30)
  top$pair_label <- factor(
    paste(top$ad_label, top$trait_label, sep = " x "),
    levels = rev(paste(top$ad_label, top$trait_label, sep = " x "))
  )
  forest <- ggplot2::ggplot(top, ggplot2::aes(x = rGM, y = pair_label, color = category_label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lb, xmax = ci_ub), height = 0.12, linewidth = 0.45) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
    ggplot2::labs(
      title = "Top stable rGM associations",
      subtitle = "Bonferroni/FWER-significant pairs with |rGM| <= 1. Error bars show 95% CI.",
      x = "rGM",
      y = NULL,
      color = "Category"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))
  save_plot(forest, out_dir, "figure_3_top_rGM_forest", 10, max(6, 0.25 * nrow(top) + 2))

  top_out <- top[, c(
    "ad_id", "ad_label", "trait_id", "trait_label", "trait_category",
    "rGM", "ci_lb", "ci_ub", "pvalue", "bonferroni_p"
  )]
  write_table(top_out, file.path(dirname(out_dir), "top_associations_bonferroni.tsv"))

  out_of_range <- records[records$out_of_range, c(
    "ad_id", "ad_label", "trait_id", "trait_label", "trait_category",
    "rGM", "ci_lb", "ci_ub", "pvalue", "bonferroni_p"
  )]
  out_of_range <- out_of_range[order(-abs(out_of_range$rGM)), ]
  write_table(out_of_range, file.path(dirname(out_dir), "out_of_range_rGM.tsv"))

  list(top = top_out, out_of_range = out_of_range)
}

make_sensitivity_plot <- function(avg_records, ukb_summary_dir, design_path, out_dir) {
  ukb_records <- load_summary(ukb_summary_dir, design_path)
  merged <- merge(
    avg_records[, c("ad_id", "trait_id", "ad_label", "trait_label", "trait_category", "category_label", "rGM")],
    ukb_records[, c("ad_id", "trait_id", "rGM")],
    by = c("ad_id", "trait_id"),
    suffixes = c("_AVERAGE", "_UKB")
  )
  merged <- merged[is.finite(merged$rGM_AVERAGE) & is.finite(merged$rGM_UKB), ]
  r <- if (nrow(merged) >= 3) cor(merged$rGM_AVERAGE, merged$rGM_UKB) else NA_real_
  label <- paste0("Pearson r = ", ifelse(is.finite(r), sprintf("%.3f", r), "NA"), "\nN = ", nrow(merged))

  plot <- ggplot2::ggplot(merged, ggplot2::aes(x = rGM_AVERAGE, y = rGM_UKB, color = category_label)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(alpha = 0.75, size = 1.8) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    ggplot2::annotate("text", x = -0.95, y = 0.95, hjust = 0, vjust = 1, label = label, size = 3.5) +
    ggplot2::labs(
      title = "Reference-panel sensitivity of rGM estimates",
      subtitle = "Each point is one AD x UKB risk-factor pair; dashed line is y = x.",
      x = "rGM using AVERAGE reference panel",
      y = "rGM using UKB reference panel",
      color = "Category"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  save_plot(plot, out_dir, "figure_4_reference_panel_sensitivity", 7.2, 6.2)
  write_table(merged, file.path(dirname(out_dir), "reference_panel_rGM_comparison.tsv"))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  require_package("ggplot2")
  dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

  records <- load_summary(opts$summary_dir, opts$design)
  plots <- make_plots(records, opts$out_dir)

  if (nzchar(opts$compare_summary_dir)) {
    make_sensitivity_plot(records, opts$compare_summary_dir, opts$design, opts$out_dir)
  }

  report <- c(
    paste("Figures written to:", opts$out_dir),
    paste("AD maps:", length(unique(records$ad_id))),
    paste("Risk-factor maps:", length(unique(records$trait_id))),
    paste("Collected pairs:", nrow(records)),
    paste("Bonferroni/FWER-significant pairs:", sum(records$fwer_significant, na.rm = TRUE)),
    paste("Out-of-range rGM estimates:", sum(records$out_of_range, na.rm = TRUE)),
    paste("Top stable associations written:", nrow(plots$top))
  )
  writeLines(report, file.path(opts$summary_dir, "figure_generation_report.txt"))
  cat(paste(report, collapse = "\n"), "\n")
}

main()
