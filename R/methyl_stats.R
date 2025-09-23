
## File: R/methyl_stats.R
#' Generate methylation & coverage statistics plots and text summary (wrapper)
#'
#' This function is a wrapper that delegates to small internal functions that
#' handle saving the methylation plot, coverage plot, and writing a text summary.
#'
#' @param obj methylKit object (list-like), usually a methylRawList or methRead output
#' @param plot_dir directory to write plots and summary
#' @param prefix file prefix for outputs
#' @param feature_totals named list with total counts for features (e.g. list(CpGs = 28000000))
#' @param methplot logical whether to save methylation stat plot
#' @param covplot logical whether to save coverage stat plot
#' @import methylKit
#' @export
plot_methylStats <- function(obj, plot_dir = "./results", prefix = "sample", feature_totals = list(CpGs = 28000000), methplot = TRUE, covplot = TRUE) {
  .ensure_dir(plot_dir)
  .ensure_dir(file.path(plot_dir, "MethylationStats"))
  .ensure_dir(file.path(plot_dir, "CoverageStats"))

  if (methplot) .save_methyl_plot(obj, file.path(plot_dir, "MethylationStats", paste0(prefix, "_MethylationStats.png")))
  if (covplot) .save_coverage_plot(obj, file.path(plot_dir, "CoverageStats", paste0(prefix, "_CoverageStats.png")))

  .write_stats_text(obj, plot_dir, prefix, feature_totals)
}

#' @noRd
.save_methyl_plot <- function(obj, filename) {
  .save_png(expr = methylKit::getMethylationStats(obj[[1]], plot = TRUE), filename = filename)
}
#' @noRd
.save_coverage_plot <- function(obj, filename) {
  .save_png(expr = methylKit::getCoverageStats(obj[[1]], plot = TRUE), filename = filename)
}
#' @noRd
.write_stats_text <- function(obj, plot_dir, prefix, feature_totals) {
  stats_file <- file.path(plot_dir, paste0(prefix, "_Stats.txt"))
  if (file.exists(stats_file)) unlink(stats_file)
  sink(stats_file)
  cat("Methylation stats:
")
  print(methylKit::getMethylationStats(obj[[1]], plot = FALSE))
  cat("

Coverage stats:
")
  print(methylKit::getCoverageStats(obj[[1]], plot = FALSE))

  n_covered <- length(as(obj[[1]], "GRanges"))
  cat("

Coverage relative to supplied feature totals:
")
  for (feat in names(feature_totals)) {
    total <- feature_totals[[feat]]
    frac <- round(100 * n_covered / total, 2)
    cat(feat, " covered: ", format(n_covered, big.mark = ","), " / ", format(total, big.mark = ","), " (", frac, "% )
", sep = "")
  }
  cat("
Industry standards for WGBS (human, ~28-29M CpGs):
")
  cat("- Mean/median coverage: ≥10x (robust: 20-30x)
")
  cat("- CpGs covered at ≥1x: ≥70-80%
")
  cat("- CpGs covered at ≥5x: ≥60-70%
")
  cat("- CpGs covered at ≥10x: ≥50-60%
")
  cat("- Bisulfite conversion rate: >99%
")
  cat("- Duplication rate: <20%
")
  sink()
}