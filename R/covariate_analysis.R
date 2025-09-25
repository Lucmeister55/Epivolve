## File: R/covariate_analysis.R

#' Compute per-sample mean methylation (or metric)
#' @param mat numeric matrix of values (rows = features, cols = samples)
#' @return numeric vector of mean per sample
#' @export
compute_sample_means <- function(mat) {
  colMeans(mat, na.rm = TRUE)
}

#' Prepare metadata dataframe with sample means
#' @param sample_means numeric vector of per-sample means
#' @param metadata data.frame containing covariates (rows = samples)
#' @return data.frame with mean metric and covariates
#' @export
prepare_metadata_df <- function(sample_means, metadata) {
  df <- data.frame(sample = names(sample_means), mean_metric = sample_means)
  df <- cbind(df, metadata[match(names(sample_means), metadata$sample), , drop = FALSE])
  df
}

#' Run ANOVA for a single covariate vs mean metric
#' Handles numeric covariates by binning into quantiles
#' @param df data.frame with columns mean_metric and covariate
#' @param covariate name of covariate (string)
#' @param bins number of quantiles to bin numeric covariates (default 4)
#' @return list with p_value and plotting variable
#' @export
run_covariate_anova <- function(df, covariate, bins = 4) {
  # Only rows with non-NA covariate
  df_cc <- df[!is.na(df[[covariate]]), , drop = FALSE]

  plot_var <- covariate

  # Bin numeric covariates into quantiles
  if (is.numeric(df_cc[[covariate]])) {
    plot_var <- paste0(covariate, "_binned")
    # cut() automatically returns a factor
    df_cc[[plot_var]] <- cut(
      df_cc[[covariate]],
      breaks = quantile(df_cc[[covariate]], probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE),
      include.lowest = TRUE
    )
  } else {
    # Ensure categorical/factor covariates are factors
    df_cc[[plot_var]] <- as.factor(df_cc[[plot_var]])
  }

  # Run ANOVA if at least 2 levels
  if (length(unique(df_cc[[plot_var]])) < 2) {
    p_value <- NA
  } else {
    anova_res <- aov(mean_metric ~ get(plot_var), data = df_cc)
    p_value <- summary(anova_res)[[1]][["Pr(>F)"]][1]
  }

  # Return complete-case binned data for plotting
  list(p_value = p_value, plot_var = plot_var, df = df_cc)
}

#' Generate boxplot of mean metric vs covariate
#' @param df data.frame
#' @param plot_var variable for x-axis (binned if numeric)
#' @param p_value ANOVA p-value to append to title
#' @param title plot title
#' @return ggplot object
#' @export
plot_covariate_boxplot <- function(df, plot_var, p_value, title = NULL) {
  # Append ANOVA p-value to title
  full_title <- if (!is.null(title)) {
    paste0(title, " (ANOVA p = ", signif(p_value, 3), ")")
  } else {
    paste0("ANOVA p = ", signif(p_value, 3))
  }

  # Use your package's box_plot with updated title
  box_plot(
    data = df,
    x = !!rlang::sym(plot_var),
    y = mean_metric,
    title = full_title
  )
}

#' Top-level wrapper: compute sample means, prepare metadata, run ANOVAs and save plots
#'
#' @param mat numeric matrix (rows = features, cols = samples)
#' @param metadata data.frame of sample covariates (must have column `sample`)
#' @param covariates character vector of covariate names to test
#' @param output_dir directory to save plots
#' @param p_threshold significance threshold for storing results
#' @return data.frame of significant covariates and p-values
#' @export
analyze_covariates <- function(mat, metadata, covariates, output_dir, p_threshold = 0.05) {
  # Step 1: compute per-sample mean metric
  sample_means <- compute_sample_means(mat)

  # Step 2: prepare metadata dataframe
  df <- prepare_metadata_df(sample_means, metadata)

  # Step 3: create output folder if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Step 4: loop over covariates, run ANOVA, save boxplots
  sig_results <- data.frame(covariate = character(), p_value = numeric(), stringsAsFactors = FALSE)

  for (cov in covariates) {
    res <- run_covariate_anova(df, cov)
    if (is.null(res)) next  # skip if NA values

    # Boxplot with ANOVA p-value in title
    p <- plot_covariate_boxplot(res$df, res$plot_var, res$p_value, title = paste("Mean metric by", cov))
    ggsave(file.path(output_dir, paste0("mean_metric_by_", cov, ".png")),
           plot = p, width = 6, height = 4, dpi = 300)

    # Store significant
    if (!is.na(res$p_value) && res$p_value < p_threshold) {
      sig_results <- rbind(sig_results, data.frame(covariate = cov, p_value = res$p_value))
    }
  }

  if (nrow(sig_results) == 0) {
    message("No significant covariates at p <", p_threshold)
  }

  sig_results
}
