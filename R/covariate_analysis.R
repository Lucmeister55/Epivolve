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
# Robust run_covariate_anova that tolerates non-unique quantile breaks
run_covariate_anova <- function(df, covariate, bins = 4) {
  # Only rows with non-NA covariate
  df_cc <- df[!is.na(df[[covariate]]), , drop = FALSE]
  plot_var <- covariate

  # If no rows remain, return NA
  if (nrow(df_cc) == 0) {
    return(list(p_value = NA, plot_var = plot_var, df = df_cc))
  }

  # Bin numeric covariates into quantiles safely
  if (is.numeric(df_cc[[covariate]])) {
    plot_var <- paste0(covariate, "_binned")
    probs <- seq(0, 1, length.out = bins + 1)
    breaks <- unique(quantile(df_cc[[covariate]], probs = probs, na.rm = TRUE))

    if (length(breaks) < 2) {
      # not enough unique values to form bins -> create single-level factor and return NA p-value
      warning(sprintf("Covariate '%s' has too few unique values to bin (unique breaks = %d). Skipping ANOVA.", covariate, length(breaks)))
      df_cc[[plot_var]] <- factor(rep("all", nrow(df_cc)))
      return(list(p_value = NA, plot_var = plot_var, df = df_cc))
    }

    # If there are fewer unique breaks than requested bins+1, cut() will still work
    df_cc[[plot_var]] <- cut(df_cc[[covariate]],
                             breaks = breaks,
                             include.lowest = TRUE,
                             dig.lab = 10)
  } else {
    # Ensure categorical/factor covariates are factors
    df_cc[[plot_var]] <- as.factor(df_cc[[plot_var]])
  }

  # Run ANOVA if at least 2 levels
  if (length(unique(df_cc[[plot_var]])) < 2) {
    p_value <- NA
  } else {
    # use tryCatch around aov to be safe
    fm <- as.formula(paste("mean_metric ~", plot_var))
    anova_res <- tryCatch(aov(fm, data = df_cc),
                          error = function(e) {
                            warning(sprintf("ANOVA failed for covariate '%s': %s", covariate, e$message))
                            NULL
                          })
    if (is.null(anova_res)) {
      p_value <- NA
    } else {
      p_value <- tryCatch({
        summary(anova_res)[[1]][["Pr(>F)"]][1]
      }, error = function(e) {
        warning(sprintf("Failed extracting p-value for '%s': %s", covariate, e$message))
        NA
      })
    }
  }

  list(p_value = p_value, plot_var = plot_var, df = df_cc)
}

# Slightly more defensive plotting function
plot_covariate_boxplot <- function(df, plot_var, p_value, title = NULL) {
  full_title <- if (!is.null(title)) {
    paste0(title, " (ANOVA p = ", signif(p_value, 3), ")")
  } else {
    paste0("ANOVA p = ", signif(p_value, 3))
  }

  # Use box_plot (assumed to exist in your package). Wrap in tryCatch so failures don't abort loop.
  p <- tryCatch({
    box_plot(
      data = df,
      x = !!rlang::sym(plot_var),
      y = mean_metric,
      title = full_title
    )
  }, error = function(e) {
    warning(sprintf("Plotting failed for '%s': %s", plot_var, e$message))
    NULL
  })

  p
}

# Top-level wrapper with robust error handling around individual covariates
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
    res <- tryCatch(run_covariate_anova(df, cov), error = function(e) {
      warning(sprintf("Error while processing covariate '%s': %s", cov, e$message))
      NULL
    })
    if (is.null(res)) next

    # Attempt to make plot (may be NULL)
    p <- tryCatch(plot_covariate_boxplot(res$df, res$plot_var, res$p_value, title = paste("Mean metric by", cov)),
                  error = function(e) {
                    warning(sprintf("Plot generation error for '%s': %s", cov, e$message))
                    NULL
                  })

    # Save plot only if it was successfully created
    if (!is.null(p)) {
      out_file <- file.path(output_dir, paste0("mean_metric_by_", cov, ".png"))
      tryCatch({
        ggsave(out_file, plot = p, width = 6, height = 4, dpi = 300)
      }, error = function(e) {
        warning(sprintf("Failed to save plot for '%s' to '%s': %s", cov, out_file, e$message))
      })
    }

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
#' PCA with covariate association heatmaps (p-value & effect size)
#'
#' Computes PCA on a data matrix, quantifies association of metadata covariates with PCs,
#' and visualizes the top PCs as heatmaps: one for -log10(p-values), one for effect sizes.
#'
#' @param data_matrix numeric matrix or data.frame, rows = features, columns = samples.
#' @param metadata data.frame with sample metadata, must include `sample_id`.
#' @param covariates character vector of column names in `metadata` to test.
#' @param top_n_pcs number of top PCs to visualize (default 20).
#' @param remove_sex logical, remove features on sex chromosomes (row names starting with "chrX"/"chrY").
#' @param filter_sd_quantile numeric between 0 and 1, keep top X% most variable features (default 1 = keep all).
#' @param output_dir optional directory to save heatmaps (default NULL, skips saving).
#' @param filename_pval filename for p-value heatmap (default "PCA_covariate_heatmap_pval.png").
#' @param filename_effect filename for effect size heatmap (default "PCA_covariate_heatmap_effect.png").
#' @return A list with: `pca` object, `assoc_df_pval`, `assoc_df_effect`, `covariate_importance`, `heatmap_pval`, `heatmap_effect`.
#' @import ggplot2 dplyr tidyr
#' @export
pca_covariate_heatmap <- function(data_matrix, metadata, covariates,
                                  top_n_pcs = 20,
                                  remove_sex = TRUE,
                                  filter_sd_quantile = 1,
                                  output_dir = NULL,
                                  filename_pval = "PCA_covariate_heatmap_pval.png",
                                  filename_effect = "PCA_covariate_heatmap_effect.png") {

  stopifnot(all(c("sample_id", covariates) %in% colnames(metadata)))

  # ---- Remove sex chromosomes ----
  if (remove_sex) {
    data_matrix <- data_matrix[!grepl("^chrX|^chrY", rownames(data_matrix)), ]
  }

  # ---- Filter features by SD ----
  if (filter_sd_quantile < 1) {
    feature_sd <- apply(data_matrix, 1, sd, na.rm = TRUE)
    threshold <- quantile(feature_sd, filter_sd_quantile)
    data_matrix <- data_matrix[feature_sd >= threshold, , drop = FALSE]
  }

  # ---- PCA ----
  pca <- prcomp(t(data_matrix), center = TRUE, scale. = TRUE)

  # ---- Merge PCA scores with metadata ----
  pca_scores <- as.data.frame(pca$x)
  pca_scores$sample_id <- rownames(pca_scores)
  pca_meta <- merge(pca_scores, metadata, by = "sample_id")

  pcs <- pca_scores[, grep("^PC", colnames(pca_scores)), drop = FALSE]
  n_pcs_available <- ncol(pcs)
  n_pcs_use <- min(top_n_pcs, n_pcs_available)

  # ---- Association function returning both p-value & effect size ----
  test_assoc <- function(pc, annotation) {
    tryCatch({
      if (is.numeric(annotation)) {
        if (var(annotation, na.rm = TRUE) == 0) return(c(pval = NA_real_, effect = NA_real_))
        cor_res <- stats::cor.test(pc, annotation)
        # cor.test() returns estimate named "cor" or "estimate" depending on method; extract numeric
        effect_val <- as.numeric(cor_res$estimate)
        c(pval = as.numeric(cor_res$p.value), effect = effect_val)
      } else {
        annotation <- as.factor(annotation)
        if (nlevels(annotation) < 2) return(c(pval = NA_real_, effect = NA_real_))
        aov_res <- stats::aov(pc ~ annotation)
        ssum <- summary(aov_res)[[1]][["Sum Sq"]]
        if (length(ssum) < 2) return(c(pval = NA_real_, effect = NA_real_))
        eta2 <- ssum[1] / sum(ssum)
        pval <- summary(aov_res)[[1]][["Pr(>F)"]][1]
        c(pval = as.numeric(pval), effect = as.numeric(eta2))
      }
    }, error = function(e) c(pval = NA_real_, effect = NA_real_))
  }

  # ---- Compute associations for all covariates & PCs ----
  assoc_list <- lapply(covariates, function(cov) {
    # sapply returns a 2 x n_pcs matrix (rows: pval,effect), transpose to n_pcs x 2
    tmp <- sapply(seq_len(n_pcs_available), function(i) test_assoc(pcs[, i], pca_meta[[cov]]))
    if (is.null(dim(tmp))) {
      tmp <- matrix(tmp, nrow = 2)
    }
    tmp_t <- t(tmp)                # now n_pcs x 2
    colnames(tmp_t) <- rownames(tmp) # names from test_assoc: c("pval","effect")
    rownames(tmp_t) <- colnames(pcs)
    tmp_t
  })

  # Build matrices: rows = PCs, cols = covariates
  assoc_matrix_pval <- do.call(cbind, lapply(assoc_list, function(mat) mat[, "pval"]))
  assoc_matrix_effect <- do.call(cbind, lapply(assoc_list, function(mat) mat[, "effect"]))

  # set dimnames correctly: rows = PC names, cols = covariate names
  rownames(assoc_matrix_pval) <- colnames(pcs)
  rownames(assoc_matrix_effect) <- colnames(pcs)
  colnames(assoc_matrix_pval) <- covariates
  colnames(assoc_matrix_effect) <- covariates

  # If top_n_pcs requested > available, trim matrices to available
  if (n_pcs_use < nrow(assoc_matrix_pval)) {
    assoc_matrix_pval <- assoc_matrix_pval[seq_len(n_pcs_use), , drop = FALSE]
    assoc_matrix_effect <- assoc_matrix_effect[seq_len(n_pcs_use), , drop = FALSE]
  }

  # ---- Convert to tidy format ----
  tidy_assoc <- function(mat) {
    as.data.frame(mat) %>%
      tibble::rownames_to_column(var = "PC") %>%
      tidyr::pivot_longer(-PC, names_to = "covariate", values_to = "value")
  }

  assoc_df_pval <- tidy_assoc(assoc_matrix_pval) %>%
    dplyr::mutate(log10p = -log10(value))

  assoc_df_effect <- tidy_assoc(assoc_matrix_effect) %>%
    dplyr::mutate(effect = value)

  # ---- Order PCs & covariates (PC1 at bottom) ----
  pc_levels <- rownames(assoc_matrix_pval)
  # ensure order PC1..PCn if those names present, else keep current order
  if (all(grepl("^PC\\d+$", pc_levels))) {
    pc_numeric <- as.integer(gsub("PC", "", pc_levels))
    ord <- order(pc_numeric)
    pc_levels <- pc_levels[ord]
  }
  assoc_df_pval$PC <- factor(assoc_df_pval$PC, levels = rev(pc_levels))
  assoc_df_effect$PC <- factor(assoc_df_effect$PC, levels = rev(pc_levels))

  # Order covariates by significance (from pval matrix)
  cov_order <- assoc_df_pval %>%
    dplyr::group_by(covariate) %>%
    dplyr::summarize(max_log10p = max(log10p, na.rm = TRUE)) %>%
    dplyr::arrange(desc(max_log10p)) %>%
    dplyr::pull(covariate)

  assoc_df_pval$covariate <- factor(assoc_df_pval$covariate, levels = cov_order)
  assoc_df_effect$covariate <- factor(assoc_df_effect$covariate, levels = cov_order)

  # ---- Heatmaps ----
  heatmap_pval <- ggplot2::ggplot(assoc_df_pval, ggplot2::aes(x = covariate, y = PC, fill = log10p)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Strength of Covariate Associations with PCs\n(High -log10(p) = stronger association)",
      x = "Covariate",
      y = "Principal Component",
      fill = "-log10(p)"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Use absolute effect size for coloring
  assoc_df_effect <- assoc_df_effect %>% dplyr::mutate(effect = abs(effect))

  heatmap_effect <- ggplot2::ggplot(assoc_df_effect, ggplot2::aes(x = covariate, y = PC, fill = effect)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "magma", na.value = "grey90") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Magnitude of Covariate Associations with PCs\n(Effect size: |correlation| for numeric, eta² for categorical)",
      x = "Covariate",
      y = "Principal Component",
      fill = "Effect (|value|)"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # ---- Save if requested ----
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(output_dir, filename_pval), plot = heatmap_pval, width = 8, height = 6, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(output_dir, filename_effect), plot = heatmap_effect, width = 8, height = 6, dpi = 300, bg = "white")
  }

  # ---- Covariate importance summary (based on -log10 p) ----
  covariate_importance <- assoc_df_pval %>%
    dplyr::group_by(covariate) %>%
    dplyr::summarize(max_log10p = max(log10p, na.rm = TRUE)) %>%
    dplyr::arrange(desc(max_log10p))

  return(list(
    pca = pca,
    assoc_df_pval = assoc_df_pval,
    assoc_df_effect = assoc_df_effect,
    covariate_importance = covariate_importance,
    heatmap_pval = heatmap_pval,
    heatmap_effect = heatmap_effect
  ))
}
