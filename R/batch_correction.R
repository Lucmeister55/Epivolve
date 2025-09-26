# batch_correction.R
# Submodules for detecting, quantifying and removing batch / covariate effects
#
# Usage: source() this file or add to your package R/ directory.
#
# NOTE: functions try to avoid assuming global state. Pass Mval and metadata explicitly.

###########
# Helpers #
###########

#' Run PCA on samples (features x samples -> prcomp on samples)
#'
#' @param Mval numeric matrix, rows = features, cols = samples
#' @return list with pca (prcomp) and pca_scores (data.frame with sample_id)
#' @importFrom stats prcomp
#' @export
compute_sample_pca <- function(Mval, center = TRUE, scale. = TRUE) {
  stopifnot(is.matrix(Mval) || is.data.frame(Mval))
  Mval <- as.matrix(Mval)
  pca <- stats::prcomp(t(Mval), center = center, scale. = scale.)
  pca_scores <- as.data.frame(pca$x, stringsAsFactors = FALSE)
  pca_scores$sample_id <- rownames(pca_scores)
  list(pca = pca, pca_scores = pca_scores)
}


#' Score covariates by weighted R^2 and Fisher-combined p-value across PCs.
#'
#' @param pca_scores data.frame with PC columns named "PC1","PC2",... and column sample_id
#' @param metadata data.frame with sample_id
#' @param covariates character vector of covariate names to test
#' @param top_n_pcs integer number of PCs to consider (use available if fewer)
#' @param min_samples minimum complete cases required to compute stats
#' @return data.frame(covariate, weighted_R2, fisher_p, n)
#' @export
rank_covariates_by_pcs <- function(pca_scores, metadata, covariates,
                                   top_n_pcs = 20, min_samples = 4) {
  stopifnot("sample_id" %in% colnames(pca_scores))
  # ensure metadata aligned to pca_scores order
  metadata <- metadata[match(pca_scores$sample_id, metadata$sample_id), , drop = FALSE]

  pc_cols <- grep("^PC\\d+$", colnames(pca_scores), value = TRUE)
  if (length(pc_cols) == 0) stop("No PC columns found in pca_scores")
  n_pcs_available <- length(pc_cols)
  top_n <- min(top_n_pcs, n_pcs_available)
  pcs_names <- pc_cols[seq_len(top_n)]
  varw_all <- (summary(stats::prcomp(t(matrix(0, nrow = 1, ncol = 1))))$importance) # placeholder to avoid NOTE
  # We'll compute var weights from a temporary prcomp if available: user should pass pca for exact weights -
  # but prefer to compute from pca_scores if present (approx).
  # Simpler: estimate variance weights from pca_scores columns variance:
  var_vals <- sapply(pca_scores[, pcs_names, drop = FALSE], function(x) var(x, na.rm = TRUE))
  varw <- var_vals / sum(var_vals)

  score_list <- lapply(covariates, function(v) {
    x <- metadata[[v]]
    ok <- stats::complete.cases(x, pca_scores[, pcs_names, drop = FALSE])
    n_ok <- sum(ok)
    if (n_ok < min_samples) {
      return(data.frame(covariate = v,
                        weighted_R2 = NA_real_,
                        fisher_p = NA_real_,
                        n = n_ok,
                        stringsAsFactors = FALSE))
    }
    if (is.character(x)) x <- factor(x)
    Rs <- sapply(pcs_names, function(pc) {
      fit <- tryCatch(stats::lm(pca_scores[ok, pc] ~ x[ok]), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      summary(fit)$r.squared
    }, USE.NAMES = FALSE)
    ps <- sapply(pcs_names, function(pc) {
      fit <- tryCatch(stats::lm(pca_scores[ok, pc] ~ x[ok]), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      if (is.factor(x) || is.character(x)) {
        a <- tryCatch(stats::anova(fit), error = function(e) NULL)
        if (is.null(a)) return(NA_real_)
        a[["Pr(>F)"]][1]
      } else {
        s <- summary(fit)
        if (nrow(s$coefficients) >= 2) s$coefficients[2, 4] else NA_real_
      }
    }, USE.NAMES = FALSE)
    ps2 <- ps[!is.na(ps) & ps > 0]
    fisher_p <- if (length(ps2) > 0) stats::pchisq(-2 * sum(log(ps2)), df = 2 * length(ps2), lower.tail = FALSE) else NA_real_
    weighted_R2 <- sum(Rs * varw, na.rm = TRUE)
    data.frame(covariate = v, weighted_R2 = weighted_R2, fisher_p = fisher_p, n = n_ok, stringsAsFactors = FALSE)
  })
  batch_ranking <- do.call(rbind, score_list)
  batch_ranking[order(-batch_ranking$weighted_R2), , drop = FALSE]
}


#' Build covariate model matrix for correction
#'
#' @param metadata data.frame already ordered to match Mval columns (sample_id order)
#' @param covariates character vector of covariates to include
#' @return matrix (samples x covariate-columns)
#' @export
build_covariate_matrix <- function(metadata, covariates) {
  stopifnot("sample_id" %in% colnames(metadata))
  mats <- lapply(covariates, function(v) {
    x <- metadata[[v]]
    if (is.character(x)) x <- factor(x)
    if (is.factor(x)) {
      mm <- stats::model.matrix(~ x - 1)
      colnames(mm) <- paste0(v, "_", sub("^x", "", colnames(mm)))
      mm
    } else {
      vec <- as.numeric(x)
      mat <- matrix(scale(vec, center = TRUE, scale = FALSE), ncol = 1)
      colnames(mat) <- v
      mat
    }
  })
  cov_mat <- do.call(cbind, mats)
  cov_mat
}


#' Robust test association for a single PC vs annotation
#'
#' Returns named vector c(pval=..., effect=...)
#' For numeric: pval from cor.test, effect = correlation estimate
#' For factor: pval from aov, effect = eta^2 (proportion explained)
#' @export
test_assoc_pc_annotation <- function(pc, annotation) {
  tryCatch({
    if (is.numeric(annotation)) {
      if (var(annotation, na.rm = TRUE) == 0) return(c(pval = NA_real_, effect = NA_real_))
      res <- stats::cor.test(pc, annotation)
      effect <- as.numeric(res$estimate)
      c(pval = as.numeric(res$p.value), effect = effect)
    } else {
      annotation <- as.factor(annotation)
      if (nlevels(annotation) < 2) return(c(pval = NA_real_, effect = NA_real_))
      fit <- stats::aov(pc ~ annotation)
      ss <- summary(fit)[[1]][["Sum Sq"]]
      if (length(ss) < 2) return(c(pval = NA_real_, effect = NA_real_))
      eta2 <- ss[1] / sum(ss)
      pval <- summary(fit)[[1]][["Pr(>F)"]][1]
      c(pval = as.numeric(pval), effect = as.numeric(eta2))
    }
  }, error = function(e) c(pval = NA_real_, effect = NA_real_))
}


#' Apply batch correction using limma or linear regression fallback
#'
#' @param Mval features x samples matrix
#' @param cov_mat matrix (samples x covariate-columns)
#' @param use_limma logical
#' @return corrected matrix (same dimnames)
#' @export
apply_batch_correction <- function(Mval, cov_mat, use_limma = TRUE) {
  Mval <- as.matrix(Mval)
  # cov_mat should have rows in same order as columns of Mval
  if (nrow(cov_mat) != ncol(Mval)) {
    stop("Rows of cov_mat must match columns of Mval (samples order)")
  }
  if (use_limma && requireNamespace("limma", quietly = TRUE)) {
    corrected <- limma::removeBatchEffect(Mval, covariates = cov_mat)
  } else {
    X <- cbind(1, cov_mat)
    qrX <- qr(X)
    corrected <- t(apply(Mval, 1, function(y) {
      b <- tryCatch(qr.coef(qrX, y), error = function(e) rep(NA_real_, ncol(X)))
      if (any(is.na(b))) return(rep(NA_real_, length(y)))
      as.numeric(y - X %*% b)
    }))
    colnames(corrected) <- colnames(Mval)
    rownames(corrected) <- rownames(Mval)
  }
  corrected
}


#' Produce before/after plots for a given covariate
#'
#' @param Mval original matrix (features x samples)
#' @param Mval_corrected corrected matrix (features x samples)
#' @param metadata data.frame matched to columns order, must contain sample_id and covariate
#' @param covariate name of covariate to plot
#' @return list of ggplot objects (mean_pre, mean_post, pca_pre, pca_post) or NULL if ggplot2 not available
#' @export
plot_before_after_covariate <- function(Mval, Mval_corrected, metadata, covariate) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping plots for ", covariate)
    return(NULL)
  }
  # compute per-sample mean M-values
  mean_pre <- colMeans(Mval, na.rm = TRUE)
  mean_post <- colMeans(Mval_corrected, na.rm = TRUE)
  covx <- metadata[[covariate]]
  df_pre <- data.frame(sample_id = colnames(Mval), mean_M = mean_pre, covariate = covx, stringsAsFactors = FALSE)
  df_post <- data.frame(sample_id = colnames(Mval_corrected), mean_M = mean_post, covariate = covx, stringsAsFactors = FALSE)

  isfac <- is.factor(covx) || is.character(covx)
  if (isfac) {
    df_pre$covariate <- factor(df_pre$covariate, exclude = NULL)
    df_post$covariate <- factor(df_post$covariate, exclude = NULL)
    p_pre <- ggplot2::ggplot(df_pre, ggplot2::aes(x = covariate, y = mean_M)) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "salmon", alpha = 0.7) +
      ggplot2::labs(title = paste("Before correction —", covariate), y = "Mean M-value", x = covariate) +
      ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    p_post <- ggplot2::ggplot(df_post, ggplot2::aes(x = covariate, y = mean_M)) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "skyblue", alpha = 0.7) +
      ggplot2::labs(title = paste("After correction —", covariate), y = "Mean M-value", x = covariate) +
      ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    df_pre$covariate_num <- as.numeric(df_pre$covariate)
    df_post$covariate_num <- as.numeric(df_post$covariate)
    p_pre <- ggplot2::ggplot(df_pre, ggplot2::aes(x = covariate_num, y = mean_M)) +
      ggplot2::geom_point(alpha = 0.6) + ggplot2::geom_smooth(method = "loess", se = FALSE) +
      ggplot2::labs(title = paste("Before correction —", covariate), y = "Mean M-value", x = covariate) +
      ggplot2::theme_bw()
    p_post <- ggplot2::ggplot(df_post, ggplot2::aes(x = covariate_num, y = mean_M)) +
      ggplot2::geom_point(alpha = 0.6) + ggplot2::geom_smooth(method = "loess", se = FALSE) +
      ggplot2::labs(title = paste("After correction —", covariate), y = "Mean M-value", x = covariate) +
      ggplot2::theme_bw()
  }

  # PCA scatter before/after colored by covariate
  pca_before <- compute_sample_pca(Mval)$pca
  pca_after <- compute_sample_pca(Mval_corrected)$pca
  df_pb <- as.data.frame(pca_before$x[, 1:2, drop = FALSE]); df_pb$sample_id <- rownames(df_pb)
  df_pa <- as.data.frame(pca_after$x[, 1:2, drop = FALSE]); df_pa$sample_id <- rownames(df_pa)
  df_pb <- merge(df_pb, metadata[, c("sample_id", covariate)], by = "sample_id", all.x = TRUE)
  df_pa <- merge(df_pa, metadata[, c("sample_id", covariate)], by = "sample_id", all.x = TRUE)
  colnames(df_pb)[2:3] <- c("PC1", "PC2"); colnames(df_pa)[2:3] <- c("PC1", "PC2")
  p_pcb <- ggplot2::ggplot(df_pb, ggplot2::aes(x = PC1, y = PC2, color = .data[[covariate]])) +
    ggplot2::geom_point(alpha = 0.8) + ggplot2::labs(title = paste("PCA before —", covariate)) + ggplot2::theme_minimal()
  p_pca <- ggplot2::ggplot(df_pa, ggplot2::aes(x = PC1, y = PC2, color = .data[[covariate]])) +
    ggplot2::geom_point(alpha = 0.8) + ggplot2::labs(title = paste("PCA after —", covariate)) + ggplot2::theme_minimal()

  list(mean_pre = p_pre, mean_post = p_post, pca_before = p_pcb, pca_after = p_pca)
}


####################
# High-level module #
####################

#' Detect and correct batch / covariate effects (high-level)
#'
#' This orchestrates PCA → ranking → covariate selection → correction → plotting.
#' It calls the smaller submodules above.
#'
#' @param Mval numeric matrix features x samples
#' @param metadata data.frame with sample_id column ordered or not (will be matched)
#' @param covariates character vector of covariate names (NULL -> all except sample_id)
#' @param tech_covars character vector of covariates considered "technical" (NULL -> use all)
#' @param top_n_pcs integer number of top PCs to score (default 20)
#' @param fisher_p_threshold numeric threshold for Fisher combined p to consider correction (default 0.05)
#' @param min_samples minimal number of samples with data to score a covariate
#' @param use_limma logical to use limma for correction if available
#' @param output_dir optional directory for saving results / plots
#' @param save_corrected logical to save corrected matrix
#' @param plot_before_after logical to create and save before/after plots
#' @param verbose logical
#' @return list with Mval_corrected, batch_ranking, corrected_covariates, pca_before, pca_after, plots
#' @export
detect_and_correct_batch <- function(Mval, metadata,
                                     covariates = NULL,
                                     tech_covars = NULL,
                                     top_n_pcs = 20,
                                     fisher_p_threshold = 0.05,
                                     min_samples = 4,
                                     use_limma = TRUE,
                                     output_dir = NULL,
                                     save_corrected = TRUE,
                                     plot_before_after = TRUE,
                                     verbose = TRUE) {
  # basic checks & reorder metadata to Mval columns
  if (!is.matrix(Mval) && !is.data.frame(Mval)) stop("Mval must be a matrix or data.frame")
  Mval <- as.matrix(Mval)
  if (!"sample_id" %in% colnames(metadata)) stop("metadata must contain 'sample_id'")
  stopifnot(!is.null(colnames(Mval)))
  if (!all(colnames(Mval) %in% metadata$sample_id)) stop("Not all Mval columns found in metadata$sample_id")
  metadata <- metadata[match(colnames(Mval), metadata$sample_id), , drop = FALSE]

  if (is.null(covariates)) covariates <- setdiff(colnames(metadata), "sample_id")
  if (length(covariates) == 0) stop("No covariates to evaluate")

  # PCA
  if (verbose) message("Computing PCA...")
  pca_res <- compute_sample_pca(Mval)
  pca <- pca_res$pca
  pca_scores <- pca_res$pca_scores

  # Rank covariates
  if (verbose) message("Ranking covariates...")
  batch_ranking <- rank_covariates_by_pcs(pca_scores = pca_scores,
                                          metadata = metadata,
                                          covariates = covariates,
                                          top_n_pcs = top_n_pcs,
                                          min_samples = min_samples)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(batch_ranking, file.path(output_dir, "BatchCorrection", "covariate_batch_ranking.csv"), row.names = FALSE)
  }
  if (verbose) {
    message("Top covariates (by weighted_R2):")
    print(utils::head(batch_ranking, 10))
  }

  # select candidates based on tech_covars and fisher p
  if (is.null(tech_covars)) {
    candidates <- batch_ranking$covariate
  } else {
    candidates <- intersect(tech_covars, batch_ranking$covariate)
  }
  selected <- candidates[batch_ranking$fisher_p[match(candidates, batch_ranking$covariate)] < fisher_p_threshold]
  # require complete values
  selected <- selected[sapply(selected, function(v) all(!is.na(metadata[[v]])))]
  if (length(selected) == 0) {
    message("No covariates selected for correction (none passed threshold/completeness). Returning original matrix.")
    return(list(Mval_corrected = Mval,
                batch_ranking = batch_ranking,
                corrected_covariates = character(0),
                pca_before = pca,
                pca_after = pca,
                plots = list()))
  }

  if (verbose) message("Covariates selected for correction: ", paste(selected, collapse = ", "))

  # Build covariate matrix (rows must match Mval columns order)
  cov_mat <- build_covariate_matrix(metadata = metadata, covariates = selected)
  if (nrow(cov_mat) != ncol(Mval)) stop("covariate matrix row count mismatch with Mval columns")

  # Apply correction
  if (verbose) message("Applying batch correction...")
  Mval_corrected <- apply_batch_correction(Mval = Mval, cov_mat = cov_mat, use_limma = use_limma)

  # Save corrected and summary info
  if (!is.null(output_dir) && save_corrected) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(batch_ranking[match(selected, batch_ranking$covariate), , drop = FALSE],
                     file.path(output_dir, "BatchCorrection", "covariates_corrected_info.csv"), row.names = FALSE)
  }

  # Generate before/after plots
  plots_list <- list()
  if (plot_before_after && length(selected) > 0) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not available; skipping plotting.")
    } else {
      for (v in selected) {
        pl <- plot_before_after_covariate(Mval = Mval, Mval_corrected = Mval_corrected, metadata = metadata, covariate = v)
        plots_list[[v]] <- pl
        if (!is.null(output_dir) && !is.null(pl)) {
          subdir <- file.path(output_dir, "BatchCorrection")
          dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
          # save each plot element separately (user can combine later)
          if (!is.null(pl$mean_pre)) ggplot2::ggsave(file.path(subdir, paste0("meanM_by_", v, "_before.png")), plot = pl$mean_pre, width = 6, height = 4, dpi = 300)
          if (!is.null(pl$mean_post)) ggplot2::ggsave(file.path(subdir, paste0("meanM_by_", v, "_after.png")),  plot = pl$mean_post, width = 6, height = 4, dpi = 300)
          if (!is.null(pl$pca_before)) ggplot2::ggsave(file.path(subdir, paste0("PCA_by_", v, "_before.png")), plot = pl$pca_before, width = 6, height = 5, dpi = 300)
          if (!is.null(pl$pca_after))  ggplot2::ggsave(file.path(subdir, paste0("PCA_by_", v, "_after.png")),  plot = pl$pca_after, width = 6, height = 5, dpi = 300)
        }
      }
    }
  }

  pca_after <- compute_sample_pca(Mval_corrected)$pca

  list(Mval_corrected = Mval_corrected,
       batch_ranking = batch_ranking,
       corrected_covariates = selected,
       pca_before = pca,
       pca_after = pca_after,
       plots = plots_list)
}
