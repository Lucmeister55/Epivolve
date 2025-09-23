
## File: R/analysis_wrappers.R
#' Wrapper to run PCA/LDA exploration, outlier detection and save plots
#'
#' Delegates to smaller internal functions: .run_pca, .plot_pca_save, .run_lda, .plot_lda_save, .run_techrep_pca
#' @param meth_united methylKit unified object (e.g. methylBase)
#' @param feature_name name string
#' @param meta sample metadata data.frame containing at least sample_id and relevant covariates
#' @param plot_path path to write plots
#' @param tech_reps data.frame with columns rep1, rep2 listing technical replicate pairs (optional)
#' @param z_threshold numeric z-score threshold for flagging PC1/PC2 outliers
#' @param lda_class optional column name in meta for LDA class
#' @import methylKit
#' @import ggplot2
#' @export
count_exploration_wrapper <- function(meth_united, feature_name, meta, plot_path,
                                      tech_reps = NULL, z_threshold = 5, lda_class = NULL) {
  .ensure_dir(plot_path)
  .ensure_dir(file.path(plot_path, "PCA"))
  .ensure_dir(file.path(plot_path, "LDA"))

  pca_main <- .run_pca(meth_united)
  pca_scores_main <- as.data.frame(pca_main$x)
  pca_scores_main$sample_id <- rownames(pca_scores_main)
  pca_meta_main <- merge(pca_scores_main, meta, by = "sample_id")

  outlier_samples <- .detect_outliers(pca_meta_main, z_threshold = z_threshold)
  message("Flagged ", length(outlier_samples), " outliers: ", paste(outlier_samples, collapse = ", "))

  p_all <- .plot_pca_save(pca_main, meta, file.path(plot_path, "PCA", "PCASamples_all.png"), title = paste0(feature_name, ": PCA (all samples)"))

  keep_samples <- setdiff(meth_united@sample.ids, outlier_samples)
  if (length(keep_samples) == 0) stop("All samples flagged as outliers")
  meth_clean <- methylKit::reorganize(meth_united, sample.ids = keep_samples, treatment = rep(0, length(keep_samples)))

  pca_clean <- .run_pca(meth_clean)
  p_clean <- .plot_pca_save(pca_clean, meta %>% dplyr::filter(sample_id %in% keep_samples), file.path(plot_path, "PCA", "PCASamples_no_outliers.png"), title = paste0(feature_name, ": PCA (no outliers)"))

  # Scree
  png(file.path(plot_path, "PCA", "PCASamples_scree.png"), width = 800, height = 600)
  methylKit::PCASamples(meth_united, screeplot = TRUE)
  dev.off()

  lda_obj <- NULL
  if (!is.null(lda_class) && lda_class %in% names(meta)) {
    lda_obj <- .run_lda(pca_clean, meta %>% dplyr::filter(sample_id %in% keep_samples), lda_class)
    .plot_lda_save(lda_obj, meta %>% dplyr::filter(sample_id %in% keep_samples), file.path(plot_path, "LDA", "LDA.png"), title = paste0(feature_name, ": LDA (", lda_class, ")"))
  }

  if (!is.null(tech_reps) && nrow(tech_reps) > 0) {
    .run_techrep_pca(meth_clean, meta, tech_reps, file.path(plot_path, "PCA", "PCASamples_techReps.png"))
  }

  return(list(pca_obj = pca_main, pca_meta = pca_scores_main, lda_obj = lda_obj, outliers = outlier_samples, meth_clean = meth_clean))
}

#' @noRd
.run_pca <- function(meth_obj) {
  suppressWarnings(methylKit::PCASamples(meth_obj, obj.return = TRUE))
}
#' @noRd
.plot_pca_save <- function(pca_obj, meta, filename, title = "PCA") {
  p <- plot_prcomp_gg(pca_obj, metadata = meta, color_by = "label", shape_by = "gender", draw_by = "label", point_size = 3, title = title)
  ggplot2::ggsave(filename = filename, plot = p, width = 10, height = 8)
  invisible(p)
}
#' @noRd
.detect_outliers <- function(pca_meta, z_threshold = 5) {
  z_pc1 <- scale(pca_meta$PC1)
  z_pc2 <- scale(pca_meta$PC2)
  pca_meta$outlier <- (abs(z_pc1) > z_threshold) | (abs(z_pc2) > z_threshold)
  pca_meta$sample_id[pca_meta$outlier]
}
#' @noRd
.run_lda <- function(pca_obj, meta, lda_class) {
  var_exp <- summary(pca_obj)$importance[2, ]
  k <- which(cumsum(var_exp) >= 0.95)[1]
  lda_data <- as.data.frame(pca_obj$x[, 1:k, drop = FALSE])
  lda_data <- lda_data[, apply(lda_data, 2, function(col) length(unique(col)) > 1), drop = FALSE]
  lda_df <- cbind(meta[lda_class], lda_data)
  colnames(lda_df)[1] <- lda_class
  MASS::lda(as.formula(paste(lda_class, "~ .")), data = lda_df)
}
#' @noRd
.plot_lda_save <- function(lda_obj, meta, filename, title = "LDA") {
  p <- plot_lda_gg(lda_obj, metadata = meta, class_col = names(meta)[1], color_by = names(meta)[1], shape_by = "gender", title = title)
  ggplot2::ggsave(filename = filename, plot = p, width = 10, height = 8)
  invisible(p)
}
#' @noRd
.run_techrep_pca <- function(meth_clean, meta, tech_reps, filename) {
  sample_ids_tech <- unique(c(tech_reps$rep1, tech_reps$rep2))
  sample_ids_tech <- sample_ids_tech[sample_ids_tech %in% meth_clean@sample.ids]
  if (length(sample_ids_tech) > 0) {
    meth_tech <- methylKit::reorganize(meth_clean, sample.ids = sample_ids_tech, treatment = rep(0, length(sample_ids_tech)))
    p_tech <- plot_prcomp_gg(methylKit::PCASamples(meth_tech, obj.return = TRUE), metadata = meta %>% dplyr::filter(sample_id %in% sample_ids_tech), color_by = "donor_id", label_by = "sample_id", shape_by = "label", title = "PCA of Technical Replicates", point_size = 3)
    ggplot2::ggsave(filename = filename, plot = p_tech, width = 8, height = 6)
  }
}

#' LDA plotting helper using ggplot2
#' @import ggplot2
#' @export
plot_lda_gg <- function(lda, metadata, class_col, dim_axes = c(1,2), color_by = NULL, shape_by = NULL, label_by = NULL, label_size = 3, title = "LDA plot", point_size = 3, alpha = 0.8, theme = ggplot2::theme_minimal()) {
  if (!inherits(lda, "lda")) stop("`lda` must be an object from MASS::lda")
  lda_scores <- as.data.frame(predict(lda)$x)
  colnames(lda_scores) <- paste0("LD", seq_len(ncol(lda_scores)))
  lda_scores$sample_id <- metadata$sample_id
  lda_scores <- merge(lda_scores, metadata, by = "sample_id")
  xlab <- sprintf("LD%d", dim_axes[1])
  ylab <- sprintf("LD%d", dim_axes[2])
  aes_map <- ggplot2::aes_string(x = colnames(lda_scores)[dim_axes[1] + 1], y = colnames(lda_scores)[dim_axes[2] + 1])
  if (!is.null(color_by)) aes_map$colour <- as.name(color_by)
  if (!is.null(shape_by)) aes_map$shape <- as.name(shape_by)
  p <- ggplot2::ggplot(lda_scores, aes_map) + ggplot2::geom_point(size = point_size, alpha = alpha) + ggplot2::labs(title = title, x = xlab, y = ylab, colour = color_by, shape = shape_by) + theme
  if (!is.null(label_by)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(label = label_by), size = label_size)
    } else {
      p <- p + ggplot2::geom_text(ggplot2::aes_string(label = label_by), size = label_size, vjust = -1)
    }
  }
  return(p)
}

#' PCA plotting helper using ggplot2
#' @import ggplot2
#' @export
plot_prcomp_gg <- function(pca, metadata = NULL, pc_axes = c(1,2), color_by = NULL, shape_by = NULL, size_by = NULL, label_by = NULL, label_size = 3, repel = TRUE, draw_by = NULL, ellipse_fill = FALSE, ellipse_alpha = 0.15, ellipse_size = 1, title = "PCA plot", point_size = 3, alpha = 0.8, theme = ggplot2::theme_minimal()) {
  if (!inherits(pca, "prcomp")) stop("`pca` must be a prcomp object")
  if (length(pc_axes) != 2) stop("`pc_axes` must be length 2 (e.g. c(1,2))")
  scores <- as.data.frame(pca$x[, pc_axes, drop = FALSE])
  colnames(scores) <- paste0("PC", pc_axes)
  var_exp <- summary(pca)$importance[2, pc_axes] * 100
  xlab <- sprintf("PC%d (%.1f%%)", pc_axes[1], var_exp[1])
  ylab <- sprintf("PC%d (%.1f%%)", pc_axes[2], var_exp[2])
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) stop("`metadata` must be a data frame")
    if (nrow(metadata) != nrow(scores)) stop("`metadata` rows must match `pca$x` rows")
    overlap <- intersect(names(metadata), names(scores))
    if (length(overlap) > 0) {
      warning("Metadata column names overlap with PC names: ", paste(overlap, collapse = ", "))
      metadata <- metadata %>% dplyr::rename_at(dplyr::vars(dplyr::one_of(overlap)), ~ paste0(., "_meta"))
    }
    scores <- cbind(scores, metadata)
  }
  for (aes_col in c(draw_by, color_by, shape_by, size_by, label_by)) {
    if (!is.null(aes_col) && !(aes_col %in% names(scores))) stop(paste0("`", aes_col, "` not found in scores/metadata"))
  }
  p <- ggplot2::ggplot(scores, ggplot2::aes_string(x = colnames(scores)[1], y = colnames(scores)[2]))
  if (!is.null(draw_by)) {
    if (ellipse_fill) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(fill = draw_by, group = draw_by), geom = "polygon", alpha = ellipse_alpha, show.legend = FALSE)
    } else {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(color = draw_by, group = draw_by), size = ellipse_size, show.legend = FALSE)
    }
  }
  aes_map <- list()
  if (!is.null(color_by)) aes_map$colour <- as.name(color_by)
  if (!is.null(shape_by)) aes_map$shape <- as.name(shape_by)
  if (!is.null(size_by)) aes_map$size <- as.name(size_by)
  geom_call <- do.call(ggplot2::geom_point, c(list(mapping = do.call(ggplot2::aes_string, aes_map), alpha = alpha), list(size = point_size)))
  p <- p + geom_call
  if (!is.null(label_by)) {
    if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(label = label_by), size = label_size, max.overlaps = Inf)
    } else {
      p <- p + ggplot2::geom_text(ggplot2::aes_string(label = label_by), size = label_size, vjust = -1, show.legend = FALSE)
      if (repel) message("ggrepel not installed — falling back to geom_text()")
    }
  }
  p <- p + ggplot2::labs(title = title, x = xlab, y = ylab, colour = color_by, shape = shape_by, size = size_by) + theme
  return(p)
}