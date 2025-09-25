
## File: R/analysis_wrappers.R
#' Wrapper to run PCA/LDA exploration, outlier detection and save plots
#'
#' Delegates to smaller internal functions: .run_pca, .plot_pca_save, .run_lda, .plot_lda_save, .run_techrep_pca
#' @param meth_united methylKit unified object (e.g. methylBase)
#' @param feature_name name string
#' @param meta sample metadata data.frame containing at least sample_id and relevant covariates
#' @param plot_path path to write plots
#' @param z_threshold numeric z-score threshold for flagging PC1/PC2 outliers
#' @param lda_class optional column name in meta for LDA class
#' @param keep_samples optional character vector of sample IDs to analyse (will reorganize to these before PCA)
#' @param color_by character: metadata column used for point color in PCA plots
#' @param shape_by character: metadata column used for point shape in PCA plots
#' @param label_by character: metadata column used for point labels in PCA plots
#' @param draw_by character: metadata column used to draw group ellipses
#' @param point_size numeric point size for PCA plots
#' @import methylKit
#' @import ggplot2
#' @export
count_exploration_wrapper <- function(meth_united, feature_name, meta, plot_path,
                                      z_threshold = 5, lda_class = NULL,
                                      keep_samples = NULL,
                                      color_by = "label", shape_by = "gender", label_by = NULL,
                                      draw_by = NULL, point_size = 3) {
  .ensure_dir(plot_path)
  .ensure_dir(file.path(plot_path, "PCA"))
  .ensure_dir(file.path(plot_path, "LDA"))

  # If user provided a subset of samples, reorganize to those up front (manual selection)
  meth_work <- meth_united
  if (!is.null(keep_samples)) {
    keep_samples0 <- intersect(keep_samples, meth_united@sample.ids)
    if (length(keep_samples0) == 0) stop("`keep_samples` does not match any sample IDs in the methyl object")
    meth_work <- methylKit::reorganize(meth_united, sample.ids = keep_samples0, treatment = rep(0, length(keep_samples0)))
  }

  # PCA on working set
  pca_main <- .run_pca(meth_work)
  pca_scores_main <- as.data.frame(pca_main$x)
  pca_scores_main$sample_id <- rownames(pca_scores_main)
  pca_meta_main <- merge(pca_scores_main, meta, by = "sample_id")

  # Detect outliers (flag only)
  outlier_samples <- .detect_outliers(pca_meta_main, z_threshold = z_threshold)
  message("Flagged ", length(outlier_samples), " outliers: ", paste(outlier_samples, collapse = ", "))

  # --- PCA plots ---
  p_all <- .plot_pca_save(pca_main, meta %>% dplyr::filter(sample_id %in% meth_work@sample.ids),
                          file.path(plot_path, "PCA", "PCASamples_all.png"),
                          title = paste0(feature_name, ": PCA (working set)"),
                          color_by = color_by, shape_by = shape_by, label_by = label_by, draw_by = draw_by, point_size = point_size)

  # PCA with flagged outliers highlighted
  meta_flagged <- (pca_scores_main %>% dplyr::select(sample_id)) %>% dplyr::left_join(meta, by = "sample_id")
  meta_flagged$outlier_flag <- meta_flagged$sample_id %in% outlier_samples
  meta_flagged$outlier_flag <- factor(meta_flagged$outlier_flag, levels = c("FALSE", "TRUE"))
  p_flagged <- .plot_pca_save(pca_main, meta_flagged %>% dplyr::filter(sample_id %in% meth_work@sample.ids),
                              file.path(plot_path, "PCA", "PCASamples_flagged.png"),
                              title = paste0(feature_name, ": PCA (outliers flagged)"),
                              color_by = "outlier_flag", shape_by = shape_by, label_by = label_by, draw_by = NULL, point_size = point_size)

  # Scree plot
  png(file.path(plot_path, "PCA", "PCASamples_scree.png"), width = 800, height = 600)
  methylKit::PCASamples(meth_work, screeplot = TRUE)
  dev.off()

  # LDA (if requested)
  lda_obj <- NULL
  if (!is.null(lda_class) && lda_class %in% names(meta)) {
    meta_for_lda <- meta %>% dplyr::filter(sample_id %in% meth_work@sample.ids)
    lda_obj <- .run_lda(pca_main, meta_for_lda, lda_class)
    .plot_lda_save(lda_obj, meta_for_lda, file.path(plot_path, "LDA", "LDA.png"),
                   title = paste0(feature_name, ": LDA (", lda_class, ")"),
                   color_by = lda_class, shape_by = shape_by, label_by = label_by)
  }

  return(list(
    pca_obj = pca_main,
    pca_meta = pca_scores_main,
    lda_obj = lda_obj,
    outliers = outlier_samples
  ))
}

#' @noRd
.run_pca <- function(meth_obj) {
  suppressWarnings(methylKit::PCASamples(meth_obj, obj.return = TRUE))
}
#' @noRd
.plot_pca_save <- function(pca_obj, meta, filename, title = "PCA",
                           color_by = "label", shape_by = "gender", label_by = NULL,
                           draw_by = NULL, point_size = 3) {
  p <- plot_prcomp_gg(
    pca_obj,
    metadata = meta,
    color_by = color_by,
    shape_by = shape_by,
    label_by = label_by,
    draw_by = draw_by,
    point_size = point_size,
    title = title
  )
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