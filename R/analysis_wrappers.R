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
#' @param detect_outliers logical; if FALSE, skip outlier detection entirely (no flagged outliers)
#' @param save_outlier_plot logical; if FALSE, don't create the PCA plot that excludes flagged outliers
#' @param run_id optional short identifier to append to output filenames (e.g. "v2", "noBatchCorrect")
#' @param return_objects logical; if TRUE, return PCA and LDA objects and outlier list
#' @param plot_scree logical; if TRUE, create scree plot
#' @import methylKit
#' @import ggplot2
#' @export
count_exploration_wrapper <- function(meth_united, feature_name, meta, plot_path,
                                      z_threshold = 5, lda_class = NULL,
                                      keep_samples = NULL,
                                      color_by = "label", shape_by = "gender", label_by = NULL,
                                      draw_by = NULL, point_size = 3,
                                      detect_outliers = TRUE,         # NEW: whether to run detection
                                      save_outlier_plot = TRUE,       # NEW: whether to save the "no outliers" PCA plot
                                      run_id = NULL, return_objects = TRUE,
                                      plot_scree = TRUE) {                # NEW: identifier appended to filenames

  message("Starting exploration for feature: ", feature_name)
  .ensure_dir(plot_path)
  .ensure_dir(file.path(plot_path, "PCA"))
  .ensure_dir(file.path(plot_path, "LDA"))

  # helper for file names
  id_suffix <- if (!is.null(run_id) && nzchar(run_id)) paste0("_", run_id) else ""

  # Manual sample selection
  meth_work <- meth_united
  if (!is.null(keep_samples)) {
    message("Applying manual sample selection...")
    keep_samples0 <- intersect(keep_samples, meth_united@sample.ids)
    if (length(keep_samples0) == 0) stop("`keep_samples` does not match any sample IDs in the methyl object")
    meth_work <- methylKit::reorganize(meth_united, sample.ids = keep_samples0, treatment = rep(0, length(keep_samples0)))
  }

  # PCA
  message("Running PCA on working set...")
  pca_main <- .run_pca(meth_work)
  pca_scores_main <- as.data.frame(pca_main$x)
  pca_scores_main$sample_id <- rownames(pca_main$x)

  # Align metadata (keep only samples used in PCA)
  meta_for_plot <- meta %>% dplyr::filter(sample_id %in% rownames(pca_main$x))
  # Align metadata rows with pca rownames (important!)
  meta_for_plot <- meta_for_plot[match(rownames(pca_main$x), meta_for_plot$sample_id), , drop = FALSE]

  # Check aesthetics
  for (aes_col in c(color_by, shape_by, label_by, draw_by)) {
    if (!is.null(aes_col) && !(aes_col %in% colnames(meta_for_plot))) {
      stop(paste0("Column '", aes_col, "' not found in metadata for PCA plot"))
    }
  }

  # Outlier detection (optional)
  outlier_samples <- character(0)
  if (isTRUE(detect_outliers)) {
    message("Detecting outliers...")
    pca_meta_main <- merge(pca_scores_main, meta_for_plot, by = "sample_id")
    outlier_samples <- .detect_outliers(pca_meta_main, z_threshold = z_threshold)
    message("Flagged ", length(outlier_samples), " outliers: ", paste(outlier_samples, collapse = ", "))
  } else {
    message("Skipping outlier detection (detect_outliers = FALSE).")
  }

  # PCA plot with all samples
  message("Creating PCA plot for all samples...")
  p_all <- .plot_pca_save(pca_main,
                          meta_for_plot,
                          file.path(plot_path, "PCA", paste0("PCASamples_all", id_suffix, ".png")),
                          title = paste0(feature_name, ": PCA (working set)"),
                          color_by = color_by, shape_by = shape_by, label_by = label_by,
                          draw_by = draw_by, point_size = point_size)

  # PCA plot without outliers (optional)
  if (isTRUE(save_outlier_plot) && length(outlier_samples) > 0) {
    message("Creating PCA plot excluding flagged outliers...")
    keep_no_out <- setdiff(rownames(pca_main$x), outlier_samples)
    if (length(keep_no_out) > 1) {
      meth_no_out <- methylKit::reorganize(meth_work, sample.ids = keep_no_out, treatment = rep(0, length(keep_no_out)))
      pca_no_out <- .run_pca(meth_no_out)
      meta_no_out <- meta %>% dplyr::filter(sample_id %in% rownames(pca_no_out$x))
      # align meta order for the no-outlier PCA as well
      meta_no_out <- meta_no_out[match(rownames(pca_no_out$x), meta_no_out$sample_id), , drop = FALSE]

      p_no <- .plot_pca_save(pca_no_out,
                             meta_no_out,
                             file.path(plot_path, "PCA", paste0("PCASamples_no_outliers", id_suffix, ".png")),
                             title = paste0(feature_name, ": PCA (no flagged outliers)"),
                             color_by = color_by, shape_by = shape_by, label_by = label_by,
                             draw_by = draw_by, point_size = point_size)
    } else {
      message("Not enough samples remaining to compute PCA plot without outliers.")
    }
  } else if (!isTRUE(save_outlier_plot) && length(outlier_samples) > 0) {
    message("Outliers were flagged but save_outlier_plot = FALSE, skipping the 'no outliers' PCA plot.")
  }

  # Scree plot (keep name consistent with run_id)
  if (plot_scree) {
    message("Generating scree plot...")
    png(file.path(plot_path, "PCA", paste0("PCASamples_scree", id_suffix, ".png")),
        width = 800, height = 600)
    invisible(capture.output(
      methylKit::PCASamples(meth_work, screeplot = TRUE)
    ))
    dev.off()
  }
  # Optional LDA
  lda_obj <- NULL
  if (!is.null(lda_class) && lda_class %in% colnames(meta_for_plot)) {
    message("Running LDA using class: ", lda_class)
    lda_obj <- .run_lda(pca_main, meta_for_plot, lda_class)
    .plot_lda_save(
      lda_obj,
      meta_for_plot,
      lda_class,  # pass the actual LDA class here
      file.path(plot_path, "LDA", paste0("LDA", id_suffix, ".png")),
      title = paste0(feature_name, ": LDA (", lda_class, ")")
    )
  } else if (!is.null(lda_class)) {
    message("Requested lda_class '", lda_class, "' not found in plotting metadata; skipping LDA.")
  }

  message("Exploration finished for feature: ", feature_name)

  if (detect_outliers) {
    pca <- pca_no_out  # for returning
  } else {
    pca <- pca_main
  }

  # prepare result but return invisibly so nothing is printed at top-level
  result_list <- list(
    pca_obj = pca,
    pca_meta = pca_scores_main,
    lda_obj = lda_obj,
    outliers = outlier_samples
  )

  if (isTRUE(return_objects)) {
    invisible(result_list)   # returned invisibly (assignable with <- but won't auto-print)
  } else {
    invisible(NULL)          # no object returned, no printing
  }
}
#' @noRd
.run_pca <- function(meth_obj) {
  res <- NULL
  suppressWarnings(
    capture.output({
      # Open a null device to swallow plots
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      
      res <- methylKit::PCASamples(meth_obj, obj.return = TRUE)
    })
  )
  return(res)
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
.plot_lda_save <- function(lda_obj, meta, lda_class, filename, title = "LDA") {
  if (!(lda_class %in% colnames(meta))) stop(paste0("lda_class '", lda_class, "' not found in metadata"))
  
  p <- plot_lda_gg(
    lda = lda_obj,
    metadata = meta,
    class_col = lda_class,  # use the actual LDA class passed
    color_by = lda_class,   # color points by LDA class
    shape_by = NULL,        # ignore shape
    label_by = NULL,        # ignore labels
    title = title,
    point_size = 3
  )
  ggplot2::ggsave(filename = filename, plot = p, width = 10, height = 8)
  invisible(p)
}
#' LDA plotting helper (fixed)
plot_lda_gg <- function(lda, metadata, class_col, dim_axes = c(1,2), color_by = NULL,
                        shape_by = NULL, label_by = NULL, label_size = 3,
                        title = "LDA plot", point_size = 3, alpha = 0.8,
                        theme = ggplot2::theme_minimal()) {
  if (!inherits(lda, "lda")) stop("`lda` must be a MASS::lda object")
  if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame")
  if (!("sample_id" %in% colnames(metadata))) stop("`metadata` must contain a 'sample_id' column")

  lda_scores <- as.data.frame(predict(lda)$x)
  colnames(lda_scores) <- paste0("LD", seq_len(ncol(lda_scores)))

  # --- IMPORTANT: assume metadata rows are already in the same order as the lda training/prediction rows
  # Attach sample_id straight from metadata (preserves order), then left_join to keep that order
  lda_scores$sample_id <- metadata$sample_id
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install/attach dplyr for plotting: install.packages('dplyr')")
  }
  lda_scores <- dplyr::left_join(lda_scores, metadata, by = "sample_id")

  # axis labels
  xlab <- paste0("LD", dim_axes[1])
  ylab <- paste0("LD", dim_axes[2])
  xcol <- paste0("LD", dim_axes[1])
  ycol <- paste0("LD", dim_axes[2])

  # Build ggplot mapping: include color only if the column exists in the merged data
  if (!is.null(color_by) && color_by %in% colnames(lda_scores)) {
    p <- ggplot2::ggplot(lda_scores, ggplot2::aes_string(x = xcol, y = ycol, colour = color_by))
  } else {
    p <- ggplot2::ggplot(lda_scores, ggplot2::aes_string(x = xcol, y = ycol))
  }

  p <- p + ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(title = title, x = xlab, y = ylab, colour = color_by) +
    theme

  # optional labels (keeps previous behaviour)
  if (!is.null(label_by) && label_by %in% colnames(lda_scores)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(label = label_by), size = label_size)
    } else {
      p <- p + ggplot2::geom_text(ggplot2::aes_string(label = label_by), size = label_size, vjust = -1)
    }
  }

  return(p)
}
#' PCA plotting helper using ggplot2 (fixed)
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
    # If metadata has sample_id, reorder it to match pca rownames
    if ("sample_id" %in% colnames(metadata)) {
      metadata <- metadata[match(rownames(pca$x), metadata$sample_id), , drop = FALSE]
    }
    if (nrow(metadata) != nrow(scores)) stop("`metadata` rows must match `pca$x` rows (after ordering by sample_id if present)")
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

  # base plot with x/y
  aes_args <- list(x = colnames(scores)[1], y = colnames(scores)[2])
  # add aesthetics as strings (aes_string expects character names)
  if (!is.null(color_by)) aes_args$colour <- color_by
  if (!is.null(shape_by)) aes_args$shape <- shape_by
  if (!is.null(size_by))  aes_args$size <- size_by

  p <- ggplot2::ggplot(scores, do.call(ggplot2::aes_string, aes_args))

  if (!is.null(draw_by)) {
    if (ellipse_fill) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(fill = draw_by, group = draw_by), geom = "polygon", alpha = ellipse_alpha, show.legend = FALSE)
    } else {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(color = draw_by, group = draw_by), size = ellipse_size, show.legend = FALSE)
    }
  }

  p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)

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
