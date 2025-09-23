#' Association between PCs and sample annotations (assocComp-like)
#' @param data_df numeric matrix/data.frame where rows = samples, cols = features
#' @param sampleAnnotation data.frame where rows = samples (same order) and columns = annotations
#' @export
assocComp_df <- function(data_df, sampleAnnotation) {
  if (nrow(data_df) != nrow(sampleAnnotation)) stop("Number of samples in data and sampleAnnotation must match")
  pca_res <- stats::prcomp(data_df, center = TRUE, scale. = TRUE)
  vars <- (pca_res$sdev)^2
  vars <- vars / sum(vars)
  pcs <- pca_res$x
  test_assoc <- function(pc, annotation) {
    if (is.numeric(annotation)) {
      return(stats::cor.test(pc, annotation)$p.value)
    } else if (is.factor(annotation) || is.character(annotation)) {
      annotation <- as.factor(annotation)
      df <- data.frame(pc = pc, group = annotation)
      aov_res <- stats::aov(pc ~ group, data = df)
      pval <- summary(aov_res)[[1]][["Pr(>F)"]][1]
      return(pval)
    } else {
      return(NA)
    }
  }
  association <- sapply(sampleAnnotation, function(annotation_col) {
    sapply(1:ncol(pcs), function(i) test_assoc(pcs[, i], annotation_col))
  })
  rownames(association) <- colnames(pcs)
  colnames(association) <- colnames(sampleAnnotation)
  return(list(pcs = pcs, vars = vars, association = association))
}