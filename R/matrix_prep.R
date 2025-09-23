
## File: R/matrix_prep.R
#' Prepare methylation percentage matrix from a methylKit object
#' Wrapper that calls small subfunctions to remove duplicates, filter samples, and impute
#' @param meth_obj methylKit object (e.g. methylBase, methylRawList)
#' @param feature_name name used in console messages
#' @importFrom impute impute.knn
#' @export
prepare_methylation_matrix <- function(meth_obj, feature_name = "feature") {
  df <- .coerce_to_df(meth_obj)
  df <- .remove_duplicate_loci(df)
  PercM <- methylKit::percMethylation(meth_obj, rowids = TRUE)
  cat(feature_name, "original matrix dimensions:", dim(PercM), "
")
  PercM <- .filter_sample_missingness(PercM, max_missing = 0.8)
  PercM_imputed <- .impute_knn(PercM, colmax = 0.8)
  return(PercM_imputed)
}

#' @noRd
.coerce_to_df <- function(meth_obj) {
  if (is.data.frame(meth_obj)) return(meth_obj)
  if (inherits(meth_obj, "methylBase") || inherits(meth_obj, "methylRawList")) {
    return(as.data.frame(meth_obj))
  }
  stop("Unsupported meth_obj class; please provide methylKit object or data.frame")
}
#' @noRd
.remove_duplicate_loci <- function(df) {
  if (all(c("chr", "start") %in% names(df))) {
    key <- paste0(df$chr, ".", df$start, ".", ifelse("end" %in% names(df), df$end, df$start))
    no_dupes <- !duplicated(key)
    df <- df[no_dupes, , drop = FALSE]
  }
  df
}
#' @noRd
.filter_sample_missingness <- function(mat, max_missing = 0.8) {
  missing_frac <- colMeans(is.na(mat))
  cols_keep <- missing_frac <= max_missing
  removed_samples <- sum(!cols_keep)
  if (removed_samples > 0) {
    cat("Removing", removed_samples, "samples with >", max_missing * 100, "% missingness:
",
        paste(colnames(mat)[!cols_keep], collapse = ", "), "
")
  }
  mat[, cols_keep, drop = FALSE]
}
#' @noRd
.impute_knn <- function(mat, colmax = 0.8) {
  missing_before <- sum(is.na(mat))
  res <- impute::impute.knn(as.matrix(mat), colmax = colmax)$data
  missing_after <- sum(is.na(res))
  cat("Missing values imputed:", missing_before - missing_after, "
",
      "Remaining missing values:", missing_after, "
",
      "Final matrix dimensions:", dim(res), "

")
  res
}
