
## File: R/matrix_prep.R
#' Prepare methylation percentage matrix from a methylKit object
#' Wrapper that calls small subfunctions to remove duplicates, filter samples, and impute
#' @param meth_obj methylKit object (e.g. methylBase, methylRawList)
#' @param feature_name name used in console messages
#' @importFrom impute impute.knn
#' @export
#' Prepare methylation percentage matrix from a methylKit object
#' Wrapper that calls small subfunctions to remove duplicates, filter samples, and impute
#' @param meth_obj methylKit object (e.g. methylBase, methylRawList)
#' @param feature_name name used in console messages
#' @param remove_samples character vector of sample IDs to remove
#' @importFrom impute impute.knn
#' @export
prepare_methylation_matrix <- function(meth_obj, feature_name = "feature", remove_samples = NULL) {
  df <- .coerce_to_df(meth_obj)
  df <- .remove_duplicate_loci(df)
  
  PercM <- methylKit::percMethylation(meth_obj, rowids = TRUE)
  cat(feature_name, "original matrix dimensions:", dim(PercM), "\n")
  
  # Remove samples if requested
  if (!is.null(remove_samples)) {
    samples_exist <- remove_samples %in% colnames(PercM)
    if (any(samples_exist)) {
      PercM <- PercM[, !colnames(PercM) %in% remove_samples, drop = FALSE]
      cat("Removed samples:", paste(remove_samples[samples_exist], collapse = ", "), "\n")
    }
    if (!all(samples_exist)) {
      cat("Warning: these samples were not found and could not be removed:", 
          paste(remove_samples[!samples_exist], collapse = ", "), "\n")
    }
  }
  
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
  total_vals <- length(mat)
  missing_before <- sum(is.na(mat))
  res <- impute::impute.knn(as.matrix(mat), colmax = colmax)$data
  missing_after <- sum(is.na(res))
  imputed_count <- missing_before - missing_after
  imputed_pct <- round(100 * imputed_count / total_vals, 2)

  cat("Missing values imputed:", imputed_count, 
      "(", imputed_pct, "%)", "\n",
      "Remaining missing values:", missing_after, "\n",
      "Final matrix dimensions:", dim(res), "\n\n")
  res
}

#' Check consistency between methylation matrix and metadata
#'
#' Verifies that all sample IDs in metadata are present in the methylation matrix,
#' all matrix columns are in the metadata, and that the order matches.
#'
#' @param mat A numeric matrix or data.frame with samples as columns
#' @param metadata A data.frame containing a column of sample IDs
#' @param sample_col Name of the column in metadata containing sample IDs (default: "sample_id")
#' @export
check_sample_consistency <- function(mat, metadata, sample_col = "sample_id") {
  if (!sample_col %in% colnames(metadata)) {
    stop("Error: '", sample_col, "' column not found in metadata!")
  }
  
  sample_ids <- metadata[[sample_col]]
  
  if (!all(sample_ids %in% colnames(mat))) {
    stop("Error: Some sample IDs in metadata are missing from the matrix!")
  }
  if (!all(colnames(mat) %in% sample_ids)) {
    stop("Error: Some sample IDs in the matrix are missing from metadata!")
  }
  if (!all(sample_ids == colnames(mat))) {
    stop("Error: Sample IDs in metadata and matrix are not in the same order!")
  }
  
  cat("Matrix and metadata contain the same sample IDs and are in the same order.\n")
  invisible(TRUE)
}
