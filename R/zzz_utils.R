## File: R/zzz_utils.R
# Utility helpers and simple internal functions

#' Null-coalescing operator
#'
#' Return `a` if not NULL, otherwise `b`.
#' @param a Any
#' @param b Any
#' @return `a` or `b`
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Ensure a directory exists
#' @noRd
.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Save PNG wrapper
#' @noRd
.save_png <- function(expr, filename, width = 800, height = 600) {
  .ensure_dir(dirname(filename))
  png(filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

.onLoad <- function(libname, pkgname) {
    required_pkgs <- c(
        "methylKit",
        "SummarizedExperiment",
        "GenomicRanges",
        "IRanges",
        "S4Vectors",
        "BiocGenerics"
    )
    
    for (p in required_pkgs) {
        if (!requireNamespace(p, quietly = TRUE)) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
            }
            BiocManager::install(p, ask = FALSE, update = FALSE)
        }
    }
}
