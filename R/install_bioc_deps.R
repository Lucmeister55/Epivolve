.onAttach <- function(libname, pkgname) {
  # Bioconductor dependencies
  bioc_pkgs <- c(
    "methylKit",
    "bsseq",
    "GenomicRanges",
    "GenomicFeatures",
    "sva",
    "limma"
  )
  
  missing_pkgs <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
  
  if(length(missing_pkgs) > 0) {
    packageStartupMessage("Installing missing Bioconductor dependencies: ", paste(missing_pkgs, collapse = ", "))
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  }
}
