#' @title Get genotypes matrix (counts of reference allele) from a GenABEL gwaa-data object
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a vcfR object
#' @param marker_names names of markers to extract
#' @return genotype matrix ind x marker with counts of the reference allele
#'
get_genotypes_gwaa_data <- function(x, marker_names) {
  tmp <- x@gtdata[, marker_names]
  G <- tmp %>%
    as.numeric() %>%
    as.matrix()
  colnames(G) <- colnames(tmp)
  rownames(G) <- rownames(tmp)
  return(G)
}
