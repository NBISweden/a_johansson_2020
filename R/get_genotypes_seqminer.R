#' @title Get genotypes matrix (counts of reference allele) from one of the seqminer
#' read functions that return genotype matrix
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x an object returned by seqminer *Matrix* function
#' @param ind_names names of the individuals to extract
#' @param marker_names names of the markers to extract
#' @return genotype matrix ind x marker with counts of the reference allele
#'
get_genotypes_seqminer_matrix <- function(x, ind_names = NULL, marker_names = NULL) {
  if (is.null(marker_names)) {
    tmp <- x@gtdata
  } else {
    tmp <- x@gtdata[, marker_names]
  }
  if (!is.null(ind_names)) {
    tmp <- tmp[ind_names, ]
  }

  G <- tmp %>%
    as.numeric() %>%
    as.matrix()
  colnames(G) <- colnames(tmp)
  rownames(G) <- rownames(tmp)
  return(G)
}
