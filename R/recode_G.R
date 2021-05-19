#' @title Re-code genotypes according to the supplied genotype-phenotype map
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  Given a genotype-phenotype (GP) map recodes genotypes in G.
#' The GP map is a matrix with each row containing new coding for one marker in G:
#' new codes for genotypes 'AA', 'aA' and 'aa' (in G encoded as 0, 1 and 2, respectively).
#' @details If only a vector GP map is provided, same re-coding is assumed for each and every
#' locus in G.
#' @param G genotypes matrix
#' @param GP_map genotype-phenotype map
#' @return weights matrix
#' @export
#'
recode_G <- function(G, GP_map) {
  if (is.vector(GP_map) & length(GP_map) == 3) {
    message('One-dimensional GP map provided, assuming same re-coding for all loci in G.')
    GP_map <- matrix(rep(GP_map, times = ncol(G)), ncol = 3, byrow = T)
  } else if (nrow(GP_map) != ncol(G) | ncol(GP_map) != 3) {
    stop(paste0("GP_map has wrong dimensions: ", dim(GP_map), " that do not match dimensions of G: ", dim(G), "!"))
  }
  G2 <- G
  for(i in 1:ncol(G)) {
    gt <- G[ ,i]
    G2[gt == 0, i] <- GP_map[i, 1]
    G2[gt == 1, i] <- GP_map[i, 2]
    G2[gt == 2, i] <- GP_map[i, 3]
  }
  return(G2)
}

