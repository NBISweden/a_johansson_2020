#' @title Impute missing genotypes using minor allele frequencies
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param G - genotype matrix (individuals in columns, markers in rows) to impute
#' @param maf - a vector with minor allele frequencies for the population from which G comes
#' @return an imputed matrix
impute_G <- function(G, maf) {
  # Impute
  for (marker in colnames(G)) {
    marker_maf <- maf[marker]
    p_0 <- (1 - marker_maf)^2
    p_1 <- marker_maf * (1 - marker_maf)
    p_2 <- marker_maf^2
    probs <- c(p_0, p_1, p_2)
    for (ind in rownames(G)) {
      if (is.na(G[ind, marker])) {
        G[ind, marker] <- sample(c(0, 1, 2), size = 1, prob = probs)
      }
    }
  }
  return(G)
}
