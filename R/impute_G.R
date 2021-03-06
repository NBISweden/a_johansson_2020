#' @title Impute missing genotypes using minor allele frequency
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param G - genotype matrix (individuals in columns, markers in rows) to impute
#' @description Imputes missing genotypes using minor allele frequency and assuming Hardy-Weinberg equilibrium
#' @return an imputed matrix
#' @export

impute_G <- function(G) {
  G <- fix_allele_encoding(G)
  # Impute
  for (marker in colnames(G)) {
    marker_maf <- (2 * sum(G[,marker] == 2, na.rm = T) + sum(G[,marker] == 1, na.rm = T)) / sum(!is.na(G[,marker]))
    # If all alleles are missing, assume monomorphic for the major allele
    if (!is.finite(marker_maf) | is.nan(marker_maf)) {
      warning(paste0("\nToo many missing genotypes in ", marker, ". Assuming monomorphic for the major allele. You may want to investigate this further."))
      marker_maf <- 0
    }
    p_0 <- (1 - marker_maf) ^ 2
    p_1 <- marker_maf * (1 - marker_maf)
    p_2 <- marker_maf^2
    probs <- c(p_0, p_1, p_2)
    for (ind in 1:dim(G)[1]) {
      if (is.na(G[ind, marker])) {
        G[ind, marker] <- sample(c(0, 1, 2), size = 1, prob = probs)
      }
    }
  }
  return(G)
}
