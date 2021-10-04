#' @title Get minor allele frequency vector
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a gwaa-data, a vcfR or a gwasim object
#' @return vector of minor allele frequencies
#' @export
#'
get_maf <- function(x) {
  x <- as.matrix(x)
  if (sum(na.omit(x) > 2) != 0) {
	stop('Currently only the 0,1,2 GP map is supported!')
  }
  n_alleles <- 2 * nrow(x)
  maf <- colSums(x, na.rm = T) / n_alleles # allele count / total number of alleles
  validated <- validate_maf(maf)
  return(validated$maf)
}
