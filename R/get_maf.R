#' @title Get minor allele frequency vector
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a gwaa-data or a vcfR object
#' @return vector of minor allele frequencies
#' @export
#'
get_maf <- function(x) {
  n_alleles <- 2 * nrow(x)
  maf <- colSums(x) / n_alleles # allele count / total number of alleles
  validated <- validate_maf(maf)
  return(validated$maf)
}
