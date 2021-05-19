#' @title Get minor allele frequency vector
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a gwaa-data or a vcfR object
#' @return vector of minor allele frequencies
#' @export
#'
get_maf <- function(x) {
  obj_type <- class(x)
  if ("gwaa.data" %in% obj_type) {
    tmp <- summary(x@gtdata)
    maf <- ((2 * tmp$P.22) + tmp$P.12) / (2 * tmp$NoMeasured)
    names(maf) <- colnames(srdta@gtdata)
  } else if ("vcfR" %in% obj_type) {
    maf <- vcfR::maf(x)[,'Frequency']
  } else {
    #stop('ERROR! Data format not supported. Expecting GenABEL gwaa.data or vcfR.')
    n_alleles <- 2 * nrows(x)
    maf <- colSums(x) / n_alleles # allele count / total number of alleles
  }
  validated <- validate_maf(maf)
  return(validated$maf)
}
