#' @title Get minor allele frequency vector
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a gwaa-data, a vcfR or a gwasim object
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
  } else if ("gwasim" %in% obj_type) {
    maf <- colSums(x) / (2 * dim(x)[1])
  } else {
    stop('ERROR! Data format not supported. Expecting GenABEL gwaa.data or vcfR.')
  }
  validated <- validate_maf(maf)
  return(validated$maf)
}
