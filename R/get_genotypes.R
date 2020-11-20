#' @title Get genotypes matrix (counts of reference allele) from a GenABEL gwaa-data or vcfR object
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a gwaa-data or a vcfR object
#' @param marker_names names of markers to extract
#' @return genotype matrix ind x marker with counts of the reference allele
#' @export
#'
get_genotypes <- function(x, marker_names) {
  obj_type <- class(x)
  if (obj_type == "gwaa.data") {
    G <- get_genotypes_gwaa_data(x, marker_names)
  } else if (obj_type == "vcfR") {
    G <- get_genotypes_vcfR(x, marker_names)
  } else {
    stop('Data format not supported, expecting GenABEL gwaa.data or vcfR')
  }
  return(G)
}
