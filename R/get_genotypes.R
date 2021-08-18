#' @title Get genotypes matrix (counts of reference allele) from a GenABEL gwaa-data, vcfR object or
#' an object imported using seqminer
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a gwaa-data, a vcfR object or a genotype matrix loaded using seqminer
#' @param ind_names names of individuals to extract
#' @param marker_names names of markers to extract
#' @return genotype matrix ind x marker with counts of the reference allele
#' @details if marker_names is not provided, all markers present in x will be returned
#' similarly, if ind_names is not present all individuals will be returned
#' NOTE: beware that the genotype matrix returned contains counts of a reference allele
#' which NOT NECESSARILY is the minor allele. Use gwasim::validate_maf() to be sure
#' you are working with the minor alllele counts. Also, note that removing individuals
#' may change allele frequencies and thus different alllele will be considered minor.
#' @export
#'
get_genotypes <- function(x, ind_names = NULL, marker_names = NULL) {
  obj_type <- class(x)
  if (obj_type == "gwaa.data") {
    G <- get_genotypes_gwaa_data(x, ind_names, marker_names)
  } else if (obj_type == "vcfR") {
    G <- get_genotypes_vcfR(x, ind_names, marker_names)
  } else if ("matrix" %in% obj_type & "array" %in% obj_type) {
    message("Assuming input data from seqminer package.")
    G <- get_genotypes_seqminer(x, ind_names, marker_names)
  } else {
    stop('Data format not supported, expecting GenABEL gwaa.data, vcfR object or a genotypes matrix from seqminer')
  }
  if (!('gwasim' %in% class(G))) class(G) <- c(class(G), 'gwasim')
  return(G)
}
