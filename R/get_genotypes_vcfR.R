#' @title Get genotypes matrix (counts of reference allele) from a vcfR object
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x a vcfR object
#' @param marker_names names of markers to extract
#' @return genotype matrix ind x marker with counts of the reference allele
get_genotypes_vcfR <- function(x, marker_names) {
  tmp <- extract.gt(x, return.alleles = F, as.numeric = F)[marker_names,]

  G <- tmp %>%
    str_replace(., pattern = '0\\|0', replacement = '0') %>%
    str_replace(., pattern = '0\\|1', replacement = '1') %>%
    str_replace(., pattern = '1\\|0', replacement = '1') %>%
    str_replace(., pattern = '1\\|1', replacement = '2') %>%
    as.numeric() %>%
    matrix(ncol = length(marker_names), byrow = T)
  colnames(G) <- rownames(tmp)
  rownames(G) <- colnames(tmp)
  return(G)
}
