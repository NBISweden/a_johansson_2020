#' @title read specific region from VCF file to a gwasim type genotype matrix
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates (locus), the function extracts genotypes for all loci in that region,
#' constructs genotypes matrix G (applies imputation and fixes allele encoding) and computes maf per marker.
#' @param locus region to extract, e.g. 22:17638600-17641201
#' @param vcf_file path to vcf file containing genotypes
#' @param force_silent silence seqminer verbosity (may cause some issues on non-*nix machines)
#' @param GP_map a genotype-phenotype map, by default AA-0, Aa-1 and aa-2
#' @return a list with genotypes matrix and maf vector or NULL if region cannot be extracted
#' @export
#'
read_region_vcf <- function(locus, vcf_file, force_silent = F, GP_map = c(0, 1, 2)) {
  if (force_silent) {
    region <- hush(seqminer::readVCFToMatrixByRange(vcf_file, range = locus))[[1]]
  } else {
    region <- seqminer::readVCFToMatrixByRange(vcf_file, range = locus)[[1]]
  }
  result <- NULL

  if (dim(region)[1] != 0) {
    G <- region %>%
      t() %>%
      impute_G() %>%
      fix_allele_encoding()

    if (!('gwasim' %in% class(G)))
      class(G) <- c(class(G), 'gwasim')

    maf <- get_maf(G)
    if (!all(GP_map == c(0,1,2))) {
      G <- recode_G(G = G, GP_map = GP_map)
    }
    result = list(G = G, maf = maf)
  }
  return(result)
}
