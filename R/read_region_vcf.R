#' @title read specific region from VCF file to a gwasim type genotype matrix
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates (locus), the function extracts genotypes for all loci in that region,
#' constructs genotypes matrix G (applies imputation and fixes allele encoding) and computes maf per marker.
#' @param locus region to extract, e.g. 22:17638600-17641201
#' @param vcf_file path to vcf file containing genotypes
#' @param force_silent silence seqminer verbosity (may cause some issues on non-*nix machines)
#' @return a list with genotypes matrix and maf vector
#' @export
#'
read_region_vcf <- function(locus, vcf_file, force_silent = T, GP_map = c(0, 1, 2)) {
  if (force_silent) {
    region <- hush(seqminer::readVCFToMatrixByRange(vcf_file, range = locus))[[1]]
  } else {
    region <- seqminer::readVCFToMatrixByRange(vcf_file, range = locus)[[1]]
  }

  G <- region %>%
    t() %>%
    gwasim::impute_G() %>%
    gwasim::fix_allele_encoding()

  if (!('gwasim' %in% class(G)))
    class(G) <- c(class(G), 'gwasim')

  maf <- gwasim::get_maf(G)
  if (GP_map != c(0,1,2)) {
    G <- recode_G(G = G, GP_map = GP_map)
  }
  result = list(G = G, maf = maf)
  return(result)
}
