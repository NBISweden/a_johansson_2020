#' @title read specific region from VCF file to a gwasim type genotype matrix
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates (locus), the function extracts genotypes for all loci in that region,
#' constructs genotypes matrix G (applies imputation and fixes allele encoding) and computes maf per marker.
#' @param locus region to extract, e.g. 22:17638600-17641201
#' @param vcf_file path to vcf file containing genotypes
#' @param force_silent silence seqminer verbosity (may cause some issues on non-*nix machines)
#' @param GP_map a genotype-phenotype map, by default AA-0, Aa-1 and aa-2
#' @param missing_thr max proportion of missing genotypes per locus (default 5 per-cent)
#' @return a list with genotypes matrix and maf vector or NULL if region cannot be extracted
#' @details Loci (columns) with too high number of missing values will be removed
#' @export
#'
read_region_vcf <- function(locus, vcf_file, force_silent = F, GP_map = c(0, 1, 2), missing_thr = 0.05) {
  if (force_silent) {
    region <- hush(seqminer::readVCFToMatrixByRange(vcf_file, range = locus))[[1]]
  } else {
    region <- seqminer::readVCFToMatrixByRange(vcf_file, range = locus)[[1]]
  }
  result <- NULL

  if (dim(region)[1] != 0) {
    G <- region %>%
      replace(. < 0 | . > 2, NA) %>%
      t() %>%
      impute_G() %>%
      fix_allele_encoding()

    if (!('gwasim' %in% class(G)))
      class(G) <- c(class(G), 'gwasim')

    if (!all(GP_map == c(0,1,2))) {
      G <- recode_G(G = G, GP_map = GP_map)
    }
    prop_missing <- colSums(is.na(G)) / nrow(G)
    valid_columns <- prop_missing <= missing_thr
    G <- G[, valid_columns]
    maf <- get_maf(G)
    result = list(G = G, maf = maf)
  }
  return(result)
}
