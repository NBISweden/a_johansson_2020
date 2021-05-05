#' @title count number of loci below and above a maf threshold
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates, the function extracts genotypes for all loci in the region,
#' computes minor allele frequency for each marker and returns the number of loci with maf <= threshold and
#' the number of loci maf > threshold.
#' @param regions list of regions as returned by the read_regions_bed function
#' @param vcf_file path to vcf file containing genotypes
#' @param maf_threshold maf threshold for counting loci
#' @return a tibble with region name, total count of loci, count of loci above and below the maf threshold
#' @export
#'
count_by_maf <- function(regions, vcf_file, maf_threshold = 0.01) {
  registerDoFuture()
  y <- foreach(x = regions$region) %dopar% {
    #geno <- seqminer::readVCFToMatrixByRange(vcf_file, range = x)[[1]]
    #n_alleles <- 2 * dim(geno)[2]
    #maf <- pmin(rowSums(geno) / n_alleles, 1 - (rowSums(geno) / n_alleles))
    tmp <- read_region_vcf(locus = x, vcf_file = vcf_file, force_silent = T)
    count_gt <- sum(tmp$maf > maf_threshold)
    count_leq <- sum(tmp$maf <= maf_threshold)
    n_markers <- length(tmp$maf)
    tmp <- list(region = x, n_loci = n_markers, n_gt = count_gt, n_leq = count_leq)
    tmp
  }
  result <- map(y, ~data.frame(.)) %>%
    map_dfr(~mutate_all(., as.character)) %>%
    as_tibble() %>%
    mutate_at(c(2,3,4), as.integer)
  return(result)
}


