#' @title count number of loci with minor allele being rare
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates, the function extracts genotypes for all loci in the region,
#' computes minor allele frequency for each marker and returns the number of loci with maf <= threshold
#' @param regions list of regions in chr:start-end format, e.g. 22:16191600-16193600
#' @param vcf_file path to vcf file containing genotypes
#' @param maf_threshold maf threshold for counting loci with rare minor allele
#' @return a tibble with region name, count of loci and count or rare variants
#'
count_rare <- function(regions, vcf_file, maf_threshold = 0.01) {
  registerDoFuture()
  y <- foreach(x = regions) %dopar% {
    geno <- seqminer::readVCFToMatrixByRange(vcf_file, range = x)[[1]]
    n_alleles <- 2 * dim(geno)[2]
    maf <- pmin(rowSums(geno) / n_alleles, 1 - (rowSums(geno) / n_alleles))
    rare_count <- sum(maf <= maf_threshold)
    n_markers <- dim(geno)[1]
    tmp <- list(region = x, n_loci = n_markers, n_count = rare_count)
    tmp
  }
  result <- map(y, ~data.frame(.)) %>%
    map_dfr(~mutate_all(., as.character)) %>%
    as_tibble() %>%
    mutate_at(c(2,3), as.integer)

  return(result)
}


