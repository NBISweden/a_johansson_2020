#' @title count number of loci with minor allele being rare
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates, the function extracts genotypes for all loci in the region,
#' computes minor allele frequency for each marker and returns the number of loci with maf <= threshold
#' @param region region in chr:start-end format, e.g. 22:16191600-16193600
#' @param vcf_file path to vcf file containing genotypes
#' @param maf_threshold maf threshold for counting loci with rare minor allele
#' @returnthe count of loci with rare variants in the region
#'
count_rare <- function(region, vcf_file, maf_threshold = 0.01) {
  geno <- seqminer::readVCFToMatrixByRange(vcf_file, range = region)[[1]]
  maf <- pmin(rowSums(geno)/dim(geno)[2], 1 - rowSums(geno)/dim(geno)[2])
  rare_count <- sum(maf <= maf_threshold)
  return(rare_count)
}
