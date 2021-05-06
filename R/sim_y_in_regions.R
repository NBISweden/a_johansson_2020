#' @title simulate phenotypes in a number of regions
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates (locus), the function extracts genotypes for all loci in that region,
#' constructs genotypes matrix G (applies imputation and fixes allele encoding) and computes maf per marker.
#' @param locus region to extract, e.g. 22:17638600-17641201
#' @param vcf_file path to vcf file containing genotypes
#' @param force_silent silence seqminer verbosity (may cause some issues on non-*nix machines)
#' @return a list with genotypes matrix and maf vector
#'
sim_y_in_regions <- function(regions, vcf_file, n_markers, get_betas_fun = dbeta,
                             get_betas_args = list(shape1 = .1, shape2 = .1), rare = T,
                             frac_negative=0, thr=0.01, force_silent = T, seed = 42,
                             e_mean = 0, e_sd = 1) {
  registerDoFuture()
  y <- foreach(x = regions) %dorng% {
    tmp <- read_region_vcf(locus = region, vcf_file = vcf_file, force_silent = T)
    effects <- gwasim::get_effects(maf=tmp$maf, thr=thr, N = n_markers, get_betas_fun = get_betas_fun,
                                   get_betas_args = get_betas_args,
                                   rare=rare,
                                   frac_negative=frac_negative,
                                   seed = seed)
    ## Simulate phenotype
    error_term <- rnorm(n = dim(tmp$G)[1], mean = e_mean, sd = e_sd)
    phenos <- pmax(0, tmp$G[, effects$marker_idx] %*% effects$effects + error_term)
    tmp <- list(region = x, y = phenos)
    tmp
  }
  return(y)
}
