#' @title simulate phenotypes in a number of regions
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given region coordinates (locus), the function extracts genotypes for all loci in that region,
#' constructs genotypes matrix G (applies imputation and fixes allele encoding) and computes maf per marker.
#' @param regions list of regions to use
#' @param vcf_file path to vcf file containing genotypes
#' @param n_markers number of markers to use in the simulation (see get_effects)
#' @param get_betas_fun see get_effects
#' @param get_betas_args see get_effects
#' @param rare see get_effects
#' @param frac_negative see get_effects
#' @param thr see get_effects
#' @param seed see get_effects
#' @param e_mean mean error
#' @param e_sd error std. dev.
#' @param GP_map genotype-phenotype map
#' @param force_silent silence seqminer verbosity (may cause some issues on non-*nix machines)
#' @return a list with genotypes matrix and maf vector
#' @export
#'
sim_y_in_regions <- function(regions, vcf_file, n_markers, get_betas_fun = dbeta,
                             get_betas_args = list(maf=T, shape1 = .1, shape2 = .1), rare = T,
                             frac_negative=0, thr=0.01, force_silent = T, seed = 42,
                             e_mean = 0, e_sd = 1, GP_map = c(0,1,2)) {

  y <- foreach(x = regions) %dorng% {
    tmp <- read_region_vcf(locus = x, vcf_file = vcf_file, force_silent = T, GP_map = GP_map)
    if (!is.null(tmp)) {
      effects <- get_effects(maf = tmp$maf, thr = thr, N = n_markers, get_betas_fun = get_betas_fun,
                                   get_betas_args = get_betas_args,
                                   rare = rare,
                                   frac_negative = frac_negative,
                                   seed = seed)
      ## Simulate phenotype
      error_term <- rnorm(n = dim(tmp$G)[1], mean = e_mean, sd = e_sd)
      phenos <- pmax(0, tmp$G[,effects$marker_idx] %*% effects$effects + error_term)
      tmp <- list(region = x, y = phenos)
      return(tmp)
    } else {
      return(NULL)
    }
  }
  return(y)
}
