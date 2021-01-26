#' @title Simulate phenotype for a single locus based on 3 sets of dieffect distribution parameters.
#' @author Marcin Kierczak <marcin.kierczak__INSERT_AT__scilifelab.se>
#' @param  x a vcf object from vcfR package
#' @param thr_common_rare maf threshold for common and rare alleles
#' @param N_rare number of rare alleles to affect the trait
#' @param beta_params_rare parameters of the Beta distribution to simulate effect of rare alleles
#' @param perc_negative_rare how many (per cent) rare alleles have negative effect
#' @param N_common number of common alleles to affect the trait
#' @param beta_params_common parameters of the Beta distribution to simulate effect of common alleles
#' @param perc_negative_common how many (per cent) common alleles have negative effect
#' @param e a vector of mean and sd for the error term. Default (0, 1).
#' @return a vector of simulated phenotypes
#' @export
#'
simulate_phenotype <- function(x, locus, means = c(1, 2, 3), std_devs = c(1,1,1)) {
    G <- get_genotypes(x = vcf, marker_names = locus) %>%
      fix_allele_encoding() %>%
      impute_G(maf = maf)

  y_rare <- G_rare %*% rare$effects
  y_common <- G_common %*% common$effects
  e <- rnorm(n = dim(G_rare)[1], mean = e[1], sd = e[2])
  y <- y_common + y_rare + e
  return(y)
}
