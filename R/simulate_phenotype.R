#' @title Simulate phenotype given a set of common and rare variants.
#' @author Marcin Kierczak <marcin.kierczak__INSERT_AT__scilifelab.se>
#' @param  x a vcf object from vcfR package
#' @param thr_common_rare maf threshold for common and rare alleles
#' @param N_rare number of rare alleles to affect the trait
#' @param beta_params_rare parameters of the Beta distribution to simulate effect of rare alleles
#' @param perc_negative_rare how many (per cent) rare alleles have negative effect
#' @param N_common number of common alleles to affect the trait
#' @param beta_params_rare parameters of the Beta distribution to simulate effect of common alleles
#' @param perc_negative_common how many (per cent) common alleles have negative effect
#' @param e a vector of mean and sd for the error term. Default (0, 1).
#' @value vector of simulated phenotypes

simulate_phenotype <- function(x,
                               thr_common_rare = 0.01,
                               N_rare,
                               beta_params_rare,
                               perc_negative_rare,
                               N_common,
                               beta_params_common,
                               perc_negative_common,
                               e = c(0, 1)) {

  maf <- vcfR::maf(x)[,'Frequency']
  if (n_rare > 0) {
    rare <- get_effects(maf = maf, thr = thr_common_rare,
                      N = n_rare,
                      shape12 = beta_params_rare,
                      below = T,
                      perc_negative = perc_negative_rare)
    G_rare <- get_genotypes(x = vcf, marker_names = names(rare$marker_idx)) %>%
      fix_allele_encoding() %>%
      impute_G(maf = maf)
  } else {
    rare <- list(effects = 0)
  }
  if (n_common > 0) {
    common <- get_effects(maf = maf, thr = thr_common_rare,
                          N = N_common,
                          shape12 = beta_params_common,
                          below = F,
                          perc_negative = perc_negative_common)
    G_common <- get_genotypes(x = vcf, marker_names = names(common$marker_idx)) %>%
      fix_allele_encoding() %>%
      impute_G(maf = maf)

  } else {
    common <- list(effects = 0)
  }

  y_rare <- G_rare %*% rare$effects
  y_common <- G_common %*% common$effects
  e <- rnorm(n = dim(G_rare)[1], mean = e[1], sd = e[2])
  y <- y_common + y_rare + e
  return(y)
}
