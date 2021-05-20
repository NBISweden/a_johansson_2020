#' @title Effect sizes $\Beta$ from a normal distribution.
#' @description Given a vector of minor-allele frequencies (maf), mean and standard deviation, sample `length(maf)`
#' effects from a normal distribution.
#' @author Marcin Kierczak <marcin.kierczak_NO_SPAM_scilifelab.se>
#' @param maf  vector of minor-allele frequencies.
#' @param mean  vector of means.
#' @param sd  vector of standard deviations.
#' @return vector of effect sizes
#' @export
#'
rnorm_betas <- function(maf, mean = 0, sd = 1) {
  return(rnorm(n = length(maf), mean = mean, sd = sd))
}
