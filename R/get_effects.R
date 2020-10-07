#' @title Get weights for markers using beta distribution.
#' @description Given a vector of minor-allele frequencies (maf) for a number of markers and a number of markers to assign an effect to,
#' the function randomly selects N markers below or above provided maf threshold (default 1%) and assigns effect of a size sampled
#' from Beta distribution with given parameters so that the magnitude of the effect depends on the allele frequency. Desired fraction
#' of the effects will be negative. Monomorphic markers with `maf = 0` will be automatically excluded from sampling!
#' @author Marcin Kierczak <marcin.kierczak_NO_SPAM_scilifelab.se>
#' @param maf - vector of minor-allele frequencies
#' @param N - number of markers to be assigned an effect
#' @param shape12 - a two element vector of Beta distribution shapes.
#' @param below - a boolean, if F markers below and equal to the `thr` will be sampled
#' @param perc_negative - percentage of effects to be set to negative
#' @return a list with two elements: `marker_idx` is a vector of indices of the markers assigned an effect and `effects` is a vector of
#' effects, one for each marker
get_effects <- function(maf, N, shape12, thr=0.01, below=T, perc_negative=.2) {
  if (!below) {
    idx <- sample(which(maf > thr), N, replace = F)
  } else {
    idx <- sample(which(maf <= thr & maf > 0), N, replace = F)
  }
  signs <- sample(c(-1,1), N, replace = T, prob = c(perc_negative, 1 - perc_negative_common))
  betas <- dbeta(x = maf[idx], shape1 = shape12[1], shape2 = shape12[2]) * signs
  output <- list(marker_idx = idx, effects = betas)
  return(output)
}
