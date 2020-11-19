#' @title Get weights for markers using beta distribution.
#' @description Given a vector of minor-allele frequencies (maf) for a number of markers and a number of markers to assign an effect to,
#' the function randomly selects N markers below or above provided maf threshold (default 1%) and assigns effect of a size sampled
#' from Beta distribution with given parameters so that the magnitude of the effect depends on the allele frequency. Desired fraction
#' of the effects will be negative. Monomorphic markers with `maf = 0` will be automatically excluded from sampling!
#' @author Marcin Kierczak <marcin.kierczak_NO_SPAM_scilifelab.se>
#' @param maf - vector of minor-allele frequencies
#' @param N - number of markers to be assigned an effect
#' @param thr - threshold for maf. All markers with maf <= thr will be treated as rare.
#' @param shape12 - a two element vector of Beta distribution shapes.
#' @param rare - a boolean, if TRUE, markers below and equal to the `thr` will be sampled, i.e. the rare variants
#' @param frac_negative - fraction of effects to be set to negative
#' @param seed - set seed for sampling, if FALSE, default `sample` seed will be used
#' @return a list with two elements: `marker_idx` is a vector of indices of the markers assigned an effect and `effects` is a vector of
#' effects, one for each marker
get_effects <- function(maf, N, shape12, thr=0.01, rare=T, frac_negative=0, seed = F) {
  if (seed) {
    set.seed(seed)
  }
  if (max(maf) > 1 | min(maf) < 0) {
    warning(paste0("Weird values of maf detected. min: ", min(maf), " max: ", max(maf), "!"))
  }
  if (rare) {
    valid_markers <- which(maf <= thr & maf > 0)
  } else {
    valid_markers <- which(maf > thr & maf < 1)
  }

  l <- length(valid_markers)
  if (l == 0) {
    warning("No markers matching criteria in the region! Returning NULL!")
    return(NULL)
  } else if ( l < N) {
    warning(paste0('Expected ', N, ' markers while only ', l, ' match maf criteria!'))
  }

  idx <- sample(valid_markers, pmin(l, N), replace = F)
  signs <- sample(c(-1,1), N, replace = T, prob = c(frac_negative, 1 - frac_negative))
  betas <- dbeta(x = maf[idx], shape1 = shape12[1], shape2 = shape12[2]) * signs
  output <- list(marker_idx = idx, effects = betas)
  return(output)
}
