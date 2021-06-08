#' @title Get weights for markers using provided distribution function.
#' @description Given a vector of minor-allele frequencies (maf) for a number of markers and a number of markers to assign an effect to,
#' the function randomly selects N markers below or above provided maf threshold (default 1%) and assigns effect of a size sampled
#' from Beta distribution with given parameters so that the magnitude of the effect depends on the allele frequency. Desired fraction
#' of the effects will be negative. Monomorphic markers with `maf = 0` will be automatically excluded from sampling!
#' @author Marcin Kierczak \email{marcin.kierczak@@scilifelab.se}
#' @param maf - vector of minor-allele frequencies
#' @param N - number of markers to be assigned an effect
#' @param thr - threshold for maf. All markers with maf <= thr will be treated as rare.
#' @param get_betas_fun - name of a function for getting Beta parameters given a vector x of MAFs. Note, if you want a maf-independent
#' function like \link{rnorm}, see \details.
#' @param get_betas_args - a list with additional function parameters for the function defined in `get_betas`
#' @param rare - a boolean, if TRUE, markers below and equal to the `thr` will be sampled, i.e. the rare variants
#' @param frac_negative - fraction of effects to be set to negative
#' @param seed - set seed for sampling, if FALSE, default `sample` seed will be used
#' @return a list with two elements: `marker_idx` is a vector of indices of the markers assigned an effect and `effects` is a vector of
#' effects, one for each marker
#' @export
#'
get_effects <- function(maf, N, get_betas_fun = dbeta, get_betas_args = list(shape1 = 1, shape2 = 25), thr=0.01, rare=T, frac_negative=0, seed = F) {
  if (seed) {
    set.seed(seed)
  }
  if (max(maf) > 1 | min(maf) < 0) {
    warning(paste0("Weird values of maf detected. min: ",
                   min(maf), " max: ", max(maf), "!"))
  }
  if (rare) {
    valid_markers <- which(maf <= thr & maf > 0)
  }
  else {
    valid_markers <- which(maf > thr & maf < 1)
  }
  l <- length(valid_markers)
  if (l == 0) {
    warning("No markers matching criteria in the region! Returning NULL!")
    return(NULL)
  }
  else if (l < N) {
    warning(paste0("Expected ", N, " markers while only ",
                   l, " match maf criteria!"))
    N <- l
  }
  idx <- sample(valid_markers, N, replace = F)
  signs <- sample(c(-1, 1), N, replace = T, prob = c(frac_negative,
                                                     1 - frac_negative))
  sel_maf <- maf[idx]
  betas <- do.call(what = get_betas_fun, args = c(as.name("sel_maf"),
                                                  get_betas_args))
  betas <- betas * signs
  output <- list(marker_idx = idx, effects = betas)
  return(output)
}
