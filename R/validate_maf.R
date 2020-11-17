#' @title Validate minor allele frequency vector
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param x maf vector
#' @return a list with problematic elements, the reason why they
#' are problematic and fixed maf values (if possible). If no problems,
#' reason and problems will be NULL.
#'
validate_maf <- function(x) {
  problems <- list(reason = NULL, problems = NULL, maf = maf)
  if (sum(x > 0.5) != 0) {
    reason <- 'Some minor allele frequencies are > 0.5! Converting to 1 - maf!'
    warning(reason)
    maf <- pmin(x, (1 - x))
    problems <- list(reason = reason, problems = x[x > 0.5], maf=maf)
  }
  return(problems)
}
