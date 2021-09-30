#' @title Count associated positives
#' @description Summarizes results from scan_CommonRare.
#' Counts number of times a non-focal p-value was lower or equal to the p-value from the focal region.
#' @param x a list of scan_CommonRare results
#' @return a vector of counts
#' @author Marcin Kierczak
#' @export
count_ap <- function(x) {
  focal <- x[x$is_focal,]
  non_focal <- x[!x$is_focal,]
  ap <- length(na.omit(which(non_focal$p_value <= focal$p_value)))
  return(ap)
}
