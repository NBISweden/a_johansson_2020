fix_encoding <- function(v) {
  # fix encoding by swaping codes if the alternative allele is not the minor one
  if (sum(v == 2) > sum(v == 0)) {
    warning('Swapping alleles. The ALT was not the minor one!')
    w <- v
    w[v == 0] <- 2
    w[v == 2] <- 0
    v <- w
  }
  return(v)
}
