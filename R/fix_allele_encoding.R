#' @title  Fixes allele encoding so that the minor allele homozygote is always encoded as 2
#' @author Marcin Kierczak <marcin.kierczak__REPLACE_WITH_AT__scilifelab.se>
#' @details The genotype matrix G should be returning counts of minor allele.
#' If the major allele is not the reference allele, one needs to fix encoding and this
#' function takes care of this task. It counts reference and alternative allele and
#' if the reference allele is the minor allele, it swaps encoding for homozygotes.
#' @param G - goenotype matrix
#' @return genotype matrix G (row - ind, col - marker) with fixed encoding:
#' 0 - major allele homozygote,
#' 1 - heterozygote,
#' 2 - minor alllele homozygote
#' being simply the count of the minor allele
#'
fix_allele_encoding <- function(G) {
  for (marker in colnames(G)) {
    #marker <- colnames(G)[1]
    ref_cnt <- sum(G[,marker] == 0, na.rm = T) + sum(G[,marker] == 1, na.rm = T)
    alt_cnt <- sum(G[,marker] == 2, na.rm = T) + sum(G[,marker] == 1, na.rm = T)
    if (ref_cnt < alt_cnt) {
      tmp <- G[, marker]
      ref <- which(tmp == 0)
      alt <- which(tmp == 2)
      tmp[ref] <- 2
      tmp[alt] <- 0
      G[,marker] <- tmp
    }
  }
  return(G)
}