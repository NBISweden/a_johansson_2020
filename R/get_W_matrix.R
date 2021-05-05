#' @title get weights matrix based on genotype-phenotype map
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  given a 3 column genotype-phenotype map matrix with row per marker and columns
#' representing weights for genotypes: 'aa', 'aA/Aa' and 'AA' encoded
#' as 0, 1 and 2 in the genotype matrix, compute weights matrix. If no GP_map is provided,
#' a matrix of ones will be returned meaning no weighting.
#' @param G genotypes matrix
#' @param GP_map genotype-phenotype map
#' @return weights matrix
#'
get_W_matrix <- function(G, GP_map = NULL) {
  W <- G
  if (is.null(GP_map)) {
    W[!is.null(W)] <- 1
  } else {
    for (i in 1:dim(G)[2]) {
      x <- G[ ,i]
      W[x == 0, i] <- GP_map[i, 1]
      W[x == 1, i] <- GP_map[i, 2]
      W[x == 2, i] <- GP_map[i, 3]
    }
  }
return(W)
}

# test
# 4 individuals, 5 markers
# betas <- rnorm(n = 5, 0, 1)
# G <- matrix(sample(c(0,1,2), size = 20, replace = T), nrow=4)
# GP_map <- matrix(sample(c(0,0.5,1), size = 15, replace = T), ncol=3)
# get_W_matrix(G, GP_map)
