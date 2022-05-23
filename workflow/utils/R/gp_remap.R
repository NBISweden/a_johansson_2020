# Remap genotype-phenotype map 
gp_remap <- function(G, map = c(0, 1, 2)) {
  # If needed, remap G to get a desired genotype-phenotype map, e.g.:
  # G[G == 1] <- 2
  # to say othat ne allele is enough (dominance model)
  G2 <- G
  G2[G == 0] <- map[1]
  G2[G == 1] <- map[2]
  G2[G == 2] <- map[3]

  return(G2) 
}
