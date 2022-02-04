impute <- function(v, mode='naive') {
  # fix encoding by swaping codes if the alternative allele is not the minor one
  if (sum(v == 2) > sum(v == 0)) {
    warning('Swapping alleles. The ALT was not the minor one!')
    w <- v
    w[v == 0] <- 2
    w[v == 2] <- 0
    v <- w
  }
  
  N <- 2 * sum(v %in% c(0,1,2))
  missing_ids <- which(!(v %in% c(0,1,2)))
  
  if (mode == 'HW') {
    f_a1 <- (2 * sum(v == 0) + sum(v == 1)) / N
    f_a2 <- (2 * sum(v == 2) + sum(v == 1)) / N
    p2 <- f_a1 * f_a1
    pq <- f_a1 * f_a2
    q2 <- f_a2 * f_a2
    probs <- c(p2,pq,q2)
    probs[is.na(probs)] <- 0
    if (all(probs) == 0) { probs = c(.33, .33, .33) }
    v[missing_idx] <- sample(c(0,1,2),
                             size = length(missing_idx),
                             replace = T,
                             prob = probs)
  } else if (mode == 'naive') {
    v[missing_ids] <- 2
  } else {
    warning('Imputation mode not supported. Assuming naive!')
    v[missing_ids] <- 2
  }
  return(v)
}
