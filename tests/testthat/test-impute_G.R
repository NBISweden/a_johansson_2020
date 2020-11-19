test_that("impution of missing genotypes works", {
  G_test <- as.matrix(
    data.frame(
      snp1 = sample(c(0,1,2,NA), size = 100, replace = T, prob = c(.15,.1,.7,.05)),
      snp2 = sample(c(0,1,2,NA), size = 100, replace = T, prob = c(.7,.1,.15,.05)),
      snp3 = rep(0, times = 100),
      snp4 = rep(1, times = 100),
      snp5 = rep(2, times = 100)
    ))

  expect_equal(sum(is.na(impute_G(G_test))), 0)
})
