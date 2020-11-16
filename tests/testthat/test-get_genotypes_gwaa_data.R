test_that("getting genotypes from gwaa-data works", {
  require('GenABEL')
  data(srdta)
  snp_names <- c("rs10","rs18","rs29","rs65","rs73")
  G <- get_genotypes_gwaa_data(srdta, marker_names = snp_names)
  expect_equal(dim(G), c(2500, 5))
  expect_equal(typeof(G), "double")
})
