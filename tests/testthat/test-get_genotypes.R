test_that("getting genotypes works for gwaa.data", {
  require('GenABEL')
  data(srdta)
  snps <- srdta@gtdata@snpnames[1:5]
  G <- get_genotypes(srdta, snps)
  expect_equal(dim(G), c(2500, 5))
  expect_true(all(unique(as.vector(G)) %in% c(0, 1, 2, NA)))
})

test_that("getting genotypes works for vcfR data", {
  ## TODO: make it proper!
  require('vcfR')
  data(vcfR_test)
  snps <- c("rs6054257", "rs6040355", "microsat1")
  expect_warning(get_genotypes(vcfR_test, snps))

})
