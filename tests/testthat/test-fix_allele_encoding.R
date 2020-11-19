test_that("fix_allele_encoding works", {
  G_test <- as.matrix(
    data.frame(
    snp1 = sample(c(0,1,2), size = 100, replace = T, prob = c(.2,.1,.7)),
    snp2 = sample(c(0,1,2), size = 100, replace = T, prob = c(.7,.1,.2))
    ))
  G_fixed <- fix_allele_encoding(G_test)
  expect_condition(length(G_fixed[snp1 == 2]) > length(G_fixed[snp1 == 0]))
  expect_condition(length(G_fixed[snp2 == 2]) < length(G_fixed[snp2 == 0]))
})
