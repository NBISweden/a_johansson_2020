test_that("check the data structure is S4 and holds", {
  test_data <<- tibble(
    chr = c(1),
    pos = c(1001,1002,1005,1017),
    snp = c('snp1', 'snp2', 'snp3', 'snp4'),
    strand = c('u'),
    ind1 = c(0,1,1,2),
    id2 = c(2,NA,1,1),
    john = c(0,0,1,0),
    jill = c(1,1,1,2),
    elon = c(2,0,NA,1)
  )
  sex <<- c(0,1,1,0,1)
  trait <<- c(1.3,3.4,4.4,2.2,1.7)
  result <<- tibble_to_gwaa(x = test_data, sex = sex, trait = trait)
  expect_type(result, "S4")
  expect_s4_class(result, 'gwaa.data')
  expect_equal(GenABEL::nids(result), 5)
  expect_equal(GenABEL::nsnps(result), 4)
  expect_equal(result@phdata[3,1], 'john')
  expect_type(result@phdata[,2], 'double')
})
