maf <- c(0.1,0.2,0.3)
test_that("validate maf works for well-behaved maf values", {
  tmp <- validate_maf(maf)
  expect_null(tmp[[1]])
  expect_null(tmp[[2]])
  expect_equal(tmp[[3]], maf)
})

maf <- c(0.1,0.2,0.7)
tmp <- suppressWarnings(validate_maf(maf))
test_that("validate maf works for ill-behaved maf values, greater than 0.5", {
  expect_warning(validate_maf(maf))
  expect_false(is.null(tmp[[1]]))
  expect_equal(tmp[[2]], 0.7)
  expect_equal(tmp[[3]], pmin(maf, 1-maf))
})
