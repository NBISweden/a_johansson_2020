test_that("getting maf works", {
  require('GenABEL')
  data(srdta)
  expect_warning(get_maf(srdta))
  maf <- suppressWarnings(get_maf(srdta))
  expect_true(near(maf[1], 0.1325503356))
  expect_error(get_maf(srdta@gtdata))
})
