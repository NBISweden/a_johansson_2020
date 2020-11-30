test_that("simulating effect sizes works for common variants", {
  require('GenABEL')
  data(srdta)
  maf <- seq(from = 0.001, to = 0.9, by=0.005)
  eff1 <- get_effects(maf,
              N = 5,
              get_betas_fun=dbeta,
              get_betas_args=list(shape1 = 0.1, shape2 = 0.1),
              thr=0.01,
              rare=F,
              frac_negative=.4,
              seed = 42)
    expect_equal(length(eff1$marker_idx), 5)
    expect_true(all(maf[eff1$marker_idx] > 0.01))
    expect_lt(max(eff1$effects), .3)
    expect_lt(min(eff1$effects), 0)
})

test_that("simulating effect sizes works for rare variants", {
  require('GenABEL')
  data(srdta)
  maf <- seq(from = 0.001, to = 0.9, by=0.005)
  eff2 <- get_effects(maf,
                      N = 5,
                      get_betas_fun=dbeta,
                      get_betas_args=list(shape1 = 1, shape2 = 25),
                      thr=0.05,
                      rare=T,
                      frac_negative=.4,
                      seed = 47)
  expect_equal(length(eff2$marker_idx), 5)
  expect_lt(max(eff2$effects), 25)
  expect_lt(min(eff2$effects), 0)
  expect_true(all(maf[eff2$marker_idx] < 0.05))
})

test_that("warnings work", {
  require('GenABEL')
  data(srdta)
  maf <- seq(from = 0.004, to = 1.1, by=0.005)

  expect_warning(get_effects(maf,
                      N = 5,
                      thr=0.005,
                      rare=T,
                      frac_negative=.4,
                      seed = 47))

  expect_warning(get_effects(maf,
                             N = 5,
                             thr=0.0001,
                             rare=T,
                             frac_negative=.4,
                             seed = 47))

  get_betas_one <- function(x) {
    return(1)
  }
  maf <- seq(from = 0.001, to = 0.5, by=0.01)
  effs <- get_effects(maf, N = 5, get_betas=get_betas_one, get_betas_args = list(), thr=0.01, rare=F, frac_negative=.4, seed = 47)
  expect_equal(range(effs$effects), c(-1, 1))
})
