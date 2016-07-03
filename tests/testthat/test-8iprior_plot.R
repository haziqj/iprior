context("Plots")

test_that("Plot for FBM kernel",{

  data("datfbm")
  mod <- iprior(y ~ x, data = datfbm, model = list(kernel = "FBM"),
                control = list(silent = TRUE))
  expect_that(mod, is_a("ipriorMod"))
  tmp <- plot(mod, plots = "fitted")
  tmp <- plot(mod, plots = "resid")
  tmp <- plot(mod, plots = "qqplot")

})

test_that("Multilevel plot",{

  data("simdat")
  mod <- iprior(y ~ . ^ 2, data = simdat,
                control = list(silent = TRUE, lambda = c(0.47012736, 0.02357745),
                               psi = 3.61407951))
  expect_that(mod, is_a("ipriorMod"))
  tmp <- plot(mod, plots = "fitted")
  tmp <- plot(mod, plots = "resid")
  tmp <- plot(mod, plots = "qqplot")

})

test_that("Unable to plot if x dim > 1",{

  mod <- iprior(stack.loss ~ ., stackloss, control = list(silent = TRUE,
                                                          maxit = 5))
  expect_message(plot(mod, plots = "fitted"))

})
