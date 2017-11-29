context("Plots")

# test_that("Plot for FBM kernel",{
#
#   data("datfbm")
#   mod <- iprior(y ~ x, data = datfbm, model = list(kernel = "FBM"),
#                 control = list(silent = TRUE))
#   expect_that(mod, is_a("ipriorMod"))
#   tmp <- plot(mod, plots = "fitted")
#   tmp <- plot(mod, plots = "resid")
#   tmp <- plot(mod, plots = "qqplot")
#
# })
#
# test_that("Multilevel plot",{
#
#   data("simdat")
#   mod <- iprior(y ~ . ^ 2, data = simdat,
#                 control = list(silent = TRUE, lambda = c(0.47012736, 0.02357745),
#                                psi = 3.61407951))
#   expect_that(mod, is_a("ipriorMod"))
#   tmp <- plot(mod, plots = "fitted")
#   tmp <- plot(mod, plots = "resid")
#   tmp <- plot(mod, plots = "qqplot")
#
# })
#
# test_that("Unable to plot if x dim > 1",{
#
#   mod <- iprior(stack.loss ~ ., stackloss, control = list(silent = TRUE,
#                                                           maxit = 5))
#   expect_message(plot(mod, plots = "fitted"))
#
# })

test_that("Plots", {

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(y ~ ., dat, kernel = "fbm", method = "em",
                control = list(silent = TRUE, maxit = 10))
  expect_silent(p <- plot(mod))
  expect_silent(p <- plot_fitted(mod))
  expect_silent(p <- plot_resid(mod))
  expect_silent(p <- plot_iter(mod))
  expect_silent(p <- plot_ppc(mod, draws = 2))

  dat <- gen_multilevel(n = 10, m = 2, seed = 123)
  mod <- iprior(y ~ ., dat, method = "em",
                control = list(silent = TRUE, maxit = 10))
  expect_silent(p <- plot_fitted_multilevel(mod))
  expect_silent(p <- plot_fitted_multilevel(mod, extrapolate = TRUE))

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(y ~ ., dat, kernel = "fbm", method = "em", nystrom = 5,
                control = list(silent = TRUE, maxit = 10))
  expect_silent(p <- plot_fitted(mod))

})
