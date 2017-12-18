context("Kernel loader")

test_that("Kernel translator",{

  expect_true(is.kern_linear(kernel_translator(1:3, kernel = "linear")))

  res <- kernel_translator(1:3, kernel = "fbm,0.7")
  expect_true(is.kern_fbm(res))
  expect_equal(get_hyperparam(res), 0.7)

  res <- kernel_translator(1:3, kernel = "se,0.7")
  expect_true(is.kern_se(res))
  expect_equal(get_hyperparam(res), 0.7)

  res1 <- kernel_translator(1:3, kernel = "poly,0.7")
  res2 <- kernel_translator(1:3, kernel = "poly3")
  res3 <- kernel_translator(1:3, kernel = "poly4,0.9")
  expect_true(all(is.kern_poly(res1), is.kern_poly(res2), is.kern_poly(res3)))
  expect_equal(get_hyperparam(res1), 0.7)
  expect_equal(get_hyperparam(res2), 0)
  expect_equal(get_hyperparam(res3), 0.9)
  expect_equal(get_polydegree(res1), 2)
  expect_equal(get_polydegree(res2), 3)
  expect_equal(get_polydegree(res3), 4)

})

test_that("Kernel to param to theta", {

  # Remember that kernels vector need to have the hyperparameters specified,
  # i.e. fbm,0.5 instead of just fbm, etc.
  # e.g. no of theta = 5 + 1 + 1 + 1 = 8
  # param table is 5 x 5
  kernels <- c("linear", "fbm,0.7", "se,2", "poly3,0.15", "pearson")
  which.pearson <- c(F, F, F, F, T)
  poly.degree <- c(NA, NA, NA, 3, NA)
  lambda <- rep(1, 5)
  psi <- 1
  names(psi) <- "psi"

  res1 <- kernel_to_param(kernels, lambda)
  res2 <- param_to_theta(res1, list(est.lambda      = c(FALSE, TRUE),
                                    est.hurst       = TRUE,
                                    est.lengthscale = FALSE,
                                    est.offset      = TRUE,
                                    est.psi         = TRUE))
  res3 <- theta_to_param(res2$theta, list(thetal        = res2,
                                          which.pearson = which.pearson,
                                          poly.degree   = poly.degree))

  expect_equal(kernels, res3$kernels)
  expect_equal(length(collapse_param(res3)$param), 8)
  expect_equal(theta_to_psi(res2$theta, list(thetal = res2)), psi)
  expect_equal(res1, res3)

})

test_that("Incorrect specification of kernels", {

  y <- x <- 1:3
  # Not enough kernels
  expect_warning(mod <- kernL(y, x, x, x, kernel = c("poly", "se")))
  # Too many kernels
  expect_warning(mod <- kernL(y, x, kernel = c("poly", "se")))

})

test_that("Interactions", {

  y <- x <- 1:3
  mod1 <- kernL(y, x, x, x, interactions = c("1:2", "1:2:3"))
  mod2 <- kernL(stack.loss ~ . ^ 2, stackloss)
  expect_true(is.ipriorKernel(mod1))
  expect_equal(mod1$no.int, 1)
  expect_equal(mod1$no.int.3plus, 1)
  expect_true(is.ipriorKernel(mod2))
  expect_equal(mod2$no.int, 3)
  expect_equal(mod2$no.int.3plus, 0)

})

test_that("Fixed hyperparameter option", {

  y <- x <- 1:3
  mod <- kernL(y, x, kernel = "fbm,0.8", fixed.hyp = TRUE)
  expect_true(is.ipriorKernel(mod))
  expect_false(all(unlist(mod$estl)))

})

test_that("Formula input", {

  dat <- gen_smooth(3)
  mod1 <- kernL(y ~ ., dat, kernel = "poly")
  mod2 <- kernL(dat$y, dat$X, kernel = "poly")
  tmp1 <- capture.output(print(mod1))
  tmp2 <- capture.output(print(mod2))
  expect_equal(tmp1[-3], tmp2[-3])  # object size does not match

})

test_that("get_Hlam()", {

  y <- x1 <- 1:3
  x2 <- factor(1:3)
  x3 <- 4:6
  mod1 <- kernL(y, x1, x2, x3, kernel = c("se", "pearson", "fbm"))
  theta1 <- mod1$thetal$theta
  mod2 <- kernL(y, x1, x2, x3, kernel = c("se", "pearson", "fbm"),
                 est.hurst = TRUE, est.lengthscale = TRUE)
  theta2 <- mod2$thetal$theta
  expect_equivalent(get_Hlam(mod1, theta1, theta.is.lambda = TRUE),
                    get_Hlam(mod2, theta2))

})

test_that("Training samples specified", {

  # Non-formula
  y <- x <- 1:10
  mod <- kernL(y, x, train.samp = 1:5)
  expect_equal(mod$n, 5)
  expect_equal(mod$y.test, mod$Xl.test[[1]])

  # Formula
  mod <- kernL(y ~ x, data.frame(y, x), train.samp = 1:5)
  expect_equal(mod$n, 5)
  expect_equivalent(mod$y.test, mod$Xl.test[[1]])

})

test_that("Test samples specified", {

  # Non-formula
  y <- x <- 1:10
  mod <- kernL(y, x, test.samp = 1:5)
  expect_equal(mod$n, 5)
  expect_equal(mod$y.test, mod$Xl.test[[1]])

  # Formula
  mod <- kernL(y ~ x, data.frame(y, x), test.samp = 1:5)
  expect_equal(mod$n, 5)
  expect_equivalent(mod$y.test, mod$Xl.test[[1]])

})

test_that("Warn if training samples not correct", {

  y <- x <- 1:10
  expect_warning(kernL(y, x, train.samp = 1:200))
  expect_warning(kernL(y, x, train.samp = 250))

})

test_that("Stop if train.samp and test.samp both specified", {

  y <- x <- 1:10
  expect_error(kernL(y, x, train.samp = 1:5, test.samp = 1:5))

})

# test_that("Mixed kernels", {
#
#   expect_warning(
#     mod <- kernL(stack.loss ~ ., data = stackloss,
#                  model = list(order = c(1,"1^2",2),
#                               kernel = c("FBM", "Canonical", "FBM")))
#   )
#
# })
#
# test_that("Mixed kernels multiple Hurst", {
#
#   mod1 <- kernL(stack.loss ~ ., stackloss,
#                 model = list(kernel = c("FBM,0.1", "Canonical", "FBM,0.1")))
#   mod2 <- kernL(stack.loss ~ ., stackloss,
#                 model = list(kernel = c("FBM", "Canonical", "FBM"), Hurst = 0.1))
#   expect_equivalent(mod1$Hl, mod2$Hl)
#   expect_warning(
#     kernL(stack.loss ~ ., stackloss,
#           model = list(kernel = c("FBM,0.9", "Canonical", "FBM,0.9"), Hurst = 0.1))
#   )
#
# })
#
# test_that("Warn when using one.lam = TRUE", {
#
#   expect_warning(
#     mod <- kernL(stack.loss ~ ., data = stackloss,
#                  model = list(one.lam = TRUE,
#                               kernel = c("FBM", "Canonical", "FBM")))
#   )
#
# })
#
# test_that("Automatic and manual are the same", {
#
#   mod1 <- kernL(len ~ ., data = ToothGrowth,
#                 model = list(one.lam = TRUE))
#   toot <- ToothGrowth
#   toot$supp <- as.numeric(ToothGrowth$supp)
#   suppressWarnings(mod2 <- kernL(len ~ ., data = toot,
#                                  model = list(one.lam = TRUE,
#                                               kernel = c("Pearson", "Canonical"))))
#   expect_equivalent(mod1$Hl, mod2$Hl)
#
# })
