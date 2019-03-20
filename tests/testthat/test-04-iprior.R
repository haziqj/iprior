context("Estimation method checker")

test_that("Fixed", {

  mod <- list(thetal = list(n.theta = 0))
  res <- iprior_method_checker(mod, "fixed")
  expect_true(res["fixed"])

  mod <- list(thetal = list(n.theta = 10))
  res <- iprior_method_checker(mod, "fixed")
  expect_true(res["fixed"])

  # mod <- list(thetal = list(n.theta = 0))
  # expect_warning(res <- iprior_method_checker(mod, "em"))
  # expect_true(res["fixed"])

})

test_that("Canonical", {

  mod <- list(thetal = list(n.theta = 10), kernels = "linear", no.int = 0)
  res <- iprior_method_checker(mod, "canonical")
  expect_true(res["canonical"])

  mod <- list(thetal = list(n.theta = 10), kernels = "se,1", no.int = 0)
  expect_warning(iprior_method_checker(mod, "canonical"))
  expect_true(suppressWarnings(
    iprior_method_checker(mod, "canonical")["direct"]
  ))

  mod <- list(thetal = list(n.theta = 10), kernels = "linear", no.int = 1)
  expect_warning(iprior_method_checker(mod, "canonical"))
  expect_true(suppressWarnings(
    iprior_method_checker(mod, "canonical")["direct"]
  ))

})

test_that("EM", {

  mod <- list(thetal = list(n.theta = 10), BlockBStuff = 1,
              estl = list(est.lambda = TRUE, est.psi = TRUE))
  res <- iprior_method_checker(mod, "em")
  expect_true(res["em.closed"])

  mod <- list(thetal = list(n.theta = 10), BlockBStuff = NULL)
  res <- iprior_method_checker(mod, "em")
  expect_true(res["em.reg"])

})

test_that("Direct", {

  mod <- list(thetal = list(n.theta = 10))
  res <- iprior_method_checker(mod, "direct")
  expect_true(res["direct"])

})

test_that("Nystrom", {

  mod <- structure(list(thetal = list(n.theta = 10), nystroml = 1),
                   class = "ipriorKernel")
  res <- iprior_method_checker(mod, "direct")
  expect_true(res["nystrom"])

})

context("Various estimation methods")

test_that("iprior_canonical", {

  suppressWarnings(RNGversion("3.5.0"))
  mod <- kernL(y ~ ., gen_smooth(n = 5, seed = 123))
  suppressWarnings({
    mod1 <- iprior(mod, method = "canonical",
                   control = list(silent = TRUE, maxit = 20, theta0 = 1:2))
    mod2 <- iprior_canonical(mod, 1:2, list(trace = 0, maxit = 19))
  })
  expect_equal(as.numeric(get_hyp(mod1)), c(0.04241322, 0.05310080),
               tolerance = 1e-5)
  expect_equal(get_hyp(mod1), mod2$param.full, tolerance = 1e-5)

})

test_that("iprior_direct", {

  y <- 1:3
  x1 <- 1:3
  x2 <- factor(7:9)
  mod <- kernL(y, x1, x2)
  suppressWarnings({
    mod1 <- iprior(mod, control = list(silent = TRUE, maxit = 0, theta0 = 1:3))
    mod2 <- iprior_direct(mod, loglik_iprior, 1:3, list(trace = 0, maxit = 0))
  })
  expect_equal(as.numeric(get_hyp(mod1)), c(-0.2031616, -2.2147428, 1.8114034),
               tolerance = 1e-5)
  expect_equal(get_hyp(mod1), mod2$param.full, tolerance = 1e-5)

})

test_that("iprior_fixed", {

  mod <- kernL(stack.loss ~ ., stackloss, kernel = "se,3", fixed.hyp = TRUE)
  mod1 <- iprior(mod, control = list(silent = TRUE))
  mod2 <- iprior_fixed(mod)
  expect_equal(as.numeric(get_hyp(mod1)), c(1, 1, 1, 3, 3, 3, 1))
  expect_equal(get_hyp(mod1), mod2$param.full, tolerance = 1e-5)

})

test_that("iprior_em_closed", {

  mod <- kernL(stack.loss ~ ., stackloss)
  suppressWarnings({
    mod1 <- iprior(mod, method = "em",
                   control = list(maxit = 2, silent = TRUE, theta0 = 1:4))
    mod2 <- iprior_em_closed(mod, maxit = 2, silent = TRUE, theta0 = 1:4)
  })
  expect_equal(as.numeric(get_hyp(mod1)),
               c(0.9994840, 1.9973584, 2.9036381, 0.4840088), tolerance = 1e-5)
  expect_equal(get_hyp(mod1), mod2$param.full, tolerance = 1e-5)

})

test_that("iprior_em_reg", {

  mod <- kernL(stack.loss ~ Air.Flow, stackloss, kernel = "poly3")
  suppressWarnings({
    mod1 <- iprior(mod, method = "em",
                   control = list(maxit = 2, silent = TRUE, theta0 = 1:2))
    mod2 <- iprior_em_reg(mod, maxit = 2, silent = TRUE, theta0 = 1:2)
  })
  # expect_equal(as.numeric(get_hyp(mod1)), c(8.726267e-09, 0, 4.338007e-03),
  #              tolerance = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tolerance = 1e-5)

})

test_that("iprior_em_mixed", {

  mod <- kernL(y ~ . ^ 2, gen_multilevel(10, 3, seed = 123))
  suppressWarnings({
    mod1 <- iprior(mod, method = "mixed",
                   control = list(maxit = 0, silent = TRUE, theta0 = 1:3))
    mod2 <- iprior_mixed(mod, silent = TRUE, theta0 = 1:3,
                         control.optim = list(maxit = 0))
  })
  expect_equal(as.numeric(get_hyp(mod1)), c(0.7378479, 1.7849915, 0.2541533),
               tolerance = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tolerance = 1e-5)

})

test_that("iprior_nystrom", {

  suppressWarnings(RNGversion("3.5.0"))
  mod <- kernL(y ~ ., gen_smooth(100, seed = 123), kernel = "fbm,0.7",
                nystrom = 10)
  suppressWarnings({
    mod1 <- iprior(mod, control = list(maxit = 0, silent = TRUE, theta0 = 1:2))
    mod2 <- iprior_direct(mod, loglik_nystrom, 1:3, list(trace = 0, maxit = 0))
  })
  expect_equal(as.numeric(get_hyp(mod1)), c(6.663954, 0.7, 4.746405),
               tolerance = 1e-5)
  expect_equal(mod1$param.full, mod2$param.full, tolerance = 1e-5)

})

# test_that("iprior_parallel", {  # Can't do this test in travis-ci or appveyor
#
#   mod <- kernL(y ~ ., gen_smooth(10, seed = 123), kernel = "fbm")
#   suppressWarnings({
#     mod1 <- iprior(mod, control = list(silent = TRUE, restarts = 1))
#     expect_message(
#       mod2 <- iprior_parallel(mod, control = list(theta0 = 1:2, silent = TRUE,
#                                                   restarts = 1))
#     )
#   })
#   expect_equal(logLik(mod1), -25.4149, tolerance = 1e-5)
#   expect_equal(logLik(mod1), logLik(mod2), tolerance = 1e-5)
#
# })

test_that("update and iprior.ipriorMod", {

  mod <- iprior(y ~ ., gen_smooth(30, seed = 123), kernel = "se", method = "em",
                control = list(maxit = 10, silent = TRUE))
  expect_message(mod <- iprior(mod, iter.update = 10),
                 "Updating iprior model with 10")
  expect_message(update(mod, iter.update = 10), "Updating iprior model with 10")
  expect_true(is.ipriorMod(mod))

})

context("Print and summary for ipriorMod")

test_that("print()", {

  y <- 1:3
  x1 <- 1:3
  x2 <- factor(7:9)
  mod <- iprior(y, x1, x2, fixed.hyp = TRUE, control = list(silent = TRUE))
  tmp <- capture.output(print(mod))
  expect_equal(tmp, tmp)
  tmp <- summary(mod)
  tmp <- capture.output(print(tmp))
  expect_equal(tmp, tmp)

})
