context("iprior2 and estimation methods")

test_that("iprior_direct", {

  y <- 1:3
  x1 <- 1:3
  x2 <- factor(7:9)
  mod <- kernL2(y, x1, x2)
  suppressWarnings(
    res <- iprior_direct(mod, loglik_iprior, 1:3, list(fnscale = -2, trace = 0,
                                                       maxit = 1))
  )
  theta <- c(0.4311141, -0.2127298, -0.3132339)
  names(theta) <- names(res$theta)
  expect_equal(res$theta, theta, tolerance = 1e-6)
  expect_equal(res$loglik[length(res$loglik)], -4.051600, tolerance = 1e-6)

})

test_that("iprior_fixed", {

  mod <- iprior2(stack.loss ~ ., stackloss, fixed.hyp = TRUE,
                 control = list(silent = TRUE))
  expect_equal(as.numeric(mod$param.full), rep(1, 4))

})

test_that("iprior_em_closed", {

  set.seed(123)
  mod <- iprior2(kernL2(stack.loss ~ ., stackloss), method = "em",
                 control = list(maxit = 3, silent = TRUE))
  expect_equal(as.numeric(mod$param.full),
               c(0.12818, 1.68720, 0.25099, 0.13146), tolerance = 1e-5)

})

test_that("print()", {

  y <- 1:3
  x1 <- 1:3
  x2 <- factor(7:9)
  mod <- iprior2(y, x1, x2, fixed.hyp = TRUE, control = list(silent = TRUE))
  tmp <- capture.output(print(mod))
  tmp <- summary(mod)
  tmp <- capture.output(print(tmp))

})
