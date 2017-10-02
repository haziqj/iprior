context("iprior2")

test_that("Hlam", {

  y <- 1:3
  x1 <- 1:3
  x2 <- 4:6
  x3 <- factor(7:9)
  mod <- kernL2(y, x1, x2, x3, kernels = c("fbm", "se", "pearson"))
  K <- get_Hlam(mod, mod$thetal$theta)
  tmp <- eigen_Hlam(K)
  expect_equal(det(K), 32.73602, tolerance = 1e-6)
  expect_equal(tmp$u, c(1.869905, 3.598763, 4.864665), tolerance = 1e-6)

})

test_that("logLik", {

  y <- 1:3
  x1 <- 1:3
  x2 <- 4:6
  x3 <- factor(7:9)
  mod <- kernL2(y, x1, x2, x3, kernel = "poly", est.offset = TRUE)
  res <- loglik_iprior(c(1, 1, 1, 0, 0, 0), object = mod)
  expect_equal(res, -8.735662, tolerance = 1e-6)

})

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

  mod <- iprior2(stack.loss ~ ., stackloss, fixed.hyp = TRUE)
  expect_equal(as.numeric(mod$param.full), rep(1, 4))

})
