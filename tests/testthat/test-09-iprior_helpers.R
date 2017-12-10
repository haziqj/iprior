context("Miscellaneous tests")

test_that("Hlam", {

  y <- 1:3
  x1 <- 1:3
  x2 <- 4:6
  x3 <- factor(7:9)
  mod <- kernL(y, x1, x2, x3, kernels = c("fbm", "se", "pearson"))
  K <- get_Hlam(mod, mod$thetal$theta)
  tmp <- eigen_Hlam(K)
  expect_equal(tmp$u, c(4.066071e-16, 3.569738e+00, 4.864665e+00), tolerance = 1e-6)

})

test_that("logLik", {

  y <- 1:3
  x1 <- 1:3
  x2 <- 4:6
  x3 <- factor(7:9)
  mod <- kernL(y, x1, x2, x3, kernel = "poly", est.offset = TRUE)
  res <- loglik_iprior(c(1, 1, 1, 0, 0, 0), object = mod)
  expect_equal(res, -8.735662, tolerance = 1e-6)

})

test_that("EM loop logical",{

  environment(em_loop_logical) <- environment()
  stop.crit <- 1e-8
  maxit <- 2
	loglik <- c(NA, NA)
	niter <- 0

	expect_true(em_loop_logical())  # should be TRUE since must complete 1 iteration
	niter <- 1
	loglik[1] <- rnorm(1)
	expect_true(em_loop_logical())  # should be TRUE since maxit > 1

	niter <- 2
	loglik[2] <- loglik[1] + stop.crit * 2
	expect_true(!em_loop_logical())  # maxit is met so this is FALSE

	maxit <- 100
	expect_true(em_loop_logical())  # stop.crit not met so is TRUE

	loglik[2] <- loglik[1] + stop.crit / 2
	expect_true(!em_loop_logical())  # stop.crit is met so returns FALSE to stop

	loglik[2] <- loglik[1] - stop.crit
	expect_warning(em_loop_logical())  # decrease in log-likelihood creates warning

})

test_that("Polynomial scale parameters multiply correctly", {

  x <- y <- rnorm(3)
  lambda <- abs(rnorm(1))
  res1 <- kern_poly(x, d = 3, lam.poly = lambda)

  mod <- kernL(y, x, kernel = "poly3", lambda = lambda)
  res2 <- get_Hlam(mod, c(log(lambda), log(1)))
  res3 <- get_Htildelam(mod, c(log(lambda), log(1)), list(x))
  expect_equal(res1, res2)
  expect_equal(res2, res3)

})
