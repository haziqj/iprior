context("ipriorMod methods and accessors")

test_that("ipriorMod methods and accessors", {

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(y ~ ., dat, kernel = "fbm", fixed.hyp = TRUE)

  # Methods
  expect_equal(as.numeric(sigma(mod)), 1)
  expect_equal(logLik(mod), -53.06776, tolerance = 1e-5)
  expect_equal(logLik(mod), logLik(mod$ipriorKernel))
  expect_equal(deviance(mod), 106.1355, tolerance = 1e-5)
  expect_equal(deviance(mod), deviance(mod$ipriorKernel))

  # Accessor functions
  expect_equal(get_intercept(mod), 17.84021, tolerance = 1e-5)
  expect_equal(as.numeric(get_y(mod)), dat$y)
  # expect_output(get_size(mod))
  expect_equal(get_hyp(mod), c(lambda = 1, hurst = 0.5, psi = 1))
  expect_equal(get_hyp(mod), c(get_lambda(mod), get_hurst(mod), get_psi(mod)))
  expect_equal(capture.output(get_lengthscale(mod)), "NA")
  expect_equal(capture.output(get_offset(mod)), "NA")
  expect_equal(capture.output(get_degree(mod)), "NA")
  expect_equal(get_se(mod), c(lambda = NA, hurst = NA, psi = NA))
  expect_equal(get_kernels(mod), c(X = "fbm,0.5"))
  expect_equal(get_kern_matrix(mod), kern_fbm(dat$X))
  expect_equal(get_kern_matrix(mod, newdata = dat[1:2, ]),
               kern_fbm(dat$X, dat$X[1:2]))
  expect_equal(get_kern_matrix(mod$ipriorKernel), kern_fbm(dat$X))
  expect_equal(get_kern_matrix(mod$ipriorKernel, newdata = dat[1:2, ]),
               kern_fbm(dat$X, dat$X[1:2]))
  expect_equal(as.numeric(get_mse(mod)), 7.817885, tolerance = 1e-5)
  expect_equal(get_estl(mod), c(est.lambda = FALSE, est.hurst = FALSE,
                                est.lengthscale = FALSE, est.offset = FALSE,
                                est.psi = FALSE))
  expect_equal(capture.output(get_method(mod)),
               "Estimation with fixed hyperparameters.")
  expect_equal(capture.output(get_convergence(mod)),
               "Convergence not assessed.")
  expect_equal(capture.output(get_niter(mod)), "Iterations: NA/100.")
  expect_silent(get_time(mod))
  expect_true(length(get_theta(mod)) == 0)

})

test_that("ipriorMod methods and accessors (non-formula)", {

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(dat$y, dat$X, kernel = "fbm", fixed.hyp = TRUE)

  # Accessor functions
  expect_equal(get_kern_matrix(mod), kern_fbm(dat$X))
  expect_equal(get_kern_matrix(mod, newdata = list(dat$X[1:2])),
               kern_fbm(dat$X, dat$X[1:2]))
  expect_equal(get_kern_matrix(mod$ipriorKernel), kern_fbm(dat$X))
  expect_equal(get_kern_matrix(mod$ipriorKernel, newdata = list(dat$X[1:2])),
               kern_fbm(dat$X, dat$X[1:2]))

})

test_that("Nystrom methods", {

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(y ~ ., dat, kernel = "fbm", nystrom = 5,
                control = list(silent = TRUE))
  expect_equal(logLik(mod), 29.37298, tolerance = 1e-5)
  expect_equal(logLik(mod, 1:2), -32.11874, tolerance = 1e-5)

})

test_that("Training samples", {

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(y ~ ., dat, kernel = "se", train.samp = 1:5,
                control = list(silent = TRUE))
  expect_equal(as.numeric(get_mse(mod)), c(0.4624304, 21.0778450),
               tolerance = 1e-5)

})

test_that("Polynomial degree", {

  dat <- gen_smooth(10, seed = 123)
  mod <- iprior(y ~ ., dat, kernel = "poly3", control = list(silent = TRUE))
  expect_equal(as.numeric(get_degree(mod)), 3)

})
