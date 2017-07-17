context("Nystrom approximation")

test_that("Reordering function is the same as manual reordering",{

  set.seed(123)
  n <- 5
  X <- runif(n, -10, 10)
  y <- 1 + 3 * X + rnorm(n, sd = 10)
  set.seed(123)
  Nys.samp <- sample(1:5)

  dat <- data.frame(y, X)
  dat2 <- data.frame(y = y[Nys.samp], X = X[Nys.samp])

  mod.manual <- kernL(y ~ X + I(X^2) + X:I(X^2), dat2)
  suppressWarnings(
    mod <- iprior(y ~ X + I(X^2) + X:I(X^2), dat,
                  control = list(Nystrom = 4, Nys.seed = 123, maxit = 1,
                                 silent = TRUE))
  )
  mod.iprior <- mod$ipriorKernel

  expect_equal(mod.manual$x, mod.iprior$x)
  expect_equal(mod.manual$Hl, mod.iprior$Hl)
  expect_equal(mod.manual$BlockBstuff, mod.iprior$BlockBstuff)

})

