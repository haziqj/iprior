context("Fitted and predict methods")

test_that("Fitted", {

  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  y <- rnorm(3)
  x1 <- rnorm(3)
  mod <- iprior(y, x1, control = list(silent = TRUE))
  tmp <- fitted(mod, intervals = TRUE)
  expect_equal(tmp$y, c(-0.4288420, -0.3579418, 1.5548390), tolerance = 1e-5)
  expect_equal(tmp$lower, c(-0.5377307, -0.4555578, 1.3483344), tolerance = 1e-5)
  expect_equal(tmp$upper, c(-0.3199534, -0.2603258, 1.7613437), tolerance = 1e-5)
  expect_that(print(tmp), prints_text("Training RMSE:"))

})

test_that("Predict (non-formula)", {

  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  y <- rnorm(3)
  x1 <- rnorm(3)
  mod <- iprior(y, x1, control = list(silent = TRUE))
  tmp <- predict(mod, list(rnorm(1)), y.test = rnorm(1), intervals = TRUE)

  expect_equal(tmp$y, -1.342379, tolerance = 1e-5)
  expect_equal(tmp$lower, -1.596514, tolerance = 1e-5)
  expect_equal(tmp$upper, -1.088243, tolerance = 1e-5)
  expect_that(print(tmp), prints_text("Test RMSE:"))
  expect_that(print(predict(mod, list(rnorm(1)))), prints_text("Test RMSE: NA"))

})

test_that("Predict (formula)", {

  dat <- gen_smooth(4, seed = 123)
  mod <- iprior(y ~ ., dat[1:3, ], kernel = "fbm", control = list(silent = TRUE))
  tmp <- predict(mod, dat[4, ], intervals = TRUE)

  expect_equal(as.numeric(tmp$y), 17.30887, tolerance = 1e-6)
  expect_equal(as.numeric(tmp$lower), 17.30073, tolerance = 1e-6)
  expect_equal(as.numeric(tmp$upper), 17.317, tolerance = 1e-6)
  expect_that(print(tmp), prints_text("Test RMSE:"))

})

test_that("one.lam = TRUE", {

  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  mod <- iprior(stack.loss ~ ., stackloss, one.lam = TRUE,
                control = list(silent = TRUE))
  tmp1 <- predict(mod, stackloss[1:3, ])
  tmp2 <- predict(mod, stackloss[1:3, -4])
  expect_equal(tmp1$y, c("1" = 38.45414, "2" = 38.59653, "3" = 32.50568),
               tolerance = 1e-5)
  expect_equal(tmp2$test.error, NA)

})
