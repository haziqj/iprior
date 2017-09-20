context("Fitted and predict methods")

test_that("Fitted", {

  set.seed(123)
  y <- rnorm(3)
  x1 <- rnorm(3)
  mod <- iprior2(y, x1)
  tmp <- fitted(mod, intervals = TRUE)
  expect_equal(tmp$y, c(-0.4288420, -0.3579418, 1.5548390), tolerance = 1e-6)
  expect_equal(tmp$lower, c(-0.4677753, -0.3956874, 1.5001975), tolerance = 1e-6)
  expect_equal(tmp$upper, c(-0.3899087, -0.3201962, 1.6094805), tolerance = 1e-6)
  expect_that(print(tmp), prints_text("Training MSE:"))

})

test_that("Predict (non-formula)", {

  set.seed(123)
  y <- rnorm(3)
  x1 <- rnorm(3)
  mod <- iprior2(y, x1)
  tmp <- predict(mod, list(rnorm(1)), y.test = rnorm(1), intervals = TRUE)

  expect_equal(tmp$y, -1.342379, tolerance = 1e-6)
  expect_equal(tmp$lower, -1.408215, tolerance = 1e-6)
  expect_equal(tmp$upper, -1.276543, tolerance = 1e-6)
  expect_that(print(tmp), prints_text("Test MSE:"))
  expect_that(print(predict(mod, list(rnorm(1)))), prints_text("Test MSE: NA"))

})

test_that("Predict (formula)", {

  dat <- gen_fbm(4, seed = 123)
  mod <- iprior2(y ~ ., dat[1:3, ], kernel = "fbm")
  tmp <- predict(mod, dat[4, ], intervals = TRUE)

  expect_equal(as.numeric(tmp$y), 17.30887, tolerance = 1e-6)
  expect_equal(as.numeric(tmp$lower), 13.09142, tolerance = 1e-6)
  expect_equal(as.numeric(tmp$upper), 21.52631, tolerance = 1e-6)
  expect_that(print(tmp), prints_text("Test MSE:"))

})
