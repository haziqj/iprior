context("Kernel loader")

test_that("Higher order non-formula",{

  y <- rnorm(3)
  x1 <- rnorm(3)
  x2 <- data.frame(rnorm(3), rnorm(3), rnorm(3))
  x3 <- matrix(rnorm(3), ncol = 1)
  x4 <- x1 ^ 2
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = 1)))
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = 1:10)))

})
