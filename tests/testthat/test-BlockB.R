context("Block B stuff")

test_that("BlockB2()", {

  H1 <- kern_canonical(1:3)
  H2 <- kern_canonical(4:6)
  n <- 3
  intr <- matrix(c(1, 2), ncol = 1)
  expect_equal(BlockB_fn(list(H1), intr = NULL, n = n, p = 1)$BB.msg,
               "Single lambda")
  expect_equal(BlockB_fn(list(H1, H2), intr = NULL, n = n, p = 2)$BB.msg,
               "Multiple lambda with no interactions")
  expect_equal(BlockB_fn(list(H1, H2), intr = intr, n = n, p = 2)$BB.msg,
               "Multiple lambda with parsimonious interactions")

})

test_that("BlockB creation conditions", {

  y <- 1:3
  x1 <- 1:3
  x2 <- 4:6
  expect_true(!is.null(kernL2(y, x1, kernel = "se")$BlockBStuff))  # single lambda
  expect_true(!is.null(kernL2(y, x1, x2, kernel = "fbm")$BlockBStuff))  # mult lambda
  expect_true(!is.null(kernL2(y, x1, x2, interactions = "1:2")$BlockBStuff))
  expect_true(
    is.null(kernL2(y, x1, x2, kernel = c("poly2", "fbm"))$BlockBStuff)
  )  # polynomial kernel
  expect_true(
    is.null(kernL2(y, x1, x2, kernel = "fbm", est.hurst = TRUE)$BlockBStuff)
  )  # est hurst
  expect_true(
    is.null(kernL2(y, x1, x2, kernel = "se", est.lengthscale = TRUE)$BlockBStuff)
  )  # est lengthscale
  expect_true(
    is.null(kernL2(y, x1, x2, est.lambda = FALSE, est.psi = FALSE)$BlockBStuff)
  )  # no estimate lambda and psi
  expect_true(
    !is.null(kernL2(y, x1, x2, est.lambda = TRUE, est.psi = FALSE)$BlockBStuff)
  )  # no estimate lambda
  expect_true(
    !is.null(kernL2(y, x1, x2, est.lambda = FALSE, est.psi = TRUE)$BlockBStuff)
  )  # no estimate psi

})
