context("Kernel matrices")

test_that("FBM with Hurst=1 equals Canonical",{

	x <- rnorm(100)
	# mat.fbm <- fn.H3a(x, gamma = 1); attr(mat.fbm, "class") <- NULL
	mat.fbm <- fnH3(x, gamma = 1); attr(mat.fbm, "class") <- NULL
	mat.can <- fnH2(x); attr(mat.can, "class") <- NULL
	expect_equivalent(mat.fbm, mat.can)

})

test_that("Pearson must take in factors",{

	expect_warning(fn.H1(1:100))

})

test_that("Linear kernel", {

  x <- 1:3
  y <- matrix(1:6, ncol = 2)
  y.cen <- scale(y, scale = FALSE)
  expect_true(is.kern_linear(kern_canonical(x)))
  expect_true(is.kern_linear(kern_canonical(y)))
  expect_error(kern_canonical(x, y))
  expect_error(kern_canonical(x, y, centre = FALSE))

  res1 <- kern_canonical(y, y)
  res2 <- kern_canonical(y.cen, y.cen, centre = FALSE)
  res3 <- kern_canonical(y.cen, y.cen)
  res4 <- kern_canonical(y)
  expect_true(identical(res1, res2, res3, res4))

})

test_that("Pearson kernel", {

  x <- factor(1:3)
  y <- factor(1:2)
  z <- factor(4)
  expect_true(is.kern_pearson(kern_pearson(x)))
  expect_true(is.kern_pearson(kern_pearson(x, y)))
  expect_error(kernel(x, z))
  expect_equivalent(kern_pearson(x), kern_pearson(x, x))

})

test_that("fBm kernel", {

  x <- 1:3
  y <- matrix(1:6, ncol = 2)
  expect_true(is.kern_fbm(kern_fbm(x)))
  expect_true(is.kern_fbm(kern_fbm(y)))
  expect_error(kern_fbm(x, y))

  res1 <- kern_canonical(y, y)
  res2 <- kern_fbm(y, y, gamma = 1)
  res3 <- kern_fbm(y, gamma = 1)
  expect_equivalent(res1, res2)
  expect_equivalent(res2, res3)
  expect_equivalent(res1, res3)

})

test_that("SE kernel", {

  x <- 1:3
  y <- matrix(1:6, ncol = 2)
  expect_true(is.kern_se(kern_se(x, l = 2)))
  expect_true(is.kern_se(kern_se(y)))
  expect_error(kern_se(x, y))

  res1 <- kern_se(y, y)
  res2 <- kern_se(y)
  expect_equivalent(res1, res2)

})

test_that("Polynomial kernel", {

  x <- 1:3
  y <- matrix(1:6, ncol = 2)
  expect_true(is.kern_poly(kern_poly(x, d = 3)))
  expect_true(is.kern_poly(kern_poly(y, c = 3)))
  expect_error(kern_poly(x, y))
  expect_error(kern_poly(x, d = -1))
  expect_error(kern_poly(x, c = -1))

  res1 <- kern_poly(y, y, lam.poly = 3)
  res2 <- kern_poly(y, lam.poly = 3)
  expect_equivalent(res1, res2)

})

