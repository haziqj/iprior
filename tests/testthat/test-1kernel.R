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
  identical(res1, res2, res3)

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
