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
