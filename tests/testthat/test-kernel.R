context("kernel matrices")

test_that("FBM with Hurst=1 equals Canonical",{
	
	x <- rnorm(100)
	expect_equal(fn.H3a(x, gamma=1), fn.H2a(x))
	
})

test_that("Pearson must take in factors",{
	
	expect_warning(fn.H1(1:100))
	
})