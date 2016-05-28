context("model fitting")

test_that("Fitted object is iprior",{
	
	x <- rnorm(100); y <- rnorm(100)
	mod <- iprior(x=x, y=y, silent=T)
	expect_that(mod, is_a("iprior"))
	
})

