context("model fitting")

test_that("Fitted object is iprior",{

	x <- rnorm(100); y <- rnorm(100)
	mod <- iprior(y=y, x=x, control=list(silent=T))
	expect_that(mod, is_a("iprior"))

})

test_that("Successfully fit interactions",{

	data(cats, package="MASS")
	mod1 <- iprior(Hwt ~ Bwt * Sex, cats, control=list(silent=T))
	mod2 <- iprior(y=cats$Hwt, Bwt=cats$Bwt, Sex=cats$Sex, model=list(interactions="1:2"), control=list(silent=T))
	expect_equal(mod1$log, mod2$log) #check if the same log-likelihood value achieved

})
