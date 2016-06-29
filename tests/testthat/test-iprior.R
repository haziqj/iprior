context("model fitting")

test_that("Fitted object is iprior",{

	mod <- iprior(y = rnorm(100), x = rnorm(100), control = list(silent = TRUE))
	expect_that(mod, is_a("ipriorMod"))

})

test_that("Successfully fit interactions",{

	mod1 <- iprior(len ~ supp * dose, ToothGrowth, control = list(silent = TRUE))
	mod2 <- iprior(y = ToothGrowth$len, supp = ToothGrowth$supp,
	               dose = ToothGrowth$dose,
	               model = list(interactions = "1:2"),
	               control = list(silent = TRUE))
	expect_equal(mod1$log, mod2$log)  # check if same log-likelihood value achieved

})
