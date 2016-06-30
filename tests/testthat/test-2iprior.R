context("model fitting")

test_that("Fit using non-formula",{

	mod <- iprior(y = rnorm(100), x = rnorm(100),
	              control = list(silent = TRUE, maxit = 5))
	expect_that(mod, is_a("ipriorMod"))

})

test_that("Fit using formula",{

  mod <- iprior(stack.loss ~ ., data = stackloss,
                control = list(silent = TRUE, maxit = 5))
  expect_that(mod, is_a("ipriorMod"))

})

test_that("Fit ipriorKernel object",{

  mod <- kernL(stack.loss ~ ., data = stackloss, model = list())
  mod <- iprior(mod, control = list(silent = TRUE, maxit = 5))
  expect_that(mod, is_a("ipriorMod"))

})

test_that("Fit/update ipriorMod object",{

  mod <- iprior(stack.loss ~ ., data = stackloss,
                control = list(silent = TRUE, maxit = 5))
  iprior(mod, control = list(silent = TRUE, maxit = 5))
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
