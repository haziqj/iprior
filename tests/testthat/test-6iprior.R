context("Model fitting")

test_that("Fit using non-formula",{

  y <- rnorm(100)
  x1 <- rnorm(100)
  x2 <- matrix(rnorm(100), ncol = 1)
  x3 <- matrix(rnorm(200), ncol = 2)
  x4 <- as.data.frame(x3)
	mod1 <- iprior(y = y, x = x1,
	              control = list(silent = TRUE, maxit = 5))
	mod2 <- iprior(y = y, x1 = x1, x2 = x2,
	               control = list(silent = TRUE, maxit = 5))
	mod3 <- iprior(y = y, x1 = x1, x2 = x2, x3 = x3,
	               control = list(silent = TRUE, maxit = 5))
	mod4 <- iprior(y = rnorm(100), x = x4,
	               control = list(silent = TRUE, maxit = 5))
	expect_that(mod1, is_a("ipriorMod"))
	expect_that(mod2, is_a("ipriorMod"))
	expect_that(mod3, is_a("ipriorMod"))
	expect_that(mod4, is_a("ipriorMod"))

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
  tmp <- summary(mod1)
  tmp <- summary(mod2)

})

test_that("Successfully fit non-parsimonious interactions",{

  mod1 <- iprior(len ~ supp * dose, ToothGrowth,
                 model = list(parsm = FALSE), control = list(silent = TRUE))
  mod2 <- iprior(y = ToothGrowth$len, supp = ToothGrowth$supp,
                 dose = ToothGrowth$dose,
                 model = list(interactions = "1:2", parsm = FALSE),
                 control = list(silent = TRUE))
  expect_equal(mod1$log, mod2$log)  # check if same log-likelihood value achieved
  tmp <- summary(mod1)
  tmp <- summary(mod2)

})

test_that("Successfully fit higher order terms",{

  mod <- iprior(stack.loss ~ Air.Flow + I(Air.Flow ^ 2) + ., data = stackloss,
                model = list(order = c(1, "1^2", 2, 3)),
                control = list(silent = TRUE,
                               lambda = c(0.001549418, 0.215, -0.0077),
                               psi = 0.110104468))
  expect_that(mod, is_a("ipriorMod"))
  expect_equivalent(mod$ipriorKernel$r, 1)
  tmp <- summary(mod)

})
