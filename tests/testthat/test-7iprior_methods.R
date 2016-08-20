context("ipriorMod methods")

test_that("Methods for regular fit",{

  mod <- iprior(stack.loss ~ ., data = stackloss, control = list(silent = TRUE))
  res1 <- predict(mod)
  res2 <- predict(mod, stackloss[1:10, ])
  expect_equivalent(res1[1:10], res2)
  expect_equal(mod$log.lik, logLik(mod))
  expect_equal(-2*mod$log.lik, deviance(mod))

})

test_that("Methods for interactions fit",{

  mod <- iprior(len ~ . ^ 2, data = ToothGrowth, control = list(silent = TRUE))
  res1 <- predict(mod)
  res2 <- predict(mod, ToothGrowth[c(9, 39, 59), ])
  expect_equivalent(res1[c(9, 39, 59)], res2)
  expect_equal(mod$log.lik, logLik(mod))
  expect_equal(-2*mod$log.lik, deviance(mod))

})

test_that("Methods for one.lam = TRUE fit",{

  data(pollution)
  mod <- iprior(Mortality ~ ., data = pollution, model = list(one.lam = TRUE),
                control = list(silent = TRUE))
  res1 <- predict(mod)
  x <- sample(1:60, 10, replace = FALSE)
  res2 <- predict(mod, pollution[x, ])
  expect_equivalent(res1[x], res2)
  expect_equal(mod$log.lik, logLik(mod))
  expect_equal(-2*mod$log.lik, deviance(mod))

})

test_that("Methods for higher order fit",{

  mod <- iprior(stack.loss ~ Air.Flow + I(Air.Flow ^ 2) + ., data = stackloss,
                model = list(order = c(1, "1^2", 2, 3)),
                control = list(silent = TRUE, maxit = 5))
  res1 <- predict(mod)
  res2 <- predict(mod, stackloss[1:10, ])
  expect_equivalent(res1[1:10], res2)
  expect_equal(mod$log.lik, logLik(mod))
  expect_equal(-2*mod$log.lik, deviance(mod))

})

test_that("Progress function",{

  mod <- iprior(len ~ . ^ 2, ToothGrowth, control = list(progress = "none"))
  tmp <- progress(mod, 50)
  expect_that(tmp, is_a("data.frame"))

})

test_that("sigma method",{

  mod <- iprior(stack.loss ~ ., data = stackloss, control = list(silent = TRUE,
                                                                 maxit = 1))
  expect_that(sigma(mod), is_a("numeric"))

})
