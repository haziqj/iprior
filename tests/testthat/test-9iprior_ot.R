context("Optim and theta")

test_that("iprior with theta and lambda and psi",{

  expect_error(
    mod <- iprior(stack.loss ~ ., stackloss,
                  control = list(lambda = rnorm(3), psi = abs(rnorm(1)),
                                 theta = abs(rnorm(4))))
  )

})

test_that("iprior with theta and lambda",{

  expect_error(
    mod <- iprior(stack.loss ~ ., stackloss,
                  control = list(lambda = rnorm(3),
                                 theta = abs(rnorm(4))))
  )

})

test_that("iprior with theta and psi",{

  expect_error(
    mod <- iprior(stack.loss ~ ., stackloss,
                  control = list(psi = abs(rnorm(1)),
                                 theta = abs(rnorm(4))))
  )

})

test_that("iprior with theta and sigma",{

  expect_error(
    mod <- iprior(stack.loss ~ ., stackloss,
                  control = list(sigma = abs(rnorm(1)),
                                 theta = abs(rnorm(4))))
  )

})

test_that("iprior with psi and sigma",{

  expect_error(
    mod <- iprior(stack.loss ~ ., stackloss,
                  control = list(psi = abs(rnorm(1)),
                                 sigma = abs(rnorm(1))))
  )

})

test_that("iprior with lambda, psi and sigma",{

  expect_error(
    mod <- iprior(stack.loss ~ ., stackloss,
                  control = list(lambda = rnorm(3), psi = abs(rnorm(1)),
                                 sigma = abs(rnorm(1))))
  )

})

test_that("iprior and optim wrapper",{

  mod <- kernL(stack.loss ~ . ^ 2, data = stackloss)
  mod.fit <- ipriorOptim(mod, control = list(silent = TRUE))
  expect_that(mod.fit, is_a("ipriorMod"))

})
