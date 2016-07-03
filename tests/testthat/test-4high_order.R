context("Fit higher order terms")

test_that("Correct kernel specification of higher order non-formula",{

  y <- rnorm(3)
  x1 <- rnorm(3)
  x2 <- data.frame(rnorm(3), rnorm(3), rnorm(3))
  x3 <- matrix(rnorm(3), ncol = 1)
  x4 <- x1 ^ 2
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = 1)))
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = 1:10)))
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = c(1, "1:2"))))
  expect_warning(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = c(2,1))))

})

test_that("Correct kernel specification of higher order formula",{

  mod <- kernL(stack.loss ~ Air.Flow + I(Air.Flow ^ 2) + ., data = stackloss,
                model = list(order = c(1, "1^2", 2, 3)))
  expect_that(mod, is_a("ipriorKernel"))
  expect_equivalent(mod$l, 3)
  tmp <- summary(mod)

})

test_that("Successfully fit iprior with higher order and parsm = FALSE",{

  mod <- kernL(stack.loss ~ Air.Flow + I(Air.Flow ^ 2) + ., data = stackloss,
               model = list(order = c(1, "1^2", 2, 3), parsm = FALSE))
  mod.fit <- iprior(mod, control = list(silent = TRUE, maxit = 5))
  expect_equivalent(mod$l, 4)

})
