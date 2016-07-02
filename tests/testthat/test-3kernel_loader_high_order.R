context("Kernel loader higher order")

test_that("Correct specification of higher order non-formula",{

  y <- rnorm(3)
  x1 <- rnorm(3)
  x2 <- data.frame(rnorm(3), rnorm(3), rnorm(3))
  x3 <- matrix(rnorm(3), ncol = 1)
  x4 <- x1 ^ 2
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = 1)))
  expect_error(kernL(y = y, x1 = x1, x1sq = x4, model = list(order = 1:10)))

})

test_that("Correct specification of higher order formula",{

  mod <- kernL(stack.loss ~ Air.Flow + I(Air.Flow ^ 2) + ., data = stackloss,
                model = list(order = c(1, "1^2", 2, 3)),
                control = list(silent = TRUE))
  expect_that(mod, is_a("ipriorKernel"))
  expect_equivalent(mod$r, 1)
  tmp <- summary(mod)

})
