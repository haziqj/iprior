context("Cross-validation")

test_that("Training samples specified", {

  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  mod1 <- iprior(stack.loss, stack.x, train.samp = 1:10,
                 control = list(silent = TRUE))
  mod2 <- iprior(stack.loss ~ ., stackloss, one.lam = TRUE, train.samp = 1:10,
                 control = list(silent = TRUE))
  expect_equal(mod1$test, mod2$test, tol = 1e-5)

  # Check summary method
  tmp <- capture.output(print(summary(mod1)))
  tmp <- tmp[length(tmp)]
  expect_true(grepl("Test", tmp))

})

test_that("k-fold cross-validation", {

  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  mod1 <- iprior_cv(stack.loss, stack.x, par.cv = FALSE,
                    control = list(silent = TRUE, theta0 = c(1, 1)))
  mod2 <- iprior_cv(stack.loss ~ ., stackloss, one.lam = TRUE, par.cv = FALSE,
                    control = list(silent = TRUE, theta0 = c(1, 1)))
  expect_equal(apply(mod1$res, 2, mean), apply(mod2$res, 2, mean),
               tol = 1e-4)
  expect_equal(capture.output(mod1)[1], capture.output(mod2)[1])

})

test_that("Folds > 0", {

  expect_warning(
    mod2 <- iprior_cv(stack.loss ~ ., stackloss, one.lam = TRUE, folds = 1,
                      par.cv = FALSE, control = list(silent = TRUE,
                                                     theta0 = c(1, 1)))
  )

})
