context("Mixed kernels")

test_that("Mixed kernels", {

  expect_warning(
    mod <- kernL(stack.loss ~ ., data = stackloss,
                 model = list(order = c(1,"1^2",2),
                              kernel = c("FBM", "Canonical", "FBM")))
  )

})

test_that("Warn when using one.lam = TRUE", {

  expect_warning(
    mod <- kernL(stack.loss ~ ., data = stackloss,
                 model = list(one.lam = TRUE,
                              kernel = c("FBM", "Canonical", "FBM")))
  )

})


