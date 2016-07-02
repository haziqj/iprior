context("I-prior methods")

test_that("Predict method works",{

  mod <- iprior(stack.loss ~ ., data = stackloss, control = list(silent = TRUE))
  res1 <- predict(mod)
  res2 <- predict(mod, stackloss[1:10, ])
  expect_equivalent(res1[1:10], res2)

})
