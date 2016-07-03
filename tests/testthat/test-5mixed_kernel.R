context("Fit higher order terms")

test_that("Mixed kernels", {

  mod <- kernL(stack.loss ~ . ^ 2, data = stackloss,
               model = list(kernel = c("FBM", "Canonical", "FBM")))
  classes <- sapply(mod$Hl, function(x) attributes(x)$class)
  expect.classes <- c("FBM,0.5", "Canonical", "FBM,0.5",
                      "FBM,0.5 x Canonical", "FBM,0.5 x FBM,0.5",
                      "Canonical x FBM,0.5")
  expect_equivalent(classes, expect.classes)

})

test_that("Mixed kernel incomplete specification", {

  expect_warning(
    mod <- kernL(stack.loss ~ ., data = stackloss,
               model = list(kernel = c("FBM", "Canonical")))
  )

})
