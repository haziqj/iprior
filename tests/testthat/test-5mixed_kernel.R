context("Mixed kernels")

test_that("Mixed kernels", {

  expect_warning(
    mod <- kernL(stack.loss ~ ., data = stackloss,
                 model = list(order = c(1,"1^2",2),
                              kernel = c("FBM", "Canonical", "FBM")))
  )

})

test_that("Mixed kernels multiple Hurst", {

  mod1 <- kernL(stack.loss ~ ., stackloss,
                model = list(kernel = c("FBM,0.1", "Canonical", "FBM,0.1")))
  mod2 <- kernL(stack.loss ~ ., stackloss,
                model = list(kernel = c("FBM", "Canonical", "FBM"), Hurst = 0.1))
  expect_equivalent(mod1$Hl, mod2$Hl)
  expect_warning(
    kernL(stack.loss ~ ., stackloss,
          model = list(kernel = c("FBM,0.9", "Canonical", "FBM,0.9"), Hurst = 0.1))
  )

})

test_that("Warn when using one.lam = TRUE", {

  expect_warning(
    mod <- kernL(stack.loss ~ ., data = stackloss,
                 model = list(one.lam = TRUE,
                              kernel = c("FBM", "Canonical", "FBM")))
  )

})

test_that("Automatic and manual are the same", {

  mod1 <- kernL(len ~ ., data = ToothGrowth,
                model = list(one.lam = TRUE))
  toot <- ToothGrowth
  toot$supp <- as.numeric(ToothGrowth$supp)
  suppressWarnings(mod2 <- kernL(len ~ ., data = toot,
                   model = list(one.lam = TRUE,
                                kernel = c("Pearson", "Canonical"))))
  expect_equivalent(mod1$Hl, mod2$Hl)

})
