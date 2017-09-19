context("Kernel loader")

test_that("Kernel loader using non-formula",{

	mod <- kernL(y = stackloss$stack.loss, Air.Flow = stackloss$Air.Flow,
	             Water.Temp = stackloss$Water.Temp,
	             Acid.Conc. = stackloss$Acid.Conc.)
	expect_is(mod, "ipriorKernel")
	expect_equal(mod$p, 3)

})

test_that("Kernel loader using formula",{

  mod <- kernL(stack.loss ~ ., data = stackloss)
  expect_is(mod, "ipriorKernel")
  expect_equal(mod$p, 3)

})

test_that("one.lam = TRUE works properly",{

  mod1 <- kernL(stack.loss ~ ., data = stackloss, model = list(one.lam = TRUE))
  mod2 <- kernL(y = stackloss$stack.loss, x = stackloss[-4])
  expect_equivalent(mod1$Hl, mod2$Hl)

})

test_that("Kernel loader using formula",{

  mod <- kernL(stack.loss ~ ., data = stackloss)
  expect_is(mod, "ipriorKernel")
  expect_equal(mod$p, 3)

})

test_that("Can't use interactions with one.lam in formula input",{

  expect_error(kernL(stack.loss ~ . ^ 2, data = stackloss,
                     model = list(one.lam = TRUE)))

})

test_that("Incorrect specification of interactions",{

  expect_error(kernL(y = stackloss$stack.loss, air = stackloss$Air.Flow,
               water = stackloss$Water.Temp, model = list(interactions = 1:2)))
  expect_error(kernL(y = stackloss$stack.loss, air = stackloss$Air.Flow,
               water = stackloss$Water.Temp, model = list(interactions = "12")))

})

test_that("Overriding Hurst coefficient warning",{

  mod1 <- kernL(stack.loss ~ ., stackloss, model = list(kernel = "FBM,0.1"))
  mod2 <- kernL(stack.loss ~ ., stackloss, model = list(kernel = "FBM",
                                                        Hurst = 0.1))
  expect_equivalent(mod1$Hl, mod2$Hl)
  expect_warning(
    mod3 <- kernL(stack.loss ~ ., stackloss,
                  model = list(kernel = "FBM,0.9", Hurst = 0.1))

  )

})


test_that("Kernel translator",{

  expect_true(is.kern_linear(kernel_translator(1:3, kernel = "linear")))

  res <- kernel_translator(1:3, kernel = "fbm,0.7")
  expect_true(is.kern_fbm(res))
  expect_equal(get_hyperparam(res), 0.7)

  res <- kernel_translator(1:3, kernel = "se,0.7")
  expect_true(is.kern_se(res))
  expect_equal(get_hyperparam(res), 0.7)

  res1 <- kernel_translator(1:3, kernel = "poly,0.7")
  res2 <- kernel_translator(1:3, kernel = "poly3")
  res3 <- kernel_translator(1:3, kernel = "poly4,0.9")
  expect_true(all(is.kern_poly(res1), is.kern_poly(res2), is.kern_poly(res3)))
  expect_equal(get_hyperparam(res1), 0.7)
  expect_equal(get_hyperparam(res2), 0)
  expect_equal(get_hyperparam(res3), 0.9)
  expect_equal(get_polydegree(res1), 2)
  expect_equal(get_polydegree(res2), 3)
  expect_equal(get_polydegree(res3), 4)

})

test_that("Kernel to param to theta", {

  # Remember that kernels vector need to have the hyperparameters specified,
  # i.e. fbm,0.5 instead of just fbm, etc.
  # e.g. no of theta = 5 + 1 + 1 + 1 = 8
  # param table is 5 x 5
  kernels <- c("linear", "fbm,0.7", "se,2", "poly3,0.15", "pearson")
  which.pearson <- c(F, F, F, F, T)
  poly.degree <- c(NA, NA, NA, 3, NA)
  lambda <- rep(1, 5)
  psi <- 1
  names(psi) <- "psi"

  res1 <- kernel_to_param(kernels, lambda)
  res2 <- param_to_theta(res1, list(est.lambda      = c(FALSE, TRUE),
                                    est.hurst       = TRUE,
                                    est.lengthscale = FALSE,
                                    est.offset      = TRUE,
                                    est.psi         = TRUE))
  res3 <- theta_to_param(res2$theta, list(param.na      = res2$param.na,
                                          theta.drop    = res2$theta.drop,
                                          theta.omitted = res2$theta.omitted,
                                          which.pearson = which.pearson,
                                          poly.degree   = poly.degree))

  expect_equal(kernels, res3$kernels)
  expect_equal(length(collapse_param(res3)$param), 8)
  expect_equal(theta_to_psi(res2$theta, list(theta.drop    = res2$theta.drop,
                                             theta.omitted = res2$theta.omitted)),
               psi)
  expect_equal(res1, res3)

})

























