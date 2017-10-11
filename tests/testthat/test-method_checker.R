context("Estimation method checker")

test_that("Fixed", {

  mod <- list(thetal = list(n.theta = 0))
  res <- iprior_method_checker(mod, "fixed")
  expect_true(res["fixed"])

  mod <- list(thetal = list(n.theta = 10))
  res <- iprior_method_checker(mod, "fixed")
  expect_true(res["fixed"])

  # mod <- list(thetal = list(n.theta = 0))
  # expect_warning(res <- iprior_method_checker(mod, "em"))
  # expect_true(res["fixed"])

})

test_that("Canonical", {

  mod <- list(thetal = list(n.theta = 10), kernels = "linear", no.int = 0)
  res <- iprior_method_checker(mod, "canonical")
  expect_true(res["canonical"])

  mod <- list(thetal = list(n.theta = 10), kernels = "se,1", no.int = 0)
  expect_warning(iprior_method_checker(mod, "canonical"))
  expect_true(suppressWarnings(
    iprior_method_checker(mod, "canonical")["direct"]
  ))

  mod <- list(thetal = list(n.theta = 10), kernels = "linear", no.int = 1)
  expect_warning(iprior_method_checker(mod, "canonical"))
  expect_true(suppressWarnings(
    iprior_method_checker(mod, "canonical")["direct"]
  ))

})

test_that("EM", {

  mod <- list(thetal = list(n.theta = 10), BlockBStuff = 1)
  res <- iprior_method_checker(mod, "em")
  expect_true(res["em.closed"])

  mod <- list(thetal = list(n.theta = 10), BlockBStuff = NULL)
  res <- iprior_method_checker(mod, "em")
  expect_true(res["em.reg"])

})

test_that("Direct", {

  mod <- list(thetal = list(n.theta = 10))
  res <- iprior_method_checker(mod, "direct")
  expect_true(res["direct"])

})
