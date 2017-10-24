context("Utilities")

test_that("eigenCpp", {

  mat <- kern_linear(1:3)
  mat <- mat %*% mat
  res1 <- eigenCpp(mat)
  res2 <- eigen(mat)
  expect_equal(res1$values, sort(res2$values), tolerance = 1e-6)

})

test_that("fastSquare", {

  mat <- kern_linear(1:3)
  res1 <- mat %*% mat
  res2 <- fastSquare(mat)
  expect_equal(res1, res2, tolerance = 1e-6)

})

test_that("fastSquareRoot", {

  mat <- kern_linear(1:3)
  mat2 <- mat %*% mat
  res1 <- fastSquareRoot(mat2)
  res2 <- fastSquareRoot2(mat2)
  tmp <- eigen(mat2)
  res3 <- tmp$vec %*% (t(tmp$vec) * sqrt(tmp$val))
  expect_equal(res1, res2, tolerance = 1e-6)
  expect_equal(res1, res3, tolerance = 1e-6)
  expect_equal(res2, res3, tolerance = 1e-6)

})

test_that("fastSquareRoot", {

  mat <- kern_linear(1:3)
  mat2 <- mat %*% mat
  res1 <- fastSquareRoot(mat2)
  res2 <- fastSquareRoot2(mat2)
  tmp <- eigen(mat2)
  res3 <- tmp$vec %*% (t(tmp$vec) * sqrt(tmp$val))
  expect_equal(res1, res2, tolerance = 1e-6)
  expect_equal(res1, res3, tolerance = 1e-6)
  expect_equal(res2, res3, tolerance = 1e-6)

})

test_that("fastVDiag", {

  # for some reason it overwrites to environment...
  mat <- kern_linear(1:3)
  mat <- mat %*% mat
  tmp <- eigenCpp(mat)
  res <- fastVDiag(tmp$vec, tmp$val)
  expect_equal(mat, res, tolerance = 1e-6)

})
