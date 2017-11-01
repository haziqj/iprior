context("Block B stuff")

test_that("BlockB2()", {

  H1 <- kern_canonical(1:3)
  H2 <- kern_canonical(4:6)
  n <- 3
  intr <- matrix(c(1, 2), ncol = 1)
  expect_equal(BlockB_fn(list(H1), intr = NULL, n = n, p = 1)$BB.msg,
               "Single lambda")
  expect_equal(BlockB_fn(list(H1, H2), intr = NULL, n = n, p = 2)$BB.msg,
               "Multiple lambda with no interactions")
  expect_equal(BlockB_fn(list(H1, H2), intr = intr, n = n, p = 2)$BB.msg,
               "Multiple lambda with parsimonious interactions")

})

test_that("BlockB creation conditions", {

  y <- 1:3
  x1 <- 1:3
  x2 <- 4:6
  expect_true(!is.null(kernL(y, x1, kernel = "se")$BlockBStuff))  # single lambda
  expect_true(!is.null(kernL(y, x1, x2, kernel = "fbm")$BlockBStuff))  # mult lambda
  expect_true(!is.null(kernL(y, x1, x2, interactions = "1:2")$BlockBStuff))
  expect_true(
    is.null(kernL(y, x1, x2, kernel = c("poly2", "fbm"))$BlockBStuff)
  )  # polynomial kernel
  expect_true(
    is.null(kernL(y, x1, x2, kernel = "fbm", est.hurst = TRUE)$BlockBStuff)
  )  # est hurst
  expect_true(
    is.null(kernL(y, x1, x2, kernel = "se", est.lengthscale = TRUE)$BlockBStuff)
  )  # est lengthscale
  expect_true(
    is.null(kernL(y, x1, x2, est.lambda = FALSE, est.psi = FALSE)$BlockBStuff)
  )  # no estimate lambda and psi
  expect_true(
    !is.null(kernL(y, x1, x2, est.lambda = TRUE, est.psi = FALSE)$BlockBStuff)
  )  # no estimate lambda
  expect_true(
    !is.null(kernL(y, x1, x2, est.lambda = FALSE, est.psi = TRUE)$BlockBStuff)
  )  # no estimate psi

})


context("Kernel loader names")
#
# test_that("Names OK using non-formula (no interactions)",{
#
#   # All names entered in call
#   mod <- kernL(y = stackloss$stack.loss, Air.Flow = stackloss$Air.Flow,
#                Water.Temp = stackloss$Water.Temp,
#                Acid.Conc. = stackloss$Acid.Conc.)
#   expect_equal(names(mod$Hl), names(stackloss)[-4])
#
#   # Some names entered in call
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                Water.Temp = stackloss$Water.Temp, stackloss$Acid.Conc.)
#   xname <- names(stackloss)[-4]
#   xname[1] <- "stackloss$Air.Flow"
#   xname[3] <- "stackloss$Acid.Conc."
#   expect_equal(names(mod$Hl), xname)
#
#   # No names entered in call
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                stackloss$Water.Temp, stackloss$Acid.Conc.)
#   xname <- c("stackloss$Air.Flow", "stackloss$Water.Temp",
#              "stackloss$Acid.Conc.")
#   expect_equal(names(mod$Hl), xname)
#
#   # All names through model option
#   xname <- names(stackloss)[-4]
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                stackloss$Water.Temp, stackloss$Acid.Conc.,
#                model = list(xname = xname))
#   expect_equal(names(mod$Hl), xname)
#
#   # Some names through model option
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                stackloss$Water.Temp, stackloss$Acid.Conc.,
#                model = list(xname = "air"))
#   expect_equal(names(mod$Hl), c("air", "stackloss$Water.Temp",
#                                    "stackloss$Acid.Conc."))
#
# })
#
# test_that("Names OK using formula (no interactions)",{
#
#   # All names entered in call
#   mod <- kernL(stack.loss ~ ., data = stackloss)
#   expect_equal(names(mod$Hl), names(stackloss)[-4])
#
#   # All names through model option will not change names
#   xname <- c("air", "water", "acid")
#   mod <- kernL(stack.loss ~ ., data = stackloss, model = list(xname = xname))
#   expect_equal(names(mod$Hl), names(stackloss)[-4])
#
#   # Some names through model option will not change names
#   xname <- c("air")
#   mod <- kernL(stack.loss ~ ., data = stackloss, model = list(xname = xname))
#   expect_equal(names(mod$Hl), names(stackloss)[-4])
#
# })
#
# test_that("Names OK using non-formula with interactions",{
#
#   # All names entered in call
#   mod <- kernL(y = stackloss$stack.loss, air = stackloss$Air.Flow,
#                water = stackloss$Water.Temp, model = list(interactions = "1:2"))
#   expect_equal(names(mod$Hl), c("air", "water", "air:water"))
#
#   # Some names entered in call
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                water = stackloss$Water.Temp, stackloss$Acid.Conc.,
#                model = list(interactions = "1:2"))
#   xname <- c("stackloss$Air.Flow", "water", "stackloss$Acid.Conc.",
#              "stackloss$Air.Flow:water")
#   expect_equal(names(mod$Hl), xname)
#
#   # No names entered in call
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                stackloss$Water.Temp, stackloss$Acid.Conc.,
#                model = list(interactions = "1:2"))
#   xname <- c("stackloss$Air.Flow", "stackloss$Water.Temp",
#              "stackloss$Acid.Conc.", "stackloss$Air.Flow:stackloss$Water.Temp")
#   expect_equal(names(mod$Hl), xname)
#
#   # All names through model option
#   xname <- c("air", "water", "acid", "air:water")
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                stackloss$Water.Temp, stackloss$Acid.Conc.,
#                model = list(interactions = "1:2",
#                             xname = xname[1:3]))
#   expect_equal(names(mod$Hl), xname)
#
#   # Some names through model option
#   xname <- c("air", "stackloss$Water.Temp", "stackloss$Acid.Conc.",
#              "air:stackloss$Water.Temp")
#   mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
#                stackloss$Water.Temp, stackloss$Acid.Conc.,
#                model = list(interactions = "1:2",
#                             xname = xname[1]))
#   expect_equal(names(mod$Hl), xname)
#
# })
#
# test_that("Names OK using formula with interactions",{
#
#   # All names entered in call
#   xname <- c("Air.Flow", "Water.Temp", "Air.Flow:Water.Temp")
#   mod <- kernL(stack.loss ~ Air.Flow * Water.Temp, data = stackloss)
#   expect_equal(names(mod$Hl), xname)
#
#   # All names through model option will not change names
#   xname <- c("air", "water", "acid")
#   mod <- kernL(stack.loss ~ Air.Flow * Water.Temp, data = stackloss,
#                model = list(xname = xname))
#   expect_equal(names(mod$Hl), c("Air.Flow", "Water.Temp",
#                                    "Air.Flow:Water.Temp"))
#
#   # Some names through model option will not change names
#   xname <- c("air")
#   mod <- kernL(stack.loss ~ Air.Flow * Water.Temp, data = stackloss,
#                model = list(xname = xname))
#   expect_equal(names(mod$Hl), c("Air.Flow", "Water.Temp",
#                                    "Air.Flow:Water.Temp"))
#
# })
#
# test_that("Names OK using formula and one.lam = TRUE",{
#
#   xname <- "Air.Flow + Water.Temp + Acid.Conc."
#   mod <- kernL(stack.loss ~ ., data = stackloss, model = list(one.lam = TRUE))
#   expect_equal(names(mod$Hl), xname)
#
# })
