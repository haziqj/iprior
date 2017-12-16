# context("LEGACY kernel loader")
#
# test_that("Kernel loader using non-formula",{
#
#   mod <- .kernL(y = stackloss$stack.loss, Air.Flow = stackloss$Air.Flow,
#                Water.Temp = stackloss$Water.Temp,
#                Acid.Conc. = stackloss$Acid.Conc.)
#   expect_is(mod, "ipriorKernel_old")
#   expect_equal(mod$p, 3)
#   tmp <- capture.output(print(mod))
#   expect_equal(tmp, tmp)
#
# })
#
# test_that("Kernel loader using formula",{
#
#   mod <- .kernL(stack.loss ~ . ^ 2, data = stackloss)
#   expect_is(mod, "ipriorKernel_old")
#   expect_equal(mod$p, 3)
#
# })
#
# test_that("one.lam = TRUE works properly",{
#
#   mod1 <- .kernL(stack.loss ~ ., data = stackloss, model = list(one.lam = TRUE))
#   mod2 <- .kernL(y = stackloss$stack.loss, x = stackloss[-4])
#   expect_equivalent(mod1$Hl, mod2$Hl)
#
# })
#
# test_that("Higher order terms",{
#
#   mod <- .kernL(stack.loss ~ . ^ 3, data = stackloss)
#   expect_is(mod, "ipriorKernel_old")
#   expect_equal(mod$p, 3)
#
# })
#
# test_that("Can't use interactions with one.lam in formula input",{
#
#   expect_error(.kernL(stack.loss ~ . ^ 2, data = stackloss,
#                      model = list(one.lam = TRUE)))
#
# })
#
# test_that("Incorrect specification of interactions",{
#
#   expect_error(.kernL(y = stackloss$stack.loss, air = stackloss$Air.Flow,
#                      water = stackloss$Water.Temp, model = list(interactions = 1:2)))
#   expect_error(.kernL(y = stackloss$stack.loss, air = stackloss$Air.Flow,
#                      water = stackloss$Water.Temp, model = list(interactions = "12")))
#
# })
#
# test_that("Overriding Hurst coefficient warning",{
#
#   mod1 <- .kernL(stack.loss ~ ., stackloss, model = list(kernel = "FBM,0.1"))
#   mod2 <- .kernL(stack.loss ~ ., stackloss, model = list(kernel = "FBM",
#                                                         Hurst = 0.1))
#   expect_equivalent(mod1$Hl, mod2$Hl)
#   expect_warning(
#     mod3 <- .kernL(stack.loss ~ ., stackloss,
#                   model = list(kernel = "FBM,0.9", Hurst = 0.1))
#
#   )
#
# })
