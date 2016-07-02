context("Kernel loader names")

test_that("Names OK using non-formula (no interactions)",{

  # All names entered in call
  mod <- kernL(y = stackloss$stack.loss, Air.Flow = stackloss$Air.Flow,
               Water.Temp = stackloss$Water.Temp,
               Acid.Conc. = stackloss$Acid.Conc.)
  expect_equal(names(mod$Hl), names(stackloss)[-4])

  # Some names entered in call
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               Water.Temp = stackloss$Water.Temp, stackloss$Acid.Conc.)
  xname <- names(stackloss)[-4]
  xname[1] <- "stackloss$Air.Flow"
  xname[3] <- "stackloss$Acid.Conc."
  expect_equal(names(mod$Hl), xname)

  # No names entered in call
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               stackloss$Water.Temp, stackloss$Acid.Conc.)
  xname <- c("stackloss$Air.Flow", "stackloss$Water.Temp",
             "stackloss$Acid.Conc.")
  expect_equal(names(mod$Hl), xname)

  # All names through model option
  xname <- names(stackloss)[-4]
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               stackloss$Water.Temp, stackloss$Acid.Conc.,
               model = list(xname = xname))
  expect_equal(names(mod$Hl), xname)

  # Some names through model option
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               stackloss$Water.Temp, stackloss$Acid.Conc.,
               model = list(xname = "air"))
  expect_equal(names(mod$Hl), c("air", "stackloss$Water.Temp",
                                   "stackloss$Acid.Conc."))

})

test_that("Names OK using formula (no interactions)",{

  # All names entered in call
  mod <- kernL(stack.loss ~ ., data = stackloss)
  expect_equal(names(mod$Hl), names(stackloss)[-4])

  # All names through model option will not change names
  xname <- c("air", "water", "acid")
  mod <- kernL(stack.loss ~ ., data = stackloss, model = list(xname = xname))
  expect_equal(names(mod$Hl), names(stackloss)[-4])

  # Some names through model option will not change names
  xname <- c("air")
  mod <- kernL(stack.loss ~ ., data = stackloss, model = list(xname = xname))
  expect_equal(names(mod$Hl), names(stackloss)[-4])

})

test_that("Names OK using non-formula with interactions",{

  # All names entered in call
  mod <- kernL(y = stackloss$stack.loss, air = stackloss$Air.Flow,
               water = stackloss$Water.Temp, model = list(interactions = "1:2"))
  expect_equal(names(mod$Hl), c("air", "water", "air:water"))

  # Some names entered in call
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               water = stackloss$Water.Temp, stackloss$Acid.Conc.,
               model = list(interactions = "1:2"))
  xname <- c("stackloss$Air.Flow", "water", "stackloss$Acid.Conc.",
             "stackloss$Air.Flow:water")
  expect_equal(names(mod$Hl), xname)

  # No names entered in call
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               stackloss$Water.Temp, stackloss$Acid.Conc.,
               model = list(interactions = "1:2"))
  xname <- c("stackloss$Air.Flow", "stackloss$Water.Temp",
             "stackloss$Acid.Conc.", "stackloss$Air.Flow:stackloss$Water.Temp")
  expect_equal(names(mod$Hl), xname)

  # All names through model option
  xname <- c("air", "water", "acid", "air:water")
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               stackloss$Water.Temp, stackloss$Acid.Conc.,
               model = list(interactions = "1:2",
                            xname = xname[1:3]))
  expect_equal(names(mod$Hl), xname)

  # Some names through model option
  xname <- c("air", "stackloss$Water.Temp", "stackloss$Acid.Conc.",
             "air:stackloss$Water.Temp")
  mod <- kernL(y = stackloss$stack.loss, stackloss$Air.Flow,
               stackloss$Water.Temp, stackloss$Acid.Conc.,
               model = list(interactions = "1:2",
                            xname = xname[1]))
  expect_equal(names(mod$Hl), xname)

})

test_that("Names OK using formula with interactions",{

  # All names entered in call
  xname <- c("Air.Flow", "Water.Temp", "Air.Flow:Water.Temp")
  mod <- kernL(stack.loss ~ Air.Flow * Water.Temp, data = stackloss)
  expect_equal(names(mod$Hl), xname)

  # All names through model option will not change names
  xname <- c("air", "water", "acid")
  mod <- kernL(stack.loss ~ Air.Flow * Water.Temp, data = stackloss,
               model = list(xname = xname))
  expect_equal(names(mod$Hl), c("Air.Flow", "Water.Temp",
                                   "Air.Flow:Water.Temp"))

  # Some names through model option will not change names
  xname <- c("air")
  mod <- kernL(stack.loss ~ Air.Flow * Water.Temp, data = stackloss,
               model = list(xname = xname))
  expect_equal(names(mod$Hl), c("Air.Flow", "Water.Temp",
                                   "Air.Flow:Water.Temp"))

})

test_that("Names OK using formula and one.lam = TRUE",{

  xname <- "Air.Flow + Water.Temp + Acid.Conc."
  mod <- kernL(stack.loss ~ ., data = stackloss, model = list(one.lam = TRUE))
  expect_equal(names(mod$Hl), xname)

})
