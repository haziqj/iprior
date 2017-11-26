iprior_cv <- function(...) UseMethod("iprior_cv")

iprior_cv.default <- function(y, ..., folds = 2, kernel = "linear",
                              method = "direct", control = list(),
                              interactions = NULL, est.lambda = TRUE,
                              est.hurst = FALSE, est.lengthscale = FALSE,
                              est.offset = FALSE, est.psi = TRUE,
                              fixed.hyp = NULL, lambda = 1, psi = 1,
                              nystrom = FALSE, nys.seed = NULL, cv.par = TRUE) {
  n <- length(y)
  if (folds > n) folds <- n
  if (folds <= 1) {
    warning("Number of folds needs to be 2 or greater. Defaulting to 2.")
    folds <- 2
  }
  samp <- seq_len(n) #sample(seq_len(n))
  n.cv <- rep(n %/% folds, folds)
  n.cv[1] <- n.cv[1] + n %% folds
  n.cv.cum <- c(0, cumsum(n.cv))
  cv.samp <- list(NULL)
  for (i in seq_len(folds)) {
    cv.samp[[i]] <- samp[(n.cv.cum[i] + 1):n.cv.cum[i + 1]]
  }
  control$silent <- FALSE

  snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
  pb <- txtProgressBar(min = 0, max = folds, style = 1)

  if (isTRUE(control$restarts)) cv.par <- FALSE

  if (!isTRUE(cv.par)) {
    res <- matrix(NA, ncol = 3, nrow = folds)
    for (k in seq_along(cv.samp)) {
      mod <- kernL(y = y, ..., test.samp = cv.samp[[k]], kernel = kernel,
                   interactions = interactions, est.lambda = est.lambda,
                   est.hurst = est.hurst, est.lengthscale = est.lengthscale,
                   est.offset = est.offset, est.psi = est.psi,
                   fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
                   nystrom = nystrom, nys.seed = nys.seed)
      mod <- iprior.ipriorKernel(mod, method = method, control = control)
      intercept <- get_intercept(mod)
      xstar <- mod$ipriorKernel$Xl.test
      y.test <- intercept + mod$ipriorKernel$y.test
      Hlam.new <- get_Htildelam(mod$ipriorKernel, mod$theta, xstar)
      y.hat.new <- intercept + Hlam.new %*% mod$w
      tmp <- predict_iprior(y.test, y.hat.new)
      res[k, ] <- c(logLik(mod), get_mse(mod), tmp$train.error)
      setTxtProgressBar(pb, k)
    }
  } else {
    cl <- parallel::makeCluster(parallel::detectCores())
    doSNOW::registerDoSNOW(cl)
    res <- foreach::`%dopar%`(
      foreach::foreach(
        k = seq_along(cv.samp),
        .combine = rbind,
        .packages = "iprior",
        .options.snow = snow.options.list
      ), {
        mod <- kernL(y = y, ..., test.samp = cv.samp[[k]], kernel = kernel)
        mod <- iprior.ipriorKernel(mod, method = method, control = control)
        intercept <- get_intercept(mod)
        xstar <- mod$ipriorKernel$Xl.test
        y.test <- mod$ipriorKernel$y.test
        Hlam.new <- get_Htildelam(mod$ipriorKernel, mod$theta, xstar)
        y.hat.new <- intercept + Hlam.new %*% mod$w
        tmp <- predict_iprior(y.test, y.hat.new)
        c(logLik(mod), get_mse(mod), tmp$train.error)
      }
    )
    if (!isTRUE(control$silent)) close(pb)
    parallel::stopCluster(cl)
  }
  close(pb)

  structure(list(loglik = res[, 1], train.mse = res[, 2], test.mse = res[, 3],
                 folds = folds, n = n),
            class = "iprior_xv")
}

iprior_cv.formula <- function(formula, data, folds = 2, one.lam = FALSE, ...) {
  list2env(formula_to_xy(formula = formula, data = data, one.lam = one.lam),
           envir = environment())
  iprior_cv.default(y, Xl.formula = Xl, interactions = interactions,
                    folds = folds, ...)
}

print.iprior_xv <- function(x, ...) {
  cv.method <- ifelse(x$folds == x$n, "(Leave-one-out Cross Validation)",
                      paste0("(", x$folds, "-fold Cross Validation)"))
  cat("Test MSE =", mean(x$test.mse), cv.method, "\n")
}
