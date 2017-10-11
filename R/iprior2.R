iprior2 <- function(...) UseMethod("iprior2")

iprior2.default <- function(y, ..., kernel = "linear", method = "direct",
                            control = list()) {
  if (is.ipriorKernel2(y)) {
    mod <- y
  } else {
    mod <- kernL2(y = y, ..., kernel = kernel)
  }

  method <- match.arg(method, c("direct", "em", "fixed", "canonical"))

  control_ <- list(
    maxit     = 500,
    stop.crit = 1e-8, # sqrt(.Machine$double.eps),  # roughly 1e-8
    theta0    = NULL,
    silent    = FALSE,
    report    = 10,
    psi.reg   = FALSE  # option for iprior_em_reg()
  )
  control.names <- names(control_)
  control_[(control.names <- names(control))] <- control
  control <- control_

  control.optim <- list(
    fnscale = -2,
    trace = ifelse(isTRUE(control$silent), 0, 1),
    maxit = control$maxit,
    REPORT = control$report
  )

  if (is.null(control$theta0)) {
    theta0 <- rnorm(mod$thetal$n.theta)  # rep(0, mod$thetal$n.theta)
  } else {
    if (length(theta.start) != mod$nt) {
      stop(paste("Incorrect number of parameters specified. Should be", nt))
    }
  }

  # Send to correct estimation method ------------------------------------------
  est.method <- iprior_method_checker(mod, method)
  if (est.method["fixed"]) {
    res <- iprior_fixed(mod)
    res$est.method <- "I-prior fixed."
    res$est.conv <- ""
  } else if (est.method["canonical"]) {
    res <- iprior_canonical(mod, theta0, control.optim)
    res$est.method <- "Direct minimisation of marginal deviance."
  } else {
    if (est.method["em.closed"]) {
      res <- iprior_em_closed(mod, control$maxit, control$stop.crit,
                              control$silent)
      res$est.method <- "Closed-form EM algorithm."
    }
    if (est.method["em.reg"]) {
      res <- iprior_em_reg(mod, control$maxit, control$stop.crit,
                           control$silent, theta0)
      res$est.method <- "Regular EM algorithm."
    }
    if (est.method["direct"]) {
      res <- iprior_direct(mod, loglik_iprior, theta0, control.optim)
      res$est.method <- "Direct minimisation of marginal deviance."
    }
    if (res$conv == 0)
      res$est.conv <- paste0("Converged to within ", control$stop.crit,
                             " tolerance.")
    else
      res$est.conv <- "Convergence criterion not met."
  }

  res$intercept <- attr(mod$y, "scaled:center")
  res$coefficients <- reduce_theta(res$param.full, mod$estl)$theta.reduced
  tmp <- predict_iprior(mod$y, get_Hlam(mod, res$theta), res$w, res$intercept)
  res$fitted.values <- tmp$y
  names(res$fitted.values) <- attr(mod$y, "dimnames")[[1]]
  res$residuals <- tmp$resid
  res$train.error <- tmp$train.error
  res$ipriorKernel <- mod
  res$maxit <- control$maxit
  res$stop.crit <- control$stop.crit

  cl <- match.call()
  res$fullcall <- cl
  cl[[1L]] <- as.name("iprior2")
  # names(cl)[2] <- ""  # get rid of "y ="
  names(cl)[-(1:2)] <- paste0("X", seq_along(names(cl)[-(1:2)]))
  res$call <- cl

  class(res) <- "ipriorMod2"
  res
}

iprior2.formula <- function(formula, data, kernel = "linear", method = "direct",
                            control = list(), ...) {
  mod <- kernL2.formula(formula, data, kernel = kernel, ...)
  res <- iprior2.default(y = mod, method = method, control = control)

  cl <- match.call()
  res$fullcall <- cl
  cl[[1L]] <- as.name("iprior2")
  res$call <- cl
  res
}

print.ipriorMod2 <- function(x, digits = 5) {
  loglik.max <- x$loglik[length(x$loglik)]
  cat("Log-likelihood value:", loglik.max, "\n")
  cat("\n")
  if (x$ipriorKernel$thetal$n.theta > 0)
    print(round(coef(x), digits))
  else
    cat("No hyperparameters estimated.")
}

logLik.ipriorMod2 <- function(object, ...) {
  res <- object$loglik
  res[length(res)]
}

summary.ipriorMod2 <- function(object) {
  resid.summ <- round(summary(residuals(object))[-4], 4)

  # need to use delta method here!
  coef <- object$param.full
  se <- expand_theta(object$se, object$ipriorKernel$thetal$theta.drop, NA)
  zval <- coef / se
  tab <- cbind(Estimate   = round(coef, 4),
               S.E.       = round(se, 4),
               z          = round(zval, 3),
               `P[|Z>z|]` = round(2 * pnorm(-abs(zval)), 3))

  # rename rownames, remove psi

  param.tab <- theta_to_param(object$theta, object$ipriorKernel)
  kernels.used <- rep(NA, nrow(param.tab))
  for (i in seq_along(param.tab$kernels)) {
    kernels.used[i] <- kernel_summary_translator(param.tab$kernels[i])
  }
  x.kern <- unique.kernels <- unique(kernels.used)
  for (i in seq_along(unique.kernels)) {
    ind <- kernels.used %in% unique.kernels[i]
    xs <- paste0(object$kernL$xname[ind], collapse = ", ")
    x.kern[i] <- paste0(unique.kernels[i], " (", xs, ")\n")
  }

  res <- list(resid.summ = resid.summ, tab = tab, loglik = logLik(object),
              error = object$train.error, call = object$call, x.kern = x.kern,
              est.method = object$est.method, est.conv = object$est.conv,
              niter = object$niter, maxit = object$maxit, time = object$time)
  class(res) <- "ipriorMod2_summary"
  res
}

kernel_summary_translator <- function(x) {
  # Notes: Not vectorised.
  if (is.kern_linear(x)) res <- "Linear"
  if (is.kern_pearson(x)) res <- "Pearson"
  else {
    hyperparam <-  signif(get_hyperparam(x), 3)
    if (is.kern_fbm(x)) {
      res <- paste0("Fractional Brownian motion with Hurst ", hyperparam)
    }
    if (is.kern_se(x)) {
      res <- paste0("Squared exponential with lengthscale ", hyperparam)
    }
    if (is.kern_poly(x)) {
      degree <- get_polydegree(x)
      res <- paste0("Polynomial degree ", degree, " with offset ", hyperparam)
    }
  }
  res
}

print.ipriorMod2_summary <- function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("RKHS used:\n")
  cat(x$x.kern)
  cat("\n")
  cat("Residuals:\n")
  print(x$resid.summ)
  cat("\n")
  cat("Hyperparameters:\n")
  tmp <- capture.output(printCoefmat(x$tab, P.values = TRUE, has.Pvalue = TRUE))
  cat(paste(gsub("NA", "  ", tmp), collapse = "\n"))
  cat("\n\n")
  cat(x$est.method)
  cat(" Iterations:", paste0(x$niter, "/", x$maxit), "\n")
  cat(x$est.conv)
  cat(" Time taken: ")
  print(x$time)
  cat("\n")
  cat("Log-likelihood value:", x$loglik, "\n")
  cat("Training mean squared error:", x$error, "\n")
  # cat("Standard deviation of errors: xxx with S.E.: xxx\n")
}

iprior_method_checker <- function(object, method) {
  res <- rep(FALSE, 5)
  names(res) <- c("fixed", "canonical", "em.closed", "em.reg", "direct")

  if (object$thetal$n.theta == 0 | method == "fixed") {
    res["fixed"] <- TRUE
    # if (method != "fixed") warning("No hyperparameters estimated. Using fixed estimation method.", call. = FALSE)
  } else if (method == "canonical") {
    if (all(is.kern_linear(object$kernels))) {
      res["canonical"] <- TRUE
    } else {
      res["direct"] <- TRUE
      warning("Non-linear kernels used. Using direct estimation method.", call. = FALSE)
    }
    if (object$no.int > 0) {
      res["direct"] <- TRUE
      warning("Canonical method not possible with interactions. Using direct estimation method.", call. = FALSE)
    }
  } else if (method == "em") {
    if (is.null(object$BlockBStuff)) {
      res["em.reg"] <- TRUE
    } else {
      res["em.closed"] <- TRUE
    }
  } else {
    res["direct"] <- TRUE
  }
  res
}
