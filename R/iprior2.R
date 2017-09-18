iprior2 <- function(y, ..., kernel = "linear", control = list()) {
  mod <- kernL2(y = y, ..., kernel = kernel)

  control_ <- list(
    maxit     = 100,
    stop.crit = sqrt(.Machine$double.eps),  # roughly 1e-8
    theta0    = NULL,
    silent    = FALSE,
    report    = 10
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
    theta0 <- rnorm(mod$nt) #rep(0, mod$nt)
  } else {
    if (length(theta.start) != mod$nt) {
      stop(paste("Incorrect number of parameters specified. Should be", nt))
    }
  }

  if (mod$nt == 0) {
    res <- iprior_fixed(mod)
  } else {
    res <- iprior_direct(mod, loglik_iprior, theta0, control.optim)
  }
  res$coefficients <- reduce_theta(res$param.full, mod$est.list)$theta.reduced
  tmp <- predict_iprior(mod$y, get_Hlam(mod, res$theta), res$w)
  res$fitted.values <- as.numeric(tmp$y.hat)
  res$residuals <- as.numeric(tmp$resid)
  res$train.error <- tmp$train.error
  res$kernL <- mod
  class(res) <- "ipriorMod2"
  res
}

print.ipriorMod2 <- function(x, digits = 5) {
  loglik.max <- x$loglik[length(x$loglik)]
  cat("Log-likelihood value:", loglik.max, "\n")
  cat("\n")
  if (x$kernL$nt > 0)
    print(round(x$theta, digits))
  else
    cat("No hyperparameters estimated.")
}

logLik.ipriorMod2 <- function(object, ...) {
  res <- object$loglik
  res[length(res)]
}

summary.ipriorMod2 <- function(object) {
  resid.summ <- round(summary(residuals(object))[-4], 4)

  se <- object$se  # need to use delta method here
  zval <- coef(object) / se
  tab <- cbind(Estimate   = round(coef(object), 4),
               S.E.       = round(se, 4),
               z          = round(zval, 3),
               `P[|Z>z|]` = round(2 * pnorm(-abs(zval)), 3))
  # rename rownames, remove psi

  res <- list(resid.summ = resid.summ, tab = tab, loglik = logLik(object),
              error = object$train.error)
  class(res) <- "ipriorMod2_summary"
  res
}

print.ipriorMod2_summary <- function(x) {
  cat("Call:\n")
  cat("~~~ the call ~~~\n")
  cat("\n")
  cat("RKHS used:\n")
  cat("~~~ the RKHS ~~~\n")
  cat("\n")
  cat("Residuals:\n")
  print(x$resid.summ)
  cat("\n")
  cat("Hyperparameters:\n")
  printCoefmat(x$tab, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  cat("Convergence and iterations.\n")
  cat("Log-likelihood value:", x$loglik, "\n")
  cat("Training error (residual sum of squares):", x$error, "\n")
  cat("Standard deviation of errors: xxx with S.E.: xxx\n")
  cat("\n")
}


