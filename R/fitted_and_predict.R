fitted.ipriorMod2 <- function(object, intervals = FALSE, alpha = 0.05, ...) {
  y.hat <- object$fitted.values
  res <- list(y = y.hat, resid = object$residuals,
              train.error = object$train.error)
  if (isTRUE(intervals)) {
    res <- c(res, predict_iprior_quantiles(object, NULL, y.hat, alpha))
  }
  class(res) <- "ipriorPredict"
  res
}

predict.ipriorMod2 <- function(object, newdata = list(), y.test = NULL,
                               intervals = FALSE, alpha = 0.05, ...) {
  if (length(newdata) == 0) {
    return(cat("No new data supplied. Use fitted() instead."))
  }
  if (!is.null(object$kernL$formula)) {
    tt <- object$kernL$terms
    Terms <- delete.response(tt)
    xstar <- model.frame(Terms, newdata)
    if (any(colnames(newdata) == object$kernL$yname))
      y.test <- model.extract(model.frame(tt, newdata), "response")
    xrownames <- rownames(xstar)
  } else {
    if (any(sapply(newdata, is.vector))) {
      newdata <- lapply(newdata, as.matrix)
    }
    xstar <- newdata
    xrownames <- rownames(do.call(cbind, newdata))
  }

  Hlam.new <- get_Hlam(object$kernL, object$theta, xstar = xstar)
  res <- predict_iprior(y.test, Hlam.new, object$w, object$intercept)
  names(res$y) <- xrownames
  names(res)[grep("train.error", names(res))] <- "test.error"
  if (isTRUE(intervals)) {
    res <- c(res, predict_iprior_quantiles(object, Hlam.new, res$y, alpha))
  }
  class(res) <- "ipriorPredict"
  res
}

print.ipriorPredict <- function(x, ...) {
  if (!is.null(x$train.error)) {
    cat("Training MSE:", x$train.error, "\n")
  } else if (!is.nan(x$test.error)) {
    cat("Test MSE:", x$test.error, "\n")
  } else {
    cat("Test data not provided.\n")
  }
  cat("\n")
  cat("Predicted values:\n")
  if (is.null(x$lower)) {
    print(x$y)
  } else {
    tab <- data.frame(x$lower, x$y, x$upper)
    lower <- paste0(x$alpha / 2 * 100, "%")
    upper <- paste0((1 - x$alpha / 2) * 100, "%")
    names(tab) <- c(lower, "Mean", upper)
    print(tab)
  }
}

predict_iprior <- function(y, Hlam, w, intercept) {
  # Args: y (data or test data for calculation of errors); Hlam the kernel
  # matrix; w is the posterior mean of I-prior random effects; and the
  # intercept.
  #
  # Output: A list containing the predicted values, residuals and MSE.
  #
  # Notes: This is the main helper function to calculate fitted or predicted
  # values. It appears in iprior2(), fitted() and predict() for ipriorMod2
  # objects.
  y.hat <- Hlam %*% w
  if (!is.null(y)) {
    resid <- y - y.hat
    train.error <- mean(resid ^ 2)
  } else {
    resid <- train.error <- NA
  }
  list(y = as.numeric(y.hat + intercept), resid = as.numeric(resid),
       train.error = train.error)
}

se_yhat <- function(Hlam, Hlam.new, psi) {
  # Args: The Hlam matrix; another Hlam matrix which may be equal to Hlam; and
  # psi the error precision.
  #
  # Output: Standard errors for predicted values y.hat (based on Hlam.new).
  #
  # Notes: Helper function used in predict_iprior_quantiles().
  list2env(eigen_Hlam(Hlam), environment())
  z <- psi * u ^ 2 + 1 / psi
  Vy.inv.Hlam <- vy_inv_a(1 / z, V, t(Hlam.new))
  sqrt(diag(Hlam.new %*% Vy.inv.Hlam) + 1 / psi)
}

predict_iprior_quantiles <- function(object, Hlam.new = NULL, y.hat, alpha) {
  # Args: an ipriorMod2 object; optional Hlam.new if calculating quantiles for
  # new predictions; y.hat are the predicted values for which the quantiles are
  # to be generated; alpha is the significance level.
  #
  # Output: a list containing the lower and upper intervals and the significance
  # level.
  #
  # Notes: Helper function used in fitted and predict for ipriorMod2 objects.
  Hlam <- get_Hlam(object$kernL, object$theta)
  if (is.null(Hlam.new)) Hlam.new <- Hlam
  se <- se_yhat(Hlam, Hlam.new, theta_to_psi(object$theta, object$kernL))
  names(se) <- NULL
  lower <- y.hat + qnorm(alpha / 2) * se
  upper <- y.hat + qnorm(1 - alpha / 2) * se
  list(lower = lower, upper = upper, alpha = alpha)
}
