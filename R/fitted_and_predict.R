################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2018  Haziq Jamil
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

#' Obtain predicted values from \code{ipriorMod} objects
#'
#' @param object,x An \code{ipriorMod} object.
#' @param intervals Logical. Calculate the credibility intervals for the fitted
#'   values. Defaults to \code{FALSE}.
#' @param alpha The significance level for the credibility intervals. This is a
#'   number between 0 and 1.
#' @param newdata Either a data frame when using formula method, or a list of
#'   vectors/matrices if using default method. Either way, the new data must be
#'   structurally similar to the original data used to fit the model.
#' @param y.test (Optional) Test data, in order to compute test error rates.
#' @param rows (Optional) The number of values/rows to display.
#' @param dp (Optional) Decimal places for the values.
#' @param ... Not used.
#'
#' @return A list of class \code{ipriorPredict} containing the fitted values,
#'   residuals (observed minus fitted), the training mean squared error, and the
#'   lower and upper intervals (if called).
#'
#' @examples
#' dat <- gen_smooth(20)
#' mod <- iprior(y ~ ., dat, kernel = "se")
#' fitted(mod)
#' fitted(mod, intervals = TRUE)
#' predict(mod, gen_smooth(5))
#'
#' with(dat, mod <<- iprior(y, X, kernel = "poly"))
#' newdat <- gen_smooth(30)
#' mod.pred <- predict(mod, list(newdat$X), y.test = newdat$y, intervals = TRUE)
#' str(mod.pred)
#' print(mod.pred, row = 5)
#'
#' @name predict
NULL

#' @rdname predict
#' @export
fitted.ipriorMod <- function(object, intervals = FALSE, alpha = 0.05, ...) {
  y.hat <- object$fitted.values
  res <- list(y = y.hat, resid = object$residuals,
              train.error = object$train.error)
  if (isTRUE(intervals)) {
    res <- c(res, predict_iprior_quantiles(object, NULL, y.hat, alpha))
  }
  class(res) <- "ipriorPredict"
  res
}

#' @rdname predict
#' @export
predict.ipriorMod <- function(object, newdata = list(), y.test = NULL,
                               intervals = FALSE, alpha = 0.05, ...) {
  if (length(newdata) == 0) {
    return(cat("No new data supplied. Use fitted() instead."))
  }
  list2env(before_predict(object, newdata), envir = environment())

  res <- predict_iprior(object, xstar, y.test, intervals = intervals,
                        alpha = alpha)
  names(res$y) <- xrownames
  class(res) <- "ipriorPredict"
  res
}

before_predict <- function(object, newdata) {
  # A helper function to 1) test whether mod is fitted using formula or not,
  # then 2) return xstar, y.test (if data frame supplied) and xrownames. I
  # separated this from the actual predict function because it is useful to say,
  # calculate the kernel matrix at new input points in get_kern_matrix().
  #
  # Args: ipriorMod object, and newdata list or data frame.
  #
  # Returns: A list containing y.test (if data frame supplied), xstar and
  # xrownames.
  if (!is.null(object$ipriorKernel$formula)) {
    tmp <- terms_to_xy(object$ipriorKernel, newdata)
    y.test <- tmp$y
    xstar <- tmp$Xl
    xrownames <- rownames(newdata)
    return(list(y.test = y.test, xstar = xstar, xrownames = xrownames))
  } else {
    if (any(sapply(newdata, is.vector))) {
      newdata <- lapply(newdata, as.matrix)
    }
    xstar <- newdata
    xrownames <- rownames(do.call(cbind, newdata))
    return(list(xstar = xstar, xrownames = xrownames))
  }
}

#' @rdname predict
#' @export
print.ipriorPredict <- function(x, rows = 10, dp = 3, ...) {
  if (!is.null(x$train.error)) {
    cat("Training RMSE:", sqrt(x$train.error), "\n")
  } else if (!is.nan(x$test.error)) {
    cat("Test RMSE:", sqrt(x$test.error), "\n")
  } else {
    cat("Test data not provided.\n")
  }
  cat("\n")
  cat("Predicted values:\n")

  rows <- min(length(x$y), rows)
  if (is.null(x$lower)) {
    y.print <- round(x$y[seq_len(rows)], dp)
    print(y.print)
  } else {
    tab <- data.frame(x$lower, x$y, x$upper)
    tab <- dec_plac(tab, dp)
    lower <- paste0(x$alpha / 2 * 100, "%")
    upper <- paste0((1 - x$alpha / 2) * 100, "%")
    names(tab) <- c(lower, "Mean", upper)
    print(tab[seq_len(rows), ])
  }
  if (length(x$y) > rows) cat("# ... with", length(x$y) - rows, "more values\n")
}

predict_iprior <- function(object, xstar, y.test, intervals = FALSE,
                           alpha = 0.05) {
  # This is the main helper function to calculate predicted values, and if
  # required, the credibility intervals for prediction.
  #
  # Args: An ipriorMod object, a list of new data points (in the style of
  # non-formula, even if fitted using formula), and data or test data for
  # calculation of residuals and prediction error.
  #
  # Returns: A list containing the predicted values, residuals and MSE.
  Hlam.new <- get_Htildelam(object$ipriorKernel, object$theta, xstar)
  y.hat.new <- object$ipriorKernel$intercept + Hlam.new %*% object$w
  res <- mse_iprior(y.test, y.hat.new)

  if (isTRUE(intervals)) {
    res <- c(res, predict_iprior_quantiles(object, Hlam.new, res$y, alpha))
  }

  names(res)[grep("train.error", names(res))] <- "test.error"
  res
}

mse_iprior <- function(y, y.hat) {
  # A helper function to calculate residuals and mean squared error from y and
  # y.hat. It appears in iprior(), fitted() and predict() for ipriorMod objects.
  #
  # Args: y (data or test data for calculation of errors); y.hat is the
  # predicted y values.
  #
  # Returns: A list containing the predicted values, residuals and MSE.
  y.hat <- as.numeric(y.hat)
  if (!is.null(y)) {
    resid <- y - y.hat
    train.error <- mean(resid ^ 2)
  } else {
    resid <- train.error <- NA
  }
  list(y = y.hat, resid = as.numeric(resid), train.error = train.error)
}

se_yhat <- function(Hlam, Hlam.new, psi, nystrom = FALSE) {
  # Obtains standard errors for predicted values y.hat (based on Hlam.new). This
  # is a helper function used in predict_iprior_quantiles().
  #
  # Args: The Hlam matrix; another Hlam matrix which may be equal to Hlam; and
  # psi the error precision.
  #
  # Returns: standard errors for predicted values.
  if (isTRUE(nystrom)) {
    eigen_Hlam_nys(Hlam, environment())  # assign u and V to environment
  } else {
    eigen_Hlam(Hlam, environment())  # assign u and V to environment
  }
  z <- psi * u ^ 2 + 1 / psi
  Vy.inv.Hlam <- A_times_a(1 / z, V, t(Hlam.new))
  sqrt(diag(Hlam.new %*% Vy.inv.Hlam) + 1 / psi)
}

predict_iprior_quantiles <- function(object, Hlam.new = NULL, y.hat, alpha) {
  # Helper function used in fitted and predict for ipriorMod objects.
  #
  # Args: An ipriorMod object; optional Hlam.new if calculating quantiles for
  # new predictions; y.hat are the predicted values for which the quantiles are
  # to be generated; alpha is the significance level.
  #
  # Returns: A list containing the lower and upper intervals and the
  # significance level.
  Hlam <- get_Hlam(object$ipriorKernel, object$theta)
  if (is.null(Hlam.new)) Hlam.new <- Hlam
  se <- se_yhat(Hlam, Hlam.new, theta_to_psi(object$theta, object$ipriorKernel),
                is.nystrom(object))
  names(se) <- NULL
  lower <- y.hat + qnorm(alpha / 2) * se
  upper <- y.hat + qnorm(1 - alpha / 2) * se
  list(lower = lower, upper = upper, alpha = alpha)
}
