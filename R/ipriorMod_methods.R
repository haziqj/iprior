################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2017  Haziq Jamil
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

#' Print and summary method for I-prior models
#'
#' @param object,x An \code{ipriorMod} object.
#' @param digits Number of decimal places for the printed coefficients.
#' @param ... Not used.
#'
#' @name summary.ipriorMod
#' @export
NULL

#' @rdname summary.ipriorMod
#' @export
print.ipriorMod <- function(x, digits = 5, ...) {
  loglik.max <- x$loglik[length(x$loglik)]
  cat("Log-likelihood value:", loglik.max, "\n")
  cat("\n")
  if (x$ipriorKernel$thetal$n.theta > 0)
    print(round(coef(x), digits))
  else
    cat("No hyperparameters estimated.")
}

#' @rdname summary.ipriorMod
#' @export
summary.ipriorMod <- function(object, ...) {
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
    xs <- paste0(object$ipriorKernel$xname[ind], collapse = ", ")
    x.kern[i] <- paste0(unique.kernels[i], " (", xs, ")\n")
  }

  if (object$method == "mixed") {
    maxit <- object$control$maxit + object$control$em.maxit
    niter <- object$niter + object$control$em.maxit
  } else {
    maxit <- object$control$maxit
    niter <- object$niter
  }

  res <- list(resid.summ = resid.summ, tab = tab, loglik = logLik(object),
              error = object$train.error, call = object$call, x.kern = x.kern,
              est.method = object$est.method, est.conv = object$est.conv,
              niter = niter, maxit = maxit, time = object$time)
  class(res) <- "ipriorMod_summary"
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

#' @export
print.ipriorMod_summary <- function(x, ...) {
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

if (getRversion() < "3.3.0") {
  sigma <- function(object, ...) UseMethod("sigma")
}

#' Obtain the standard deviation of the residuals 'sigma'
#'
#' Extract the standard deviation of the residuals. For I-prior models, this is
#' \code{sigma = 1 / sqrt(psi)}.
#'
#' @param object An object of class \code{ipriorMod}.
#' @param ... Not used.
#'
#' @rawNamespace if (getRversion() >= "3.3.0") importFrom(stats,sigma)
#' @rawNamespace if (getRversion() < "3.3.0") export(sigma)
#' @name sigma
#' @export
sigma.ipriorMod <- function(object, ...) {
  psi <- theta_to_psi(object$theta, object$ipriorKernel)
  res <- 1 / sqrt(psi)
  names(res) <- "sigma"
  res
}

#' Update an I-prior model
#'
#' @param object An \code{ipriorMod} object.
#' @param method An optional method. See \link[=iprior]{here} for details.
#' @param control An optional list of controls for the estimation procedure. See
#'   \link[=iprior]{here} for details.
#' @param iter.update The number of additional iterations to update the I-prior
#'   model.
#' @param ... Not used.
#'
#' @export
update.ipriorMod <- function(object, method = NULL, control = list(),
                             iter.update = 100, ...) {
  res <- iprior.ipriorMod(object, method, control, iter.update, ...)
  assign(deparse(substitute(object)), res, envir = parent.frame())
}

#' Extract the kernel matrix from I-prior models
#'
#' @param object An \code{ipriorMod} or \code{ipriorKernel2} object.
#' @param theta (Optional) Value of hyperparameters to evaluate the kernel
#'   matrix.
#' @param xstar (Optional) If not supplied, then a square, symmetric kernel
#'   matrix is returned using the data as input points. Otherwise, the kernel
#'   matrix is evaluated with respect to this set of data as well. It must be a
#'   list of vectors/matrices with similar dimensions to the original data.
#'
#' @return A kernel matrix.
#'
#' @export
get_kern_matrix <- function(object, theta = NULL, xstar = list(NULL)) {
  if (is.ipriorMod(object)) {
    # estl <- object$ipriorKernel$estl
    # til.cond <- (
    #   !isTRUE(estl$est.hurst) & !isTRUE(estl$est.lengt) & !isTRUE(estl$est.offs)
    # )
    res <- get_Hlam(object$ipriorKernel, object$theta, FALSE)
    return(res)
  } else if (is.ipriorKernel2(object)) {
    # estl <- object$estl
    # til.cond <- (
    #   !isTRUE(estl$est.hurst) & !isTRUE(estl$est.lengt) & !isTRUE(estl$est.offs)
    # )
    res <- get_Hlam(object, object$theta, FALSE)
    return(res)
  }
}