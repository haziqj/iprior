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

#' Obtain the log-likelihood and deviance of an I-prior model
#'
#' This function calculates the log-likelihood value or deviance (twice the
#' negative log-likelihood) for I-prior models. It works for both
#' \code{ipriorMod} and \code{ipriorKernel} class objects.
#'
#' For \code{ipriorKernel} objects, the log-likelihood or deviance is calculated
#' at the default parameter values: scale parameters and error precision are
#' equal to one, while hyperparameters of the kernels (e.g. Hurst index,
#' lengthscale, etc.) are the default values (see \link[=kernel]{here} for
#' details) or ones that has been specified. For \code{ipriorMod} objects, the
#' log-likelihood or deviance is calculated at the last obtained value from the
#' estimation method.
#'
#' For both types of objects, it is possible to supply parameter values at which
#' to calculate the log-likelihood/deviance. This makes estimating an I-prior
#' model more flexible, by first loading the variables into an
#' \code{ipriorKernel} object, and then using an optimiser such as
#' \code{\link[stats]{optim}}. Parameters have been transformed so that they can
#' be optimised unconstrained.
#'
#' @param object An object of class \code{ipriorMod} or \code{ipriorKernel}.
#' @param theta (Optional) Evaluates the log-likelihood at \code{theta}.
#' @param ... Not used.
#'
#' @seealso \link[=check_theta]{check_theta}.
#'
#' @export
logLik.ipriorMod <- function(object, theta = NULL, ...) {
  if (is.null(theta)) {
    res <- object$loglik
    return(res[length(res)])
  } else {
    logLik.ipriorKernel(object$ipriorKernel, theta)
  }
}

#' @rdname logLik.ipriorMod
#' @export
deviance.ipriorMod <- function(object, theta = NULL, ...) {
  -2 * logLik.ipriorMod(object, theta, ...)
}

#' @rdname logLik.ipriorMod
#' @export
logLik.ipriorKernel <- function(object, theta = NULL, ...) {
  if (is.null(theta)) theta <- object$thetal$theta
  if (is.nystrom(object)) {
    return(loglik_nystrom(theta, object))
  } else {
    return(loglik_iprior(theta, object))
  }
}

#' @rdname logLik.ipriorMod
#' @export
deviance.ipriorKernel <- function(object, theta = NULL, ...) {
  -2 * logLik.ipriorKernel(object, theta, ...)
}
