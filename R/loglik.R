################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2016  Haziq Jamil
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

#' Obtain the log-likelihood and deviance of an \code{ipriorMod} or
#' \code{ipriorKernel} object
#'
#' This function calculates the log-likelihood value or deviance (twice the
#' negative log-likelihood) for I-prior models. It works for both
#' \code{ipriorMod} and \code{ipriorKernel} class objects.
#'
#' For \code{ipriorKernel} objects, the log-likelihood or deviance is calculated
#' at \emph{random parameter values} by default. This is because it uses the
#' (random) initial values for the parameters from the EM algorithm. For
#' \code{ipriorMod} objects, the default value is the last obtained set of
#' parameter values of the EM algorithm. If the model has converged, then this
#' should be the maximum likelihood (ML) value.
#'
#' For both types of object, it is possible to supply a parameter value at which
#' to calculate the log-likelihood/deviance. This makes estimating an I-prior
#' model more flexible, by first loading the variables into an
#' \code{ipriorKernel} object, and then using an optimiser such as
#' \code{\link[stats]{optim}}. Any other optimiser that can implement bounded
#' constraints can be used. The only constraint requirement is on the parameter
#' \code{psi} which has to be greater than zero (although setting the constraint
#' to a small number such as \code{1e-9} will avoid an ill conditioned
#' variance).
#'
#' Direct maximisation of the log-likelihood of I-prior models can be unreliable
#' some times, but other times it could also lead to a better set of parameters
#' if the EM does not converge fully. Perhaps in the future, this package can
#' make use of both methods combined as standard.
#'
#' As a side note, the ML estimate for the intercept of an I-prior model is the
#' mean of the responses \code{y}, and therefore not included as a parameter of
#' the log-likelihood/deviance function.
#'
#' @param object An object of class \code{ipriorMod} or \code{ipriorKernel}.
#' @param theta (optional) Evaluates the log-likelihood at \code{theta} which is
#'   of the form \code{theta = c(lambda, psi)}.
#' @param ... No further arguments required, so this is not used for
#'   \code{ipriorMod} or \code{ipriorKernel} class objects.
#'
#' @examples
#' mod.iprior <- iprior(stack.loss ~ ., data=stackloss)
#' logLik(mod.iprior)
#'
#' # Using optim to find ML estimates
#' optim(abs(rnorm(4)), deviance, object = mod.iprior,
#'       method="L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 1e-9))
#' coef(mod.iprior)
#' deviance(mod.iprior)
#'
#' @name logLik
#' @rdname logLik
#' @export
logLik.ipriorMod <- function(object, theta = NULL, ...) {
	tmp <- with(object, {
		if (!is.null(theta)) {
			lambda <- theta[-length(theta)]
			psi <- theta[length(theta)]
		}
		ipriorEM(ipriorKernel, maxit = 0, silent = TRUE, lambda.init = lambda,
		         psi.init = psi, clean = TRUE)
	} )
	return(tmp$log.lik)
}

#' @rdname logLik
#' @export
logLik.ipriorKernel <- function(object, theta = NULL, ...) {
	lambda <- theta[-length(theta)]
	psi <- theta[length(theta)]
	tmp <- ipriorEM(object, maxit = 0, silent = TRUE, lambda.init = lambda,
	                psi.init = psi, clean = TRUE)
	return(tmp$log.lik)
}

#' @name deviance
#' @rdname logLik
#' @export
deviance.ipriorMod <- function(object, theta = NULL, ...) {
	return(-2 * logLik(object, theta))
}

#' @rdname logLik
#' @export
deviance.ipriorKernel <- function(object, theta = NULL, ...) {
	return(-2 * logLik(object, theta))
}
