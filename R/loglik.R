#' Function which returns the log-likelihood value of an \code{ipriorMod} or
#' \code{ipriorKernel} object.
#'
#' Description.
#'
#' Details.
#'
#' @param object an object of class \code{ipriorMod} or \code{ipriorKernel}.
#' @param theta (optional) evaluates the log-likelihood at \code{theta} which is
#'   of the form \code{theta = c(lambda, psi)}.
#' @param ... This is not used here.
#'
#' @examples
#' mod.iprior <- iprior(stack.loss ~ ., data=stackloss)
#' logLik(mod.iprior)
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

#' @export
deviance.ipriorMod <- function(object, theta = NULL, ...) {
	return(-2 * logLik(object, theta))
}

#' @export
deviance.ipriorKernel <- function(object, theta = NULL, ...) {
	return(-2 * logLik(object, theta))
}
