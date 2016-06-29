#' Function which returns the log-likelihood value of an \code{iprior} or
#' \code{ipriorKernel} object.
#'
#' \code{loglik} returns the log-likelihood value.
#'
#' \code{deviance} returns twice the negative log-likelihood value.
#'
#' Test
#'
#' @param x Object of class "iprior"
#' @param alpha (Optional) specify intercept value.
#' @param lambda (Optional) specify lambda value.
#' @param psi (Optional) specify psi value.
#' @return List of length 2.
#' @examples
#' mod.iprior <- iprior(stack.loss ~ ., data=stackloss)
#' loglik(mod.iprior)

#' @export
loglik <- function (object, ...) UseMethod("loglik")

#' @export
loglik.ipriorMod <- function (object, theta=NULL) {
	tmp <- with(object, {
		if (!is.null(theta)) {
			lambda <- theta[-length(theta)]
			psi <- theta[length(theta)]
		}
		ipriorEM(ipriorKernel, maxit=0, silent=T, lambda.init=lambda, psi.init=psi, clean=T)
	} )
	return(tmp$log.lik)
}

#' @export
loglik.ipriorKernel <- function (object, theta=NULL) {
	lambda <- theta[-length(theta)]
	psi <- theta[length(theta)]
	tmp <- ipriorEM(object, maxit=0, silent=T, lambda.init=lambda, psi.init=psi, clean=T)
	return(tmp$log.lik)
}

#' @export
deviance.ipriorMod <- function (object, theta=NULL) {
	return(-2*loglik(object, theta))
}

#' @export
deviance.ipriorKernel <- function (object, theta=NULL) {
	return(-2*loglik(object, theta))
}
