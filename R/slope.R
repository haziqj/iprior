#' Function which returns the Variance of Y and log-likelihood value.
#'
#' \code{slope} returns the slopes.
#'
#' Test
#'
#' @param x Object of class "iprior"
#' @param alpha (Optional) specify intercept value.
#' @param lambda (Optional) specify lambda value.
#' @param psi (Optional) specify psi value.
#' @return List of length 2.
#' @examples
#' mod.iprior <- iprior(stack.loss ~ Air.Flow, data=stackloss)
#' slope(mod.iprior)

slope <- function(object){
	if (!is(object, "iprior")) stop("Input iprior class models only.", call.=F)

	y <- object$fitted.values
	x <- object$ipriorKernel$x
	whichPearson <- object$ipriorKernel$whichPearson
	dat <- cbind(y, as.data.frame(x))
	if (any(whichPearson)) {
		grp <- dat[[which(whichPearson)+1]]
		dat <- split(dat, grp)
		res <- sapply(dat, function(x) coef(lm(x[,1]~x[,2]))[2] )
		names(res) <- levels(grp)
	} else {
		res <- coef(lm(dat[,1]~dat[,2]))[2]
		names(res) <- "slope"
	}
	res
}

