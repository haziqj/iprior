#' Function which returns the Variance of Y and log-likelihood value.
#'
#' \code{betas} returns the betas.
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
#' betas(mod.iprior)

betas <- function(object){
	lambda <- object$lambda
	w.hat <- object$w.hat
	p <- object$p
	x <- object$xval
	xnames <- colnames(x)

	# cts.vars <- which(!object$whichPearson)
	# ctg.vars <- which(object$whichPearson)
	# x.cts <- x[,cts.vars]
	# if(sum(object$whichPearson) == 1){
		# x.ctg <- x[,ctg.vars]
		# plotlvl <- levels(x.ctg)
		# grp <- as.numeric(x.ctg)
	# }
	
	Beta <- as.vector(lambda * (t(x) %*% w.hat))
	names(Beta) <- xnames
	Beta
}

