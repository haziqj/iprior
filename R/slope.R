#' Function which returns the Variance of Y and log-likelihood value.
#'
#' \code{slope} returns the slopes.
#'
#' Details.
#'
#' @param ipriorMod Object of class "ipriorMod"
#'
#' @return What.
#'
#' @examples
#' mod.iprior <- iprior(stack.loss ~ Air.Flow, data=stackloss)
#' slope(mod.iprior)
#'
#' @export
slope <- function(ipriorMod){
	if (!is.ipriorMod(ipriorMod)) stop("Input iprior class models only.", call.=F)

	y <- ipriorMod$fitted.values
	x <- ipriorMod$ipriorKernel$x
	whichPearson <- ipriorMod$ipriorKernel$whichPearson
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

