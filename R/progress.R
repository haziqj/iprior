###
### Function to check progress of EM
###

#' @export
progress <- function(x, interval=c("auto", "all", "input any number")){
	if(class(x) != "iprior") stop("Input iprior class models only.", call.=F)
	if(!x$converged) warning("The EM has not yet converged.", call.=F)
	rn <- rownames(x$res.loglik)

	### Parameters
	res <- x$res.loglik
	cn <- c("Log-likelihood", "Pred.log-lik.", "Delta(i,i-1)", colnames(x$res.param)[-1])
	res <- cbind(res, x$res.param[,-1])
	res <- as.data.frame(res, row.names=rn)
	colnames(res) <- cn

	### Trim the table
	no.iter <- x$no.iter
	if(!is.numeric(interval)) interval <- match.arg(interval)
	if(interval == "auto") interval <- max(no.iter %/% 8, 1)
	if(interval == "all") interval <- 1
	if(interval == "input any number") stop("No, what I meant was that interval should be numeric! Try interval=10.", call.=F)
	trim <- no.iter %/% interval
	first <- res[1,]; last <- res[no.iter+1,]
	res <- res[seq(from=interval+1, to=(trim*interval)+1, by=interval),]
	if(no.iter %% interval != 0) res <- rbind(first, res, last)
	else res <- rbind(first, res)
	res
}
