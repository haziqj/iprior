#' Function which returns the Variance of Y and log-likelihood value.
#'
#' \code{ipriorloglik} returns the the variance of Y and the log-likelihood value.
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
#' ipriorloglik(mod.iprior)

ipriorloglik <- function(x, alpha=NULL, lambda=NULL, psi=NULL){
	require(mvtnorm, quietly=T)
	if(class(x) != "iprior") stop("Input iprior class models only.", call.=F)
	N <- length(x$yval)
	if(is.null(alpha)) alpha <- x$alpha
	if(is.null(lambda)) lambda <- x$lambda
	else if(length(lambda) != length(x$lambda)) stop(paste("Incorrect dimension of lambda values. vector of length", length(x$lambda), "required."), call.=F)
	if(is.null(psi)) psi <- x$psi
	
	H.mat.lam <- Reduce('+', mapply('*', x$H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	log.lik <- dmvnorm(x$yval, mean=rep(alpha,N), sigma=Var.Y, log=T)
	list(vary=Var.Y, loglik=log.lik)
}

ipriorloglik2 <- function(x, alpha=NULL, lambda=NULL, psi=NULL){
	#require(mvtnorm, quietly=T)
	if(class(x) != "iprior") stop("Input iprior class models only.", call.=F)
	N <- length(x$yval)
	if(is.null(alpha)) alpha <- x$alpha
	if(is.null(lambda)) lambda <- x$lambda
	else if(length(lambda) != length(x$lambda)) stop(paste("Incorrect dimension of lambda values. vector of length", length(x$lambda), "required."), call.=F)
	if(is.null(psi)) psi <- x$psi
	
	H.mat.lam <- Reduce('+', mapply('*', x$H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	tmp <- chol(Var.Y)
	log.lik <- -N/2*log(2*pi) - sum(log(diag((tmp2 + t(tmp2)/2 )))) - 1/2*as.vector(crossprod(x$yval-alpha, solve(Var.Y, x$yval-alpha)))
	log.lik <- -N/2*log(2*pi) - 1/2*log(det(Var.Y)) - 1/2*t(x$yval-alpha) %*% solve(Var.Y) %*% (x$yval-alpha)
	list(vary=Var.Y, loglik=log.lik)
}

ipriorloglik3 <- function(x, alpha=NULL, lambda=NULL, psi=NULL){
	require(mvnfast, quietly=T)
	if(class(x) != "iprior") stop("Input iprior class models only.", call.=F)
	N <- length(x$yval)
	if(is.null(alpha)) alpha <- x$alpha
	if(is.null(lambda)) lambda <- x$lambda
	else if(length(lambda) != length(x$lambda)) stop(paste("Incorrect dimension of lambda values. vector of length", length(x$lambda), "required."), call.=F)
	if(is.null(psi)) psi <- x$psi
	
	H.mat.lam <- Reduce('+', mapply('*', x$H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	tmp <- chol(Var.Y)
	log.lik <- dmvn(x$yval, mu=rep(alpha, N), sigma=Var.Y, log=T)
	list(vary=Var.Y, loglik=log.lik)
}