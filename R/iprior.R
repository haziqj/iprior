

## classic Fisher cats data from package MASS
data(cats, package="MASS")

## Canonical kernel function
fn.H2 <- function(x, y=NA){ #takes in vector of covariates
	if(is.na(sum(y))) y <- x
	z <- c(x,y)
	tmp <- tcrossprod(y, x)
	tmp
}

## (centred) Canonical kernel function
fn.H2a <- function(x, y=NA){ #takes in vector of covariates
	if(is.na(sum(y))) y <- x
	z <- c(x,y)
	xbar <- mean(z)
	tmp <- tcrossprod(y-xbar, x-xbar)
	tmp
}

###
### EM ALGORITHM
###

ipriorEM <- function(x, y, maxit=50000, delt=0.00001, report.int=1000, silent=F){
	### Library packages
	require(Matrix)			#to create diagonal matrices
	require(MASS)			#to sample from MVN dist.
	require(mvtnorm)
	require(numDeriv)
	
	X <- x
	Y <- y
	N <- length(Y)
	p <- ncol(X)
	x0 <- rep(1, N)
	lambda <- abs(rnorm(1))
	alpha <- rnorm(1)
	psi <- abs(rnorm(1))
	
	### Define the kernel matrix
	H.mat <- 0
	for(j in 1:p){
		H.mat <- H.mat + fn.H2a(X[,j]) 
	}
	H.matsq <- H.mat %*% H.mat

	Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
	Var.Y.inv <- chol2inv(chol(Var.Y))
	#Var.Y.inv <- solve(Var.Y)
	log.lik0 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
	if(!silent) cat("START iter", 0, log.lik0, "\n")
	log.lik1 <- log.lik0 + 2*delt
	i <- 0
	
	while((i != maxit) && (abs(log.lik0 - log.lik1) > delt)){
	
	i <- i + 1
	log.lik0 <- log.lik1
	
	### Estimating alpha
	tmp.alpha <- crossprod(x0, Var.Y.inv)
	alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))
	
	### Estimating lambda using EM
	w.hat <- psi * lambda * H.mat %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
	Var.w.hat <- Var.Y.inv + w.hat %*% t(w.hat)
	T1 <- sum(H.matsq * Var.w.hat)
	T2 <- crossprod(Y-alpha, crossprod(H.mat, w.hat))
	lambda <- as.vector(T2/T1)
	
	### Estimating psi using EM	
	Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
	Var.Y.inv <- chol2inv(chol(Var.Y))
	w.hat <- psi * lambda * H.mat %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
	Var.w.hat <- Var.Y.inv + w.hat %*% t(w.hat)
	# T4 <- tr(Var.Y.inv) + crossprod(w.hat)
	# psi <- as.vector(T4/N)
	T3 <- crossprod(Y-alpha) + lambda^2 * sum(H.matsq * Var.w.hat) - 2*lambda * crossprod(Y-alpha, crossprod(H.mat, w.hat))
	psi <- as.vector(N/T3)

	### New value of log-likelihood
	Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
	log.lik1 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
	check <- i %% report.int
	if( !is.na(check) && check==0 && !silent ){
		if(log.lik1 >= log.lik0) cat("INCREASE iter", i, log.lik1, "\n") 
		else cat("DECREASE iter", i, log.lik1, "\n")
	}
	
	}

	if(!silent) cat("EM complete.\n", "\nNumber of iterations =", i, "\n")
	if(!silent) cat("Log-likelihood = ", log.lik1, "\n")
	
	list(alpha=alpha, lambda=lambda, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, H.matsq=H.matsq)
}

###
### The generic function is a standard R function with a special body, usually containing 
### only a call to UseMethod
###

iprior <- function(x, y, ...) UseMethod("iprior")

## The default method
iprior.default <- function(x, y, maxit=50000, delt=0.00001, report.int=1000, silent=F, ...){
	x <- as.matrix(x)
	y <- as.numeric(y)
	n <- length(y)
	
	est <- ipriorEM(x, y, maxit=maxit, delt=delt, report.int=report.int, silent=silent)
	param <- c(est$alpha, est$lambda, est$psi)
	names(param) <- c("alpha", "lambda", "psi")
	
	##calculate fitted values
	Var.Y <- est$lambda * est$lambda * est$psi * est$H.matsq + (1/est$psi) * diag(n)
	Var.Y.inv <- solve(Var.Y)
	w.hat <- (est$psi * est$lambda) * est$H.mat %*% (Var.Y.inv %*% matrix(y - est$alpha, nc=1))
	Var.w.hat <- Var.Y.inv + tcrossprod(w.hat)
	Y.hat <- est$alpha + est$lambda * as.vector(crossprod(est$H.mat, w.hat))
	
	est$call <- match.call()
	est$coefficients <- param
	est$fitted.values <- Y.hat
	est$w.hat <- w.hat
	est$yval <- y
	est$xval <- x
	
	class(est) <- "iprior"
	est
}

print.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	cat("\nParameter estimates:\n")
	print(x$coefficients)
	cat("\n")
}

## The summary screen
summary.iprior <- function(object, ...){
	## Calculate standard errors through Inverse (observed) Fisher information matrix
	ipriorloglik <- function(theta){
		alpha <- theta[1]; lambda <- theta[-c(1,length(theta))]; psi <- theta[length(theta)]
		n <- length(fitted(object))
		y <- object$yval
		Var.Y <- lambda * lambda * psi * object$H.matsq + (1/psi) * diag(n)
		loglik <- dmvnorm(y-alpha, rep(0,n), Var.Y, log=T)
		loglik
	}
	FisherInformation <- -hessian(ipriorloglik, coef(object))
	InverseFisher <- solve(FisherInformation)
	se <- sqrt(diag(InverseFisher))
	
	## Z values to compare against (standard) Normal distribution
	zval <- coef(object) / se
	
	## Create table
	tab <- cbind(	Estimate=coef(object),
					S.E.=se,
					z=zval,
					"P[|Z>z|]"=2*pnorm(-abs(zval)) )
					
	res <- list(call=object$call, coefficients=tab)
	
	class(res) <- "summary.iprior"
	res
}

print.summary.iprior <- function(x, ...){
	cat("Call:\n")
	print(x$call)
	cat("\n")
	printCoefmat(x$coefficients, P.value=T, has.Pvalue=T)
}

## Formulas
iprior.formula <- function(formula, data=list(), ...){
	mf <- model.frame(formula=formula, data=data)
	attr(attr(mf, "terms"), "intercept") <- 0
	x <- model.matrix(attr(mf, "terms"), data=mf)
	y <- model.response(mf)
	
	est <- iprior.default(x, y, ...)
	est$call <- match.call()
	est$formula <- formula
	names(est$fitted.values) <- row.names(mf)
	est
}

## Prediction
predict.iprior <- function(object, newdata=NULL, ...){
	if(is.null(newdata)) ystar <- fitted(object)
	else{
		if(!is.null(object$formula)){
			## model has been fitted using formula interface
			xstar <- model.matrix(object$formula, newdata)
			xcolnames <- colnames(xstar); xrownames <- rownames(xstar)
			wheres.int <- (colnames(xstar) == "(Intercept)")
			xstar <- xstar[,!wheres.int]
			x <- as.matrix(xstar)
			#colnames(xstar) <- xnames[!wheres.int]
		}
		else{
			xstar <- as.matrix(newdata)
		}
		
		## Define new kernel matrix
		Y <- object$y; X <- object$x; p <- ncol(X)
		H.mat <- 0
		for(j in 1:p){
			H.mat <- H.mat + fn.H2a(X[,j], xstar) 
		}
		ystar <- as.vector(object$alpha + object$lambda * (H.mat %*% object$w.hat))
		names(ystar) <- xrownames
	}
	ystar
}

###
### Comparison of I-prior and OLS
### 
# mod.iprior <- iprior(Hwt~Bwt*Sex, data=cats)
# mod.lm <- lm(Hwt~Bwt*Sex, data=cats)

# ### Compare value of sigma
# sigma.iprior <-  1/sqrt(mod.iprior$psi)
# sigma.lm <- summary(mod.lm)$sigma
# sigma.iprior; sigma.lm

# ### classical simple regression model
# fit.iprior <- fitted(mod.iprior)
# fit.lm <- fitted(mod.lm)

# ### comparing fitted values
# plot(fit.lm, fit.iprior, type="n", xlab="Classical regression model estimates", ylab="I-prior estimates", main="Comparison between I-prior and classical regression predicted values")
# text(fit.lm, fit.iprior, pch=as.character(1:length(fit.lm)), col=1:length(fit.lm), cex=0.7)
# abline(a=0, b=1)

# ### getting the betas back
# beta.iprior <- as.vector(lambda*crossprod(X[,1:p], w.hat))
# beta.iprior
# beta.ols[-1]

## second test
# beta.true <- matrix(c(rep(0,10), rep(1,40)), nc=1)
# n <- 200
# p <- length(beta.true)
# X <- matrix(rnorm(n*p), nr=n)
# Y <- X %*% beta.true + rnorm(n, mean=0, sd=2); Y <- as.vector(Y)
# mod.iprior <- iprior(Y~X)
# mod.lm <- lm(Y~1+X)
