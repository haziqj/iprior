###
### EM ALGORITHM
###

## This is used mainly for parsimonious interactions

ipriorEM3 <- function(x, y, whichkernel=NULL, interactions=NULL, maxit=50000, delt=0.001, report.int=100, silent=F){
	### Library packages
	require(Matrix, quietly=T)			#to create diagonal matrices
	require(MASS, quietly=T)			#to sample from MVN dist.
	require(mvtnorm, quietly=T)
	require(numDeriv, quietly=T)

	X <- x
	Y <- y
	N <- length(Y)
	p <- ncol(X)
	x0 <- rep(1, N)
	lambda <- abs(rnorm(p, sd=0.01))
	alpha <- rnorm(1)
	psi <- abs(rnorm(1))
	if(is.null(whichkernel)) whichkernel <- rep(F, p)
	
	### Define the kernel matrix
	H.mat <- NULL#; H.matsq <- NULL
	for(j in 1:p){
		if(whichkernel[j])  H.mat[[j]] <- fn.H1(X[,j]) 
		else H.mat[[j]] <- fn.H2a(X[,j]) 
	}
	
	## interactions (parsimonious way)
	## this should be a list of 2: vector or "orders" and matrix of "factors"
	if(!is.null(interactions)){
		Tmpo <- interactions[[1]]
		Tmpf <- interactions[[2]]
		no.int <- sum(Tmpo==2)
		p1 <- p + no.int; lambda.int <- lambda
		for(j in (p1-no.int+1):p1){
			H.mat[[j]] <- H.mat[[ Tmpf[1, j-p1+no.int] ]] * H.mat[[ Tmpf[2, j-p1+no.int] ]]
			lambda.int <- c(lambda.int, lambda[Tmpf[1, j-p1+no.int]]*lambda[Tmpf[2, j-p1+no.int]])
		}
	}
	else lambda.int <- lambda
	
	H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda.int, SIMPLIFY=F))
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	#Var.Y.inv <- chol2inv(chol(Var.Y))
	Var.Y.inv <- solve(Var.Y)
	log.lik0 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
	
	## This is the E-step for lambda which needs to be maximised
	Q.lam <- function(lam){
		lambda.int <- lam		
		if(!is.null(interactions)){
			for(j in (p1-no.int+1):p1) lambda.int <- c(lambda.int, lam[Tmpf[1, j-p1+no.int]]*lam[Tmpf[2, j-p1+no.int]])
		} 
		H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda.int, SIMPLIFY=F))
		H.mat.lamsq <- H.mat.lam %*% H.mat.lam
		Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
		res <- psi/2 * crossprod(Y - alpha) - psi * crossprod(Y - alpha, crossprod(H.mat.lam, w.hat)) + 1/2 * sum(Var.Y * W.hat)
		res
	}
	
	## initialise
	if(!silent) cat("START iter", 0, log.lik0, "\t")
	log.lik1 <- log.lik0 + 2*delt
	i <- 0
	pb <- txtProgressBar(min=0, max=report.int*10, style=1, char=".") #progress bar

	while((i != maxit) && (abs(log.lik0 - log.lik1) > delt)){
	
		i <- i + 1
		log.lik0 <- log.lik1
		
		### Estimating alpha
		tmp.alpha <- crossprod(x0, Var.Y.inv)
		alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))
		
		### Estimating lambda using EM
		w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		W.hat <- Var.Y.inv + tcrossprod(w.hat)
		lambda <- optim(lambda, fn=Q.lam, method="Nelder-Mead")$par
		lambda.int <- lambda
		if(!is.null(interactions)){
			for(j in (p1-no.int+1):p1) lambda.int <- c(lambda.int, lambda[Tmpf[1, j-p1+no.int]]*lambda[Tmpf[2, j-p1+no.int]])
		} 
		H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda.int, SIMPLIFY=F))		
		H.mat.lamsq <- H.mat.lam %*% H.mat.lam
		Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)	
		Var.Y.inv <- solve(Var.Y)	
		
		### Estimating psi using EM	
		w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		W.hat <- Var.Y.inv + tcrossprod(w.hat)
		T3 <- crossprod(Y-alpha) + sum(H.mat.lamsq * W.hat) - 2*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat))
		psi <- as.vector(N/T3)
	
		### New value of log-likelihood
		lambda.int <- lambda
		if(!is.null(interactions)){
			for(j in (p1-no.int+1):p1) lambda.int <- c(lambda.int, lambda[Tmpf[1, j-p1+no.int]]*lambda[Tmpf[2, j-p1+no.int]])
		} 
		H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda.int, SIMPLIFY=F))	
		H.mat.lamsq <- H.mat.lam %*% H.mat.lam
		Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)	
		log.lik1 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
		
		### Report
		check <- i %% report.int
		if(log.lik1 < log.lik0){
			cat("\nDECREASE iter", i, log.lik1, "\t")
		}
		else{
			if( !is.na(check) && check==0 && !silent ) cat("\nINCREASE iter", i, log.lik1, "\t") 
		} 
		setTxtProgressBar(pb, i)
		if(i %% report.int*10 == 0) pb <- txtProgressBar(min=i, max=report.int*10+i, style=1, char=".") 
			#reset progress bar
	}
	
	close(pb)
	converged <- !(abs(log.lik0 - log.lik1) > delt)
	if(!silent && converged) cat("EM complete.\n", "\nNumber of iterations =", i, "\n")
	else if(!silent) cat("EM NOT CONVERGED!\n", "\nNumber of iterations =", i, "\n")
	if(!silent) cat("Log-likelihood = ", log.lik1, "\n")
	
	list(alpha=alpha, lambda=lambda, lambda.int=lambda.int, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, kernel=whichkernel, converged=converged, delt=delt)
}