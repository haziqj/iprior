###
### EM ALGORITHM
###

ipriorEM1 <- function(x, y, whichkernel=NULL, interactions=NULL, maxit=50000, delt=0.001, report.int=1000, silent=F){	
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
	lambda <- abs(rnorm(1, sd=0.01))
	alpha <- rnorm(1)
	psi <- abs(rnorm(1, sd=0.01))
	if(is.null(whichkernel)) whichkernel <- rep(F, p)
	
	### Define the kernel matrix
	H.mat <- NULL
	for(j in 1:p){
		if(whichkernel[j])  H.mat[[j]] <- fn.H1(X[,j]) 
		else H.mat[[j]] <- fn.H2a(X[,j]) 
	}
	## interactions
	if(!is.null(interactions)){
		Tmpo <- interactions[[1]]
		Tmpf <- interactions[[2]]
		no.int <- sum(Tmpo==2)
		p <- p + no.int
		for(j in (p-no.int+1):p){
			H.mat[[j]] <- H.mat[[ Tmpf[1, j-p+no.int] ]] * H.mat[[ Tmpf[2, j-p+no.int] ]]
		}
	}
	H.mat <- Reduce('+', mapply('*', H.mat, 1, SIMPLIFY=F))	
	H.matsq <- H.mat %*% H.mat

	Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
	#Var.Y.inv <- chol2inv(chol(Var.Y))
	Var.Y.inv <- solve(Var.Y)
	log.lik0 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
	
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
		w.hat <- psi * lambda * H.mat %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		Var.w.hat <- Var.Y.inv + w.hat %*% t(w.hat)
		T1 <- sum(H.matsq * Var.w.hat)
		T2 <- crossprod(Y-alpha, crossprod(H.mat, w.hat))
		lambda <- as.vector(T2/T1)
		
		### Estimating psi using EM	
		Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
		#Var.Y.inv <- chol2inv(chol(Var.Y))
		Var.Y.inv <- solve(Var.Y)
		w.hat <- psi * lambda * H.mat %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		Var.w.hat <- Var.Y.inv + w.hat %*% t(w.hat)
		T3 <- crossprod(Y-alpha) + lambda^2 * sum(H.matsq * Var.w.hat) - 2*lambda * crossprod(Y-alpha, crossprod(H.mat, w.hat))
		psi <- as.vector(N/T3)
	
		### New value of log-likelihood
		Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
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
	
	list(alpha=alpha, lambda=lambda, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, H.matsq=H.matsq, kernel=whichkernel, converged=converged, delt=delt)
}
