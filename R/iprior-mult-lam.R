###
### EM ALGORITHM
###

ipriorEM2 <- function(x, y, whichkernel=NULL, maxit=50000, delt=0.001, report.int=100, silent=F){
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
	lambda <- abs(rnorm(p))
	alpha <- rnorm(1)
	psi <- abs(rnorm(1))
	if(is.null(whichkernel)) whichkernel <- rep(F, p)
	
	### Define the kernel matrix
	H.mat <- NULL; H.matsq <- NULL; J.mat <- NULL
	for(j in 1:p){
		if(whichkernel[j])  H.mat[[j]] <- fn.H1(X[,j]) 
		else H.mat[[j]] <- fn.H2a(X[,j]) 
		H.matsq[[j]] <- H.mat[[j]] %*% H.mat[[j]]
	}
	J.mat <- function(k){
		tmp <- 0
		for(i in (1:p)[-k]){
			tmp <- tmp + lambda[i] * (H.mat[[i]] %*% H.mat[[k]] + H.mat[[k]] %*% H.mat[[i]])
		}
		tmp
	}
	
	H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	#Var.Y.inv <- chol2inv(chol(Var.Y))
	Var.Y.inv <- solve(Var.Y)
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
		for(j in 1:p){
			w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
			W.hat <- Var.Y.inv + tcrossprod(w.hat)
			T1 <- sum(H.matsq[[j]] * W.hat)
			T2 <- 2*crossprod(Y-alpha, crossprod(H.mat[[j]], w.hat)) - sum(J.mat(j) * W.hat)
			lambda[j] <- as.vector(T2/(2*T1))
			H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
			H.mat.lamsq <- H.mat.lam %*% H.mat.lam
			Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
			#Var.Y.inv <- chol2inv(chol(Var.Y))
			Var.Y.inv <- solve(Var.Y)	
		}	
	
		### Estimating psi using EM	
		w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		W.hat <- Var.Y.inv + tcrossprod(w.hat)
		T3 <- crossprod(Y-alpha) + sum(H.mat.lamsq * W.hat) - 2*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat))
		psi <- as.vector(N/T3)
	
		### New value of log-likelihood
		H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
		H.mat.lamsq <- H.mat.lam %*% H.mat.lam
		Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)	
		log.lik1 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
		check <- i %% report.int
		if(log.lik1 < log.lik0){
			cat("DECREASE iter", i, log.lik1, "\n")
		}
		else{
			if( !is.na(check) && check==0 && !silent ) cat("INCREASE iter", i, log.lik1, "\n") 
		} 
		
	}
	
	converged <- !(abs(log.lik0 - log.lik1) > delt)
	if(!silent && converged) cat("EM complete.\n", "\nNumber of iterations =", i, "\n")
	else if(!silent) cat("EM NOT CONVERGED!\n", "\nNumber of iterations =", i, "\n")
	if(!silent) cat("Log-likelihood = ", log.lik1, "\n")
	
	list(alpha=alpha, lambda=lambda, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, H.matsq=H.matsq, kernel=whichkernel, converged=converged, delt=delt)
}