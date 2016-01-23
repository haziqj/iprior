###
### EM ALGORITHM
###

ipriorEM2 <- function(x, y, whichkernel=NULL, interactions=NULL, maxit=50000, stop.crit=0.001, report.int=1000, silent=F){
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
	lambda <- abs(rnorm(p, sd=0.1))
	alpha <- rnorm(1)
	if(is.null(whichkernel)) whichkernel <- rep(F, p)
	
	### Define the kernel matrix
	H.mat <- NULL; H.matsq <- NULL
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
	
	## interactions (not parsimonious)
	## this should be a list of 2: vector or "orders" and matrix of "factors"
	if(!is.null(interactions)){
		Tmpo <- interactions[[1]]
		Tmpf <- interactions[[2]]
		no.int <- sum(Tmpo==2)
		p <- p + no.int; lambda <- abs(rnorm(p))
		for(j in (p-no.int+1):p){
			H.mat[[j]] <- H.mat[[ Tmpf[1, j-p+no.int] ]] * H.mat[[ Tmpf[2, j-p+no.int] ]]
			H.matsq[[j]] <- H.mat[[j]] %*% H.mat[[j]]
		}
	}
	
	H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	#Var.Y.inv <- chol2inv(chol(Var.Y))
	Var.Y.inv <- solve(Var.Y)
	log.lik0 <- dmvnorm(Y-alpha, rep(0,N), Var.Y, log=T)
	
	## initialise
	if(!silent) cat("START iter", 0, log.lik0, "\t")
	log.lik1 <- log.lik0 + 2*stop.crit
	i <- 0
	if(!silent) pb <- txtProgressBar(min=0, max=report.int*10, style=1, char=".") #progress bar
	
	while((i != maxit) && (abs(log.lik0 - log.lik1) > stop.crit)){
	
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
		Var.Y.inv <- solve(Var.Y)
		log.lik1 <- dmvnorm(Y, mean=rep(alpha,N), sigma=Var.Y, log=T)
		
		### Report
		check <- i %% report.int
		if(log.lik1 < log.lik0 && !silent){
			cat("\nDECREASE iter", i, log.lik1, "\t")
		}
		else{
			if(!is.na(check) && check==0 && !silent) cat("\nINCREASE iter", i, log.lik1, "\t") 
		} 
		if(!silent) setTxtProgressBar(pb, i)
		if(i %% report.int*10 == 0 && !silent) pb <- txtProgressBar(min=i, max=report.int*10+i, style=1, char=".") 
			#reset progress bar
	}
	
	if(!silent) close(pb)
	converged <- !(abs(log.lik0 - log.lik1) > stop.crit)
	if(!silent && converged) cat("EM complete.\n", "\nNumber of iterations =", i, "\n")
	else if(!silent) cat("EM NOT CONVERGED!\n", "\nNumber of iterations =", i, "\n")
	if(!silent) cat("Log-likelihood = ", log.lik1, "\n")
	
	list(alpha=alpha, lambda=lambda, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, H.matsq=H.matsq, kernel=whichkernel, converged=converged, stop.crit=stop.crit)
}