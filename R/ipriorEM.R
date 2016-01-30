###
### The iprior EM
###

ipriorEM <- function(x, y, whichkernel, interactions, one.lam, parsm, kernel, maxit, stop.crit, report.int, silent, alpha.init, lambda.init, psi.init, invmethod){
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
	if(report.int == 0)	report.int <- maxit

	### Check for interactions and determine q = length of lambda
	if(!is.null(interactions)){
		Tmpo <- interactions[[1]]
		Tmpf <- interactions[[2]]
		no.int <- sum(Tmpo==2)
		if(parsm) q <- p		#parsimonious interactions
		else q <- p + no.int	#non-parsimonious interactions
	}
	else{ #no interactions
		no.int <- 0
		q <- p
	} 			
	if(one.lam){			#one lambda
		q <- 1
		if(parsm) message("Parsimonious interactions ignored because one.lam=T", call.=F)
	} 		

	### Initialise parameters
	alpha <- alpha.init
	if(is.null(lambda.init)) lambda <- abs(rnorm(q, sd=0.1))
	else{
		if(length(lambda.init) != q) stop(paste("Incorrect dimension of lambda initial values. vector of length", q, "required."), call.=F)
		else lambda <- lambda.init
	}
	psi <- psi.init
	
	### Define the kernel matrix
	H.mat <- NULL; H.matsq <- NULL
	for(j in 1:p){
		if(whichkernel[j]) H.mat[[j]] <- fn.H1(X[,j])		#Pearson
		else{
			if(kernel=="FBM") H.mat[[j]] <- fn.H3(X[,j])	#FBM
			else H.mat[[j]] <- fn.H2a(X[,j])				#Canonical
		}
		if(q > 1)  H.matsq[[j]] <- H.mat[[j]] %*% H.mat[[j]]
	}
	if(!is.null(interactions)){
		for(j in 1:no.int){
			H.mat[[p+j]] <- H.mat[[ Tmpf[1,j] ]] * H.mat[[ Tmpf[2,j] ]]
			H.matsq[[p+j]] <- H.mat[[p+j]] %*% H.mat[[p+j]]
		}
	}
	if(q == 1){	#just to save some time, if q=1 then no need to loop
		H.matsq <- Reduce('+', mapply('*', H.mat, 1, SIMPLIFY=F))
		H.matsq <- list(H.matsq %*% H.matsq)
	}
	##these are the two-way terms in the H^2 matrix
	w <- 1:(p+no.int)
	ind1 <- rep(w, times=(length(w)-1):0)
	ind2 <- unlist(lapply(2:length(w), function(x) c(NA,w)[-(0:x)]))
	H.mat2 <- NULL
	for(j in 1:length(ind1)){
		H.mat2[[j]] <- H.mat[[ ind1[j] ]] %*% H.mat[[ ind2[j] ]] + H.mat[[ ind2[j] ]] %*% H.mat[[ ind1[j] ]]
	}
	##this function is required to calculate the sufficient statistics
	J.mat <- function(k){	#always for non-parsimonious method q = p + no.int
		ind <- which(ind1==k | ind2==k)
		tmp <- Reduce('+', mapply('*', H.mat2[ind], lambda[-k], SIMPLIFY=F))
		tmp
	}
	
	### Initialise the EM
	if(!is.null(interactions) && parsm){
		for(j in 1:no.int) lambda <- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]])
		### This is the E-step for lambda which needs to be maximised in the parsimonious method
		Q.lam <- function(lambda){
			for(j in 1:no.int) lambda <- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]])
			H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
			H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
			Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
			res <- psi/2*crossprod(Y-alpha) - psi*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat)) + 1/2*sum(Var.Y*W.hat)
			res
		}
	}
	H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
	lambda <- lambda[1:q]
	Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
	ifelse(invmethod == "eigen", Var.Y.inv <- solve(Var.Y), Var.Y.inv <- chol2inv(chol(Var.Y)) )
	log.lik0 <- dmvnorm(Y, rep(alpha,N), Var.Y, log=T)
	
	if(!silent && report.int != maxit) cat("Iteration 0:\t Log-likelihood =", round(log.lik0, 4), " ")
	if(!silent) pb <- txtProgressBar(min=0, max=report.int*10, style=1, char=".") #progress bar
	log.lik1 <- log.lik0 + 2*stop.crit
	i <- 0
	
	while((i != maxit) && (abs(log.lik0 - log.lik1) > stop.crit)){
	
		i <- i + 1
		log.lik0 <- log.lik1
				
		### Estimating alpha
		tmp.alpha <- crossprod(x0, Var.Y.inv)
		alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))
		
		### Estimating lambda using EM
		if(!is.null(interactions) && parsm){ #parsimonious method
			w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
			W.hat <- Var.Y.inv + tcrossprod(w.hat)
			lambda <- optim(lambda, fn=Q.lam, method="Nelder-Mead")$par
			for(j in 1:no.int) lambda <- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]])
			H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))		
			H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
			lambda <- lambda[1:q]
			Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)	
			ifelse(invmethod == "eigen", Var.Y.inv <- solve(Var.Y), Var.Y.inv <- chol2inv(chol(Var.Y)) )
		}
		else{ #non-parsimonious method
			for(j in 1:q){
				w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
				W.hat <- Var.Y.inv + tcrossprod(w.hat)
				T1 <- sum(H.matsq[[j]] * W.hat)
				T2 <- 2*crossprod(Y-alpha, crossprod(H.mat[[j]], w.hat)) - sum(J.mat(j) * W.hat)
				lambda[j] <- as.vector(T2/(2*T1))
				H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
				H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
				Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
				ifelse(invmethod == "eigen", Var.Y.inv <- solve(Var.Y), Var.Y.inv <- chol2inv(chol(Var.Y)) )
			}	
		}

		### Estimating psi using EM	
		w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		W.hat <- Var.Y.inv + tcrossprod(w.hat)
		T3 <- crossprod(Y-alpha) + sum(H.mat.lamsq * W.hat) - 2*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat))
		psi <- as.vector(N/T3)
	
		### New value of log-likelihood
		Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)	
		ifelse(invmethod == "eigen", Var.Y.inv <- solve(Var.Y), Var.Y.inv <- chol2inv(chol(Var.Y)) )
		log.lik1 <- dmvnorm(Y, mean=rep(alpha,N), sigma=Var.Y, log=T)
		
		### Report
		check <- i %% report.int
		if(log.lik1 < log.lik0) warning(paste("Log-likelihood decreased at iteration", i), call.=F)
		if(!is.na(check) && check==0 && !silent) cat("\nIteration", paste0(i, ":"), "\t Log-likelihood =", round(log.lik1, 4), " ")  
		if(!silent && report.int != maxit) setTxtProgressBar(pb, i)		
		if(i %% report.int*10 == 0 && !silent) pb <- txtProgressBar(min=i, max=report.int*10+i, style=1, char=".") 
			#reset progress bar
	}
	
	if(!silent && report.int != maxit) close(pb)
	converged <- !(abs(log.lik0 - log.lik1) > stop.crit)
	if(!silent && converged) cat("EM complete.\n", "\nNumber of iterations =", i, "\n")
	else if(!silent) cat("EM NOT CONVERGED!\n", "\nNumber of iterations =", i, "\n")
	if(!silent) cat("Log-likelihood = ", log.lik1, "\n")
	
	list(alpha=alpha, lambda=lambda, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, H.matsq=H.matsq, H.mat.lam=H.mat.lam, VarY=Var.Y, w.hat=w.hat, kernel=kernel, whichPearson=whichkernel, converged=converged, stop.crit=stop.crit, one.lam=one.lam, parsm=parsm, interactions=interactions, q=q)
}