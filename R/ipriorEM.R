###
### The iprior EM
###

ipriorEM <- function(x, y, whichkernel, interactions, one.lam, parsm, kernel, gamfbm, maxit, stop.crit, report.int, silent, alpha.init, lambda.init, psi.init, invmethod, clean, paramprogress){	
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
		no.int <- sum(Tmpo==2)	#counts the number of 2nd order interactions
		if(parsm) q <- p		#parsimonious interactions
		else q <- p + no.int	#non-parsimonious interactions
	}
	else{ #no interactions
		no.int <- 0
		q <- p
	} 			
	if(one.lam){ #single lambda
		if(q == 1) message("Option one.lam=T used with a single covariate anyway.")
		q <- 1
		#if(parsm) message("Parsimonious interactions ignored because one.lam=T.")
	} 		
	
	### Initialise parameters
	alpha <- as.numeric(mean(Y))
	if(is.null(lambda.init)) lambda <- abs(rnorm(q, sd=0.1))
	else{
		if(length(lambda.init) != q) stop(paste("Incorrect dimension of lambda initial values. vector of length", q, "required."), call.=F)
		else lambda <- lambda.init
	}
	psi <- psi.init
	lambda.fn <- function(){
		lambda <<- lambda[1:q]
		if(parsm && no.int > 0){ for(j in 1:no.int) lambda <<- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]]) }
	}
	lambda.fn()

	### Store results
	res.loglik <- matrix(NA, nr=maxit+1, nc=3)	#loglik, predloglik, delta
	res.param <- matrix(NA, nr=maxit+1, nc=2+q)
	rownames(res.loglik) <- paste0("Iteration ", 0:maxit, ":"); colnames(res.loglik) <- c("Log-lik.", "Pred.log-l.", "Delta_i,i-1")
	rownames(res.param) <- paste0("Iteration ", 0:maxit, ":")
	if(q == 1) colnames(res.param) <- c("(Intercept)", "lambda", "psi")	
	else colnames(res.param) <- c("(Intercept)", paste0("lambda", 1:q), "psi")
	
	### Setup the kernel matrix
	H.mat <- NULL
	for(j in 1:p){
		if(whichkernel[j]) H.mat[[j]] <- fn.H1(X[,j])					#Pearson
		else{
			if(kernel=="FBM") H.mat[[j]] <- fn.H3a(X[,j], gamma=gamfbm)	#FBM
			else H.mat[[j]] <- fn.H2a(X[,j])							#Canonical
		}
	}
	if(!is.null(interactions)){ #add in interactions, if any
		for(j in 1:no.int) H.mat[[p+j]] <- H.mat[[ Tmpf[1,j] ]] * H.mat[[ Tmpf[2,j] ]]
	}

	### Indexer function
	indx.fn <- function(k){	#indices for H.mat2
		ind.int1 <- Tmpf[1,]==k; ind.int2 <- Tmpf[2,]==k	#locating var/kernel matrix
		ind.int <- which(ind.int1 | ind.int2)				#of interactions (out of 1:no.int)
		k.int <- ind.int+p	#which kernel matrix has interactions involves k 
		k.int.lam <- c(Tmpf[1,][ind.int2], Tmpf[2,][ind.int1])	#which lambdas has interaction with k	
		nok <- (1:p)[-k]	#all variables excluding k
		k.noint <- which(!(ind.int1 | ind.int2)) + p	#the opposite of k.int
	
		find.H2 <- function(z){ #this function finds position of H2
			x <- z[1]; y <- z[2]
			which( (ind1==x & ind2==y) | (ind2==x & ind1==y) )
		}
		
		### P.mat %*% R.mat + R.mat %*% P.mat indices
		za <- which((ind1 %in% k & ind2 %in% nok) | (ind2 %in% k & ind1 %in% nok))
		grid.PR <- expand.grid(k.int, nok)
		zb <- which(	(ind1 %in% grid.PR[,1] & ind2 %in% grid.PR[,2]) |
						(ind2 %in% grid.PR[,1] & ind1 %in% grid.PR[,2])
		)
		grid.PR.lam <- expand.grid(k.int.lam, nok)
	
		### P.mat %*% U.mat + U.mat %*% P.mat indices	
		grid.PU1 <- expand.grid(k, k.noint)
		zc <- which(	(ind1 %in% grid.PU1[,1] & ind2 %in% grid.PU1[,2]) |
						(ind2 %in% grid.PU1[,1] & ind1 %in% grid.PU1[,2])
		)
		grid.PU2 <- expand.grid(k.int, k.noint)
		zd <- apply(grid.PU2, 1, find.H2)
		grid.PU.lam <- expand.grid(k.int.lam, k.noint)
		
		### P.mat %*% P.mat indices
		grid.Psq <- t(combn(c(k, k.int), 2))
		ze <- apply(grid.Psq, 1, find.H2)
		grid.Psq.lam <- NULL
		if(length(k.int.lam) > 0) grid.Psq.lam <- t(combn(c(0, k.int.lam), 2))
		
		list(	k.int=k.int, k.int.lam=k.int.lam, 
				PRU=c(za,zc,zb,zd), 
				PRU.lam1=c(	rep(0, length(nok)+length(k.noint)), 
							grid.PR.lam[,1], 
							grid.PU.lam[,1]	),
				PRU.lam2=c(nok, k.noint, grid.PR.lam[,2], grid.PU.lam[,2]),
				Psq=c(k, k.int), Psq.lam=k.int.lam,
				P2=ze, P2.lam1=grid.Psq.lam[,1], P2.lam2=grid.Psq.lam[,2]
		)
	}
	
	### Linear solver and inverse
	linsolvinv <- function(b=NULL){
		if(is.null(b)) a <- FastVdiag2(V,1/{u+s}) #a C++ alternative
		else a <- V %*% ( diag(1/(u+s)) %*% (t(V) %*% b) )
		a
	}
	
	### Calculating H.mat.lam
	if(one.lam | q==1){ 
		H.mat.lam.fn <- function() H.mat.lam <<- lambda[1] * P.mat[[1]]
	}
	else{
		H.mat.lam.fn <- function(){
			H.mat.lam <<- Reduce('+', mapply('*', H.mat[1:(p+no.int)], lambda[1:(p+no.int)], SIMPLIFY=F))
		}
	}

	### Block A update function
	BlockA <- function(){
		lambda.fn()
		H.mat.lam.fn()
		A <- H.mat.lam
		s <<- 1/psi
		tmp <- EigenCpp(A) #a C++ alternative
		u <<- psi*tmp$val^2
		V <<- tmp$vec
		isVarYneg <<- F; isVarYneg <<- any(u + s < 0)	#checks if Var.Y is negative
	}
	
	### Log-likelihood function
	log.lik.fn <- function(){
		a <- linsolvinv(Y-alpha)
		logdet <- Re(sum(log( (u + s)[u + s > 0]) ))
		log.lik <- -N/2*log(2*pi) - 1/2*logdet - 1/2*crossprod(Y-alpha, a)
		as.numeric(log.lik)
	}
	
	### Block B update function
	H.mat2 <- H.matsq <- P.mat <- P.matsq <- R.mat <- U.mat <- S.mat <- ind <- NULL
	if(one.lam | q==1){
		P.mat <- Reduce('+', mapply('*', H.mat, 1, SIMPLIFY=F))
		P.matsq <- list(P.mat %*% P.mat)
		P.mat <- list(P.mat)
		BlockB <- function(k){}
		S.mat <- list(matrix(0, nr=N, nc=N))
	}
	else{
		### Prepare the indices (also required for indx.fn)
		z <- 1:(p+no.int)
		ind1 <- rep(z, times=(length(z)-1):0)
		ind2 <- unlist(lapply(2:length(z), function(x) c(NA,z)[-(0:x)]))
		
		### Cross-product terms of square kernel matrices
		for(j in 1:length(ind1)){	#this is a list of (p+no.int)C2
			H.mat2[[j]] <- H.mat[[ ind1[j] ]] %*% H.mat[[ ind2[j] ]] + 
							H.mat[[ ind2[j] ]] %*% H.mat[[ ind1[j] ]]
		}
		
		if(!is.null(interactions) && parsm){ #CASE: parsimonious interactions only
			for(k in z){ #these do not depend on lambda, so setup once
				H.matsq[[k]] <- FastSquare(H.mat[[k]])
				if(k <= p) ind[[k]] <- indx.fn(k) 			
			}
			BlockB <- function(k){								
				indB <- ind[[k]]
				P.mat[[k]] <<- Reduce('+', mapply('*', H.mat[c(k,indB$k.int)], c(1,lambda[indB$k.int.lam]), SIMPLIFY=F))
				P.matsq[[k]] <<- Reduce('+', mapply('*', H.matsq[indB$Psq], 
									c(1, lambda[indB$Psq.lam]^2), SIMPLIFY=F))
				if(!is.null(indB$P2.lam1)) P.matsq[[k]] <<- P.matsq[[k]] +
											Reduce('+', mapply('*', H.mat2[indB$P2], 
												c(rep(1, sum(indB$P2.lam1==0)), lambda[indB$P2.lam1])*lambda[indB$P2.lam2], 
												SIMPLIFY=F))
				S.mat[[k]] <<- Reduce('+', mapply('*', H.mat2[indB$PRU], 
								c(rep(1, sum(indB$PRU.lam1==0)), lambda[indB$PRU.lam1])*lambda[indB$PRU.lam2],
								SIMPLIFY=F))												
			}
		}
		else{ #CASE: multiple lambda, no interactions, or non-parsimonious interactions
			for(k in 1:q){	#these do not depend on lambda, so setup once
				P.mat[[k]] <- H.mat[[k]]
				P.matsq[[k]] <- FastSquare(P.mat[[k]])			
			}
			BlockB <- function(k){	
				ind <- which(ind1==k | ind2==k)
				S.mat[[k]] <<- Reduce('+', mapply('*', H.mat2[ind], lambda[-k], SIMPLIFY=F))
			}
		}
	} 	
		
	### Block C update function
	BlockC <- function(){
		Var.Y.inv <<- linsolvinv()
		w.hat <<- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		W.hat <<- Var.Y.inv + tcrossprod(w.hat)		
	}
	
	### Checks and begin iterations
	i <- 0; check.naught <- 0
	H.mat.lam <- isVarYneg <- s <- u <- V <- Var.Y.inv <- w.hat <- W.hat <- 0	#reserve variables
	BlockA()
	log.lik0 <- log.lik.fn()
	log.lik1 <- log.lik0 + 2*stop.crit
	res.loglik[1,1] <- log.lik0
	res.param[1,] <- c(alpha, lambda[1:q], psi)
	if(!silent){
		if(clean) cat(format(paste0("Iteration " , 0, ":"), width=16, just="left"), "Log-likelihood = ", ipriorEMprettyLoglik(log.lik0), " ", sep="" )
		else{
			head.tab <- format(" ", width=16, just="right")
			if(paramprogress) head.tab <- c(head.tab, format(c(colnames(res.loglik), colnames(res.param)[-1]), width=11, just="right"))
			else head.tab <- c(head.tab, format(c(colnames(res.loglik)), width=11, just="right"))
			cat(head.tab, "\n")		#prints the table headers
			if(paramprogress) ipriorEMprettyIter(c(res.loglik[1,], res.param[1,-1]), 0)
			else ipriorEMprettyIter(res.loglik[1,], 0)
		}
	}
	if(!silent) pb <- txtProgressBar(min=0, max=report.int*10, style=1, char=".") #progress bar
	if(isVarYneg) warning(paste("Variance of Y is not positive definite at iteration", i), call.=F)
	
	### The EM algorithm
	while((i != maxit) && (abs(log.lik0 - log.lik1) > stop.crit)){
		i <- i + 1
		log.lik0 <- log.lik1
		
		### Update for lambda 
		BlockC() #Var.Y.inv through linear solver, and calculation of w.hat and W.hat
		for(k in 1:q){
			BlockB(k)
			T1 <- sum(P.matsq[[k]] * W.hat)
			T2 <- 2*crossprod(Y-alpha, crossprod(P.mat[[k]], w.hat)) - sum(S.mat[[k]] * W.hat)
			lambda[k] <- as.vector(T2/(2*T1))
		}

		### Update for psi
		H.mat.lamsq <- FastSquare(H.mat.lam) #a C++ alternative
		T3 <- crossprod(Y-alpha) + sum(H.mat.lamsq * W.hat) - 2*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat))	
		psi <- sqrt(max(0, as.numeric(sum(diag(W.hat))/T3)))

		### Estimating alpha
		# tmp.alpha <- crossprod(x0, Var.Y.inv)
		# alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))

		### New value of log-likelihood		
		BlockA() #H.mat.lam, and eigendecomposition
		log.lik1 <- log.lik.fn()

		### Storage
		dloglik <- log.lik1-log.lik0		
		dloglikold <- res.loglik[i,3]
		
		a <- ifelse((0 < dloglik) & (dloglik < dloglikold), dloglik/dloglikold, 0)
		predloglik <- log.lik1 - dloglik + dloglik/(1-a)
		res.loglik[i+1,] <- c(log.lik1, predloglik, dloglik)	
		res.param[i+1,] <- c(alpha, lambda[1:q], psi)
		
		### Report and conclusion
		check.naught <- max(0, i %% report.int)
		if(log.lik1 < log.lik0) warning(paste("Log-likelihood decreased at iteration", i), call.=F)
		if(!is.na(check.naught) && check.naught==0 && !silent){
			if(clean) cat("\n", format(paste0("Iteration " , i, ":"), width=16, just="left"), "Log-likelihood = ", ipriorEMprettyLoglik(log.lik1), " ", sep="")
			else{
				cat("\n")
				if(paramprogress) ipriorEMprettyIter(c(res.loglik[i+1,], res.param[i+1,-1]), i)
				else ipriorEMprettyIter(res.loglik[i+1,], i) 
			}
		}
		if(!silent) setTxtProgressBar(pb, i)		
		if(i %% report.int*10 == 0 && !silent) pb <- txtProgressBar(min=i, max=report.int*10+i, style=1, char=".")	#reset progress bar
		if(isVarYneg) warning(paste("Variance of Y is not positive definite at iteration", i), call.=F)
	}
	
	### Final report
	if(!silent && check.naught!=0){
		if(clean) cat("\n", format(paste0("Iteration " , i, ":"), width=16, just="left"), "Log-likelihood = ", ipriorEMprettyLoglik(log.lik1), " ", sep="")
		else{
			cat("\n")
			if(paramprogress) ipriorEMprettyIter(c(res.loglik[i+1,], res.param[i+1,-1]), i)
			else ipriorEMprettyIter(res.loglik[i+1,], i) 
		}
	}
	
	res.loglik <- res.loglik[1:(i+1),]; res.param <- res.param[1:(i+1),]
	
	if(!silent) close(pb)
	converged <- !(abs(log.lik0 - log.lik1) > stop.crit)
	if(!silent && converged) cat("EM complete.\n")#, "\nNumber of iterations =", i, "\n")
	else if(!silent) cat("EM NOT CONVERGED!\n")#, "\nNumber of iterations =", i, "\n")	
	
	list(alpha=alpha, lambda=lambda[1:q], psi=psi, log.lik=log.lik1, no.iter=i, P.matsq=P.matsq, S.mat=S.mat, H.mat.lam=H.mat.lam, Var.Y.inv=Var.Y.inv, w.hat=w.hat, kernel=kernel, whichPearson=whichkernel, converged=converged, stop.crit=stop.crit, one.lam=one.lam, parsm=parsm, interactions=interactions, q=q, p=p, res.loglik=res.loglik, res.param=res.param)
}