###
### The iprior EM
###

ipriorEM.old <- function(x, y, whichkernel, interactions, one.lam, parsm, kernel, gamfbm, maxit, stop.crit, report.int, silent, alpha.init, lambda.init, psi.init, invmethod, clean, paramprogress, force.regular){
	### Library packages
	require(Matrix, quietly=T)			#to create diagonal matrices
	require(MASS, quietly=T)			#to sample from MVN dist.
	#require(mvtnorm, quietly=T)
	require(numDeriv, quietly=T)
	require(mvnfast, quietly=T)
	
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
	
	### Store results
	res.loglik <- matrix(NA, nr=maxit+1, nc=3)	#loglik, predloglik, delta
	res.param <- matrix(NA, nr=maxit+1, nc=2+q)
	rownames(res.loglik) <- paste0("Iteration ", 0:maxit, ":"); colnames(res.loglik) <- c("Log-lik.", "Pred.log-l.", "Delta_i,i-1")
	rownames(res.param) <- paste0("Iteration ", 0:maxit, ":")
	if(q == 1) colnames(res.param) <- c("(Intercept)", "lambda", "psi")	
	else colnames(res.param) <- c("(Intercept)", paste0("lambda", 1:q), "psi")
	
	### Define the kernel matrix
	H.mat <- NULL; H.matsq <- NULL
	for(j in 1:p){
		if(whichkernel[j]) H.mat[[j]] <- fn.H1(X[,j])					#Pearson
		else{
			if(kernel=="FBM") H.mat[[j]] <- fn.H3(X[,j], gamma=gamfbm)	#FBM
			else H.mat[[j]] <- fn.H2a(X[,j])							#Canonical
		}
		if(q > 1) H.matsq[[j]] <- H.mat[[j]] %*% H.mat[[j]]
	}
	if(!is.null(interactions)){
		for(j in 1:no.int){
			H.mat[[p+j]] <- H.mat[[ Tmpf[1,j] ]] * H.mat[[ Tmpf[2,j] ]]
			H.matsq[[p+j]] <- H.mat[[p+j]] %*% H.mat[[p+j]]
		}
	}
	
	if(q == 1){	#just to save some time, if q=1 then no need to loop and define H2 and J.mat
		H.mat <- Reduce('+', mapply('*', H.mat, 1, SIMPLIFY=F))
		H.matsq <- list(H.mat %*% H.mat)
		H.mat <- list(H.mat)
		H.mat2 <- list(matrix(0, nr=N, nc=N))
		J.mat <- function(k) matrix(0, nr=N, nc=N)
		ind1 <- 1; ind2 <- 1
	}
	else{
		w <- 1:(p+no.int)
		ind1 <- rep(w, times=(length(w)-1):0)
		ind2 <- unlist(lapply(2:length(w), function(x) c(NA,w)[-(0:x)]))
		H.mat2 <- NULL
		for(j in 1:length(ind1)){
			H.mat2[[j]] <- H.mat[[ ind1[j] ]] %*% H.mat[[ ind2[j] ]] + H.mat[[ ind2[j] ]] %*% H.mat[[ ind1[j] ]]
		}
		J.mat <- function(k){	#always for non-parsimonious method q = p + no.int
			ind <- which(ind1==k | ind2==k)
			tmp <- Reduce('+', mapply('*', H.mat2[ind], lambda[-k], SIMPLIFY=F))
			tmp
		}
	}
	
	### Initialise the EM
	if((!is.null(interactions) && parsm) | force.regular){
		if(no.int > 0){for(j in 1:no.int) lambda <- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]])}
		### This is the E-step for lambda which needs to be maximised in the parsimonious method
		Q.lam <- function(lambda){
			if(no.int > 0){for(j in 1:no.int) lambda <- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]])}
			H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
			H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
			Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N)
			res <- psi/2*crossprod(Y-alpha) - psi*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat)) + 1/2*sum(Var.Y*W.hat)
			res
		}
	}
	
	vary.inv <- function(b=diag(N)){	#function to calculate inverse of Var.Y multiplied by b
		A <- psi*H.mat.lamsq
		s <- 1/psi
		tmp <- eigen(A)
		isVarYneg <<- F; isVarYneg <<- any(tmp$val + s < 0)
		Q <-  tmp$vectors
		x <- Q %*% ( diag(1/(tmp$val + s)) %*% (t(Q) %*% b) )
		list(x=Re(x), logdet=Re(sum(log((tmp$val + s)[tmp$val + s > 0]))))
	}
	
	fn.loglik <- function(){
		tmp <- vary.inv(b=matrix(Y-alpha, nc=1))
		log.lik <- -N/2*log(2*pi) - 1/2*tmp$logdet - 1/2*matrix(Y-alpha, nr=1) %*% tmp$x
		as.numeric(log.lik)
	}
	
	H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))
	H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
	lambda <- lambda[1:q]
	ifelse(invmethod == "eigen", Var.Y.inv <- vary.inv()$x, 
		{Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N); Var.Y.inv <- chol2inv(chol(Var.Y))} )		
	log.lik0 <- fn.loglik()

	log.lik1 <- log.lik0 + 2*stop.crit
	res.loglik[1,1] <- log.lik0
	res.param[1,] <- c(alpha, lambda, psi)
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
	i <- 0; check.naught <- 0
	
	if(isVarYneg) warning(paste("Variance of Y is not positive definite at iteration", i), call.=F)
	
	while((i != maxit) && (abs(log.lik0 - log.lik1) > stop.crit)){
	
		i <- i + 1
		log.lik0 <- log.lik1
				
		### Estimating alpha
		tmp.alpha <- crossprod(x0, Var.Y.inv)
		alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))
		
		### Estimating lambda using EM
		if((!is.null(interactions) && parsm) | force.regular){ #parsimonious method
			w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
			W.hat <- Var.Y.inv + tcrossprod(w.hat)
			lambda <- optim(lambda, fn=Q.lam, method="Nelder-Mead")$par
			if(no.int > 0){for(j in 1:no.int) lambda <- c(lambda, lambda[Tmpf[1,j]]*lambda[Tmpf[2,j]])}
			H.mat.lam <- Reduce('+', mapply('*', H.mat, lambda, SIMPLIFY=F))		
			H.mat.lamsq <- Reduce('+', mapply('*', H.matsq, lambda^2, SIMPLIFY=F)) + Reduce('+', mapply('*', H.mat2, lambda[ind1]*lambda[ind2], SIMPLIFY=F))
			lambda <- lambda[1:q]
			ifelse(invmethod == "eigen", Var.Y.inv <- vary.inv()$x, 
				{Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N); Var.Y.inv <- chol2inv(chol(Var.Y))} )		
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
			ifelse(invmethod == "eigen", Var.Y.inv <- vary.inv()$x, 
				{Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N); Var.Y.inv <- chol2inv(chol(Var.Y))} )		
			}	
		}

		### Estimating psi using EM	
		w.hat <- psi*H.mat.lam %*% (Var.Y.inv %*% matrix(Y - alpha, nc=1))
		W.hat <- Var.Y.inv + tcrossprod(w.hat)
		T3 <- crossprod(Y-alpha) + sum(H.mat.lamsq * W.hat) - 2*crossprod(Y-alpha, crossprod(H.mat.lam, w.hat))
		# psi <- as.vector(N/T3)
		# T4 <- sum(W.hat); psi <- as.vector(T4/N)	#The other sufficient statistic.		
		T5 <- sum(diag(W.hat))	#trace of W.hat
		psi <- sqrt(max(0,as.numeric(T5/T3)))
		ifelse(invmethod == "eigen", Var.Y.inv <- vary.inv()$x, 
			{Var.Y <- psi*H.mat.lamsq + (1/psi)*diag(N); Var.Y.inv <- chol2inv(chol(Var.Y))} )		
	
		### New value of log-likelihood		
		log.lik1 <- fn.loglik()
		
		### Storage
		dloglik <- log.lik1-log.lik0
		dloglikold <- res.loglik[i,3]
		a <- ifelse((0 < dloglik) & (dloglik < dloglikold), dloglik/dloglikold, 0)
		predloglik <- log.lik1 - dloglik + dloglik/(1-a)
		res.loglik[i+1,] <- c(log.lik1, predloglik, dloglik)	
		res.param[i+1,] <- c(alpha, lambda, psi)
		
		### Report
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
		if(i %% report.int*10 == 0 && !silent) pb <- txtProgressBar(min=i, max=report.int*10+i, style=1, char=".") 
			#reset progress bar
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
	#if(!silent) cat("Log-likelihood = ", log.lik1, "\n")
	
	list(alpha=alpha, lambda=lambda, psi=psi, log.lik=log.lik1, no.iter=i, H.mat=H.mat, H.matsq=H.matsq, H.mat.lam=H.mat.lam, VarY=psi*H.mat.lamsq + 1/psi*diag(N), w.hat=w.hat, kernel=kernel, whichPearson=whichkernel, converged=converged, stop.crit=stop.crit, one.lam=one.lam, parsm=parsm, interactions=interactions, q=q, p=p, res.loglik=res.loglik, res.param=res.param)
}