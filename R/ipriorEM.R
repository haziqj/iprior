###
### The iprior EM
###

ipriorEM <- function(ipriorKernel, maxit=10, stop.crit=1e-7, report.int=1, silent=F, lambda.init=NULL, psi.init=NULL, clean=F, paramprogress=F){

	list2env(ipriorKernel, environment())
	list2env(BlockBstuff, environment())
	list2env(model, environment())
	environment(BlockB) <- environment()

	if(report.int == 0)	report.int <- maxit

	### Initialise parameters
	alpha <- as.numeric(mean(Y))
	if (is.null(psi.init)) psi <- abs(rnorm(1)) else psi <- psi.init
	if (is.null(lambda.init)) lambda <- abs(rnorm(q, sd=0.1))
	else{
		if(length(lambda.init) != q) stop(paste("Incorrect dimension of lambda initial values. vector of length", q, "required."), call.=F)
		else lambda <- lambda.init
	}
	lambda.fn <- function(){
		lambda <<- lambda[1:q]
		if(parsm && no.int > 0){ for(j in 1:no.int) lambda <<- c(lambda, lambda[intr[1,j]]*lambda[intr[2,j]]) }
	}
	lambda.fn()

	### Store results
	res.loglik <- matrix(NA, nr=maxit+1, nc=3)	#loglik, predloglik, delta
	res.param <- matrix(NA, nr=maxit+1, nc=2+q)
	rownames(res.loglik) <- paste0("Iteration ", 0:maxit, ":"); colnames(res.loglik) <- c("Log-lik.", "Pred.log-l.", "Delta_i,i-1")
	rownames(res.param) <- paste0("Iteration ", 0:maxit, ":")
	if(q == 1) colnames(res.param) <- c("(Intercept)", "lambda", "psi")
	else colnames(res.param) <- c("(Intercept)", paste0("lambda", 1:q), "psi")

	### Linear solver and inverse
	linsolvinv <- function(b=NULL){
		if(is.null(b)) a <- FastVdiag(V, 1/{u+s}) #a C++ alternative
		else a <- V %*% ( diag(1/(u+s)) %*% (t(V) %*% b) )
		a
	}

	### Calculating H.mat.lam
	if (q==1) {
		H.mat.lam.fn <- function() H.mat.lam <<- lambda[1] * P.mat[[1]]
	}
	else {
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

	list(alpha=alpha, lambda=lambda[1:q], psi=psi, log.lik=log.lik1, no.iter=i, P.matsq=P.matsq, S.mat=S.mat, H.mat.lam=H.mat.lam, Var.Y.inv=Var.Y.inv, w.hat=w.hat, converged=converged, res.loglik=res.loglik, res.param=res.param)
}
