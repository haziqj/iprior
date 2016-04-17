###
### The generic function is a standard R function with a special body, usually containing 
### only a call to UseMethod
###

iprior <- function(formula, data, one.lam, parsm, progress=c("lite", "none", "full", "predloglik"), ...) UseMethod("iprior")

### The default method
iprior.default <- function(x, y, interactions=NULL, parsm=T, one.lam=F, kernel=c("Canonical", "FBM"), gamfbm=NULL, maxit=50000, stop.crit=1e-7, report.int=100, alpha=rnorm(1), lambda=NULL, psi=10, invmethod=c("eigen", "chol"), progress=c("lite", "none", "full", "predloglik"), ...){
	kernel <- match.arg(kernel)
	invmethod <- match.arg(invmethod)
	progress <- match.arg(progress)
	if(progress == "lite"){ clean <- T; silent <- F; paramprogress <- F }
	if(progress == "none"){ clean <- T; silent <- T; paramprogress <- F }
	if(progress == "full"){ clean <- F; silent <- F; paramprogress <- T }
	if(progress == "predloglik"){ clean <- F; silent <- F; paramprogress <- F }
	ifelse(is.null(ncol(x)), Whichkernel <- is.factor(x), Whichkernel <- sapply(x, is.factor))
	x <- as.data.frame(x)
	y <- as.numeric(y)
	n <- length(y)
	
	est <- ipriorEM(x, y, whichkernel=Whichkernel, interactions=interactions, one.lam=one.lam, parsm=parsm, kernel=kernel, gamfbm=gamfbm, maxit=maxit, stop.crit=stop.crit, report.int=report.int, silent=silent, alpha.init=alpha, lambda.init=lambda, psi.init=psi, invmethod=invmethod, clean=clean, paramprogress=paramprogress)
	param <- c(est$alpha, est$lambda, est$psi)
	if(length(param) == 3) names(param) <- c("(Intercept)", "lambda", "psi")	
	else names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)), "psi")
	colnames(est$res.param) <- names(param)
	if(is.null(gamfbm)) gamma <- 0.5 else gamma <- gamfbm
	
	### Calculate fitted values
	Y.hat <- est$alpha + as.vector(crossprod(est$H.mat.lam, est$w.hat))
	est$fitted.values <- Y.hat
	est$residuals <- y-Y.hat
	
	### Changing the call to simply iprior
	cl <- match.call(); est$fullcall <- cl
    cl[[1L]] <- as.name("iprior")
    m <- match(c("x", "y", "one.lam", "parsm"), names(cl), 0L)
    cl <- cl[c(1L, m)]
	est$call <- cl
	
	### Other things to return	
	est$gamma <- gamma
	est$coefficients <- param
	est$yval <- y
	est$xval <- x
	est$sigma <- 1/sqrt(est$psi)
	est$T2 <- as.numeric(crossprod(est$w.hat) / est$psi)
	
	class(est) <- "iprior"
	est
}

print.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	if(x$kernel == "Canonical") CanOrFBM <- "Canonical" else CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ", x$gamma)
	kerneltypes <- c(CanOrFBM, "Pearson", paste(CanOrFBM, "& Pearson"))
	if(all(x$whichPearson)) cat("\nRKHS used:", kerneltypes[2])
	else{
		if(!all(x$whichPearson) && !any(x$whichPearson)) cat("\nRKHS used:", kerneltypes[1])
		else cat("\nRKHS used:", kerneltypes[3])
	} 
	if(x$q == 1) cat(", with a single scale parameter.\n")
	else cat(", with multiple scale parameters.\n")
	cat("\n")
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
		if(!object$one.lam){
			lambda.int <- lambda
			if(!is.null(object$interactions) && object$parsm){
				Tmpo <- object$interactions[[1]]
				Tmpf <- object$interactions[[2]]
				no.int <- sum(Tmpo==2)
				for(j in 1:no.int) lambda.int <- c(lambda.int, lambda[Tmpf[1, j]]*lambda[Tmpf[2, j]])
				H.mat.lam <- Reduce('+', mapply('*', object$H.mat, lambda.int, SIMPLIFY=F))
			} 
			else H.mat.lam <- Reduce('+', mapply('*', object$H.mat, lambda.int, SIMPLIFY=F))
		} 
		if(object$one.lam) H.mat.lam <- lambda * object$H.mat
		H.mat.lamsq <- H.mat.lam %*% H.mat.lam	
		Var.Y <- psi*H.mat.lamsq + (1/psi) * diag(n)
		loglik <- dmvn(y-alpha, rep(0,n), Var.Y, log=T)
		loglik
	}
	FisherInformation <- -hessian(ipriorloglik, coef(object))
	InverseFisher <- solve(FisherInformation)
	se <- sqrt(diag(InverseFisher))
	
	## Z values to compare against (standard) Normal distribution
	zval <- coef(object) / se
	
	## Create table for summary screen
	tab <- cbind(	Estimate=round(coef(object), digits=4),
					S.E.=round(se, digits=4),
					z=round(zval, digits=3),
					"P[|Z>z|]"=round(2*pnorm(-abs(zval)), digits=3) )
	if(object$q == 1){ #only rename rows when using multiple lambdas
		lamnames <- "lam"
		lamnames <- c("(Intercept)", paste(lamnames, attr(object$terms, "term.labels")[1:length(lamnames)], sep="."), "psi")
		rownames(tab) <- lamnames
	}
	else{
		lamnames <- paste0("lam", 1:(length(coef(object))-2))
		lamnames <- c("(Intercept)", paste(lamnames, attr(object$terms, "term.labels")[1:length(lamnames)], sep="."), "psi")
		rownames(tab) <- lamnames
	}
	#tab <- tab[-length(coef(object)),]	#removes the psi from the table
	
	res <- list(call=object$call, coefficients=tab, whichPearson=object$whichPearson, kernel=object$kernel, resid=object$residuals, log.lik=object$log.lik, no.iter=object$no.iter, converged=object$converged, stop.crit=object$stop.crit, one.lam=object$one.lam, T2=object$T2, q=object$q, gamma=object$gamma)
	class(res) <- "summary.iprior"
	res
}

print.summary.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	xPearson <- names(x$whichPearson)[x$whichPearson]
	xCanOrFBM <- names(x$whichPearson)[!x$whichPearson]
	if(x$kernel == "Canonical") CanOrFBM <- "Canonical" else CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ", x$gamma)
	printPearson <-	paste0("Pearson (", paste(xPearson, collapse=", "), ")")
	printCanOrFBM <- paste0(CanOrFBM, " (", paste(xCanOrFBM, collapse=", "), ")")
	cat("\n")
	cat("RKHS used:\n")
	if(!(length(xCanOrFBM) == 0)) cat(printCanOrFBM, "\n")
	if(!(length(xPearson) == 0)) cat(printPearson, "\n")
	if(x$q == 1) cat("with a single scale parameter.\n")
	else cat("with multiple scale parameters.\n")
	cat("\n")
	cat("Residuals:\n")
	print(summary(x$resid)[-4])
	cat("\n")
	tab <- x$coefficients
	psi.and.se <- tab[length(rownames(tab)),]
	sigma <- 1 / sqrt(psi.and.se[1])
	sesigma <- psi.and.se[2] * sigma^3 / 2
	tab <- tab[-length(rownames(tab)),]
	printCoefmat(tab, P.value=T, has.Pvalue=T)
	cat("\n")
	if(x$converged) cat("EM converged to within", x$stop.crit, "tolerance.")
	else cat("EM failed to converge.")
	cat(" No. of iterations:", x$no.iter)
	cat("\nStandard deviation of errors:", signif(sigma, digits=4), "with S.E.:", round(sesigma, digits=4))
	cat("\nT2 statistic:", signif(x$T2, digits=4), "on ??? degrees of freedom.")
	cat("\nLog-likelihood value:", x$log.lik, "\n")
	cat("\n")
}

## Formulas
iprior.formula <- function(formula, data=list(), ...){
	mf <- model.frame(formula=formula, data=data)
	tt <- terms(mf)
	Terms <- delete.response(tt)
	#attr(attr(mf, "terms"), "intercept") <- 0
	x <- model.frame(Terms, mf)
	y <- model.response(mf)
	
	## for interactions
	tmpo <- attr(tt, "order")
	tmpf <- attr(tt, "factors")
	tmpf2 <- as.matrix(tmpf[-1, tmpo==2])	#this obtains 2nd order interactions
	int2 <- apply(tmpf2, 2, function(x) which(x == 1))
	interactions <- list(tmpo=tmpo, tmpf=int2)
	
	ifelse(max(tmpo) > 1, 
		est <- iprior(x, y, interactions=interactions, ...),
		est <- iprior(x, y, ...)
	)
	
	## changing the call to simply iprior
	cl <- match.call(); est$fullcall <- cl
    cl[[1L]] <- as.name("iprior")
    m <- match(c("formula", "data", "one.lam", "parsm"), names(cl), 0L)
    cl <- cl[c(1L, m)]
	est$call <- cl
	est$formula <- formula
	names(est$fitted.values) <- row.names(mf)
	names(est$residuals) <- row.names(mf)
	est$terms <- tt
	est
}

## Prediction
predict.iprior <- function(object, newdata=NULL, ...){
	Y <- object$y; X <- object$x; p <- ncol(X)
	if(is.null(newdata)) ystar <- fitted(object)
	else{
		if(!is.null(object$formula)){
			## model has been fitted using formula interface
			tt <- terms(object)
			Terms <- delete.response(tt)
			xstar <- model.frame(Terms, newdata) #previously m
			# xstar <- model.matrix(Terms, m)
			# xstar <- model.matrix(object$formula, newdata)
			xcolnames <- colnames(xstar); xrownames <- rownames(xstar)
			# wheres.int <- (colnames(xstar) == "(Intercept)")
			# xstar <- matrix(xstar[,!wheres.int], nc=p)
			# colnames(xstar) <- xcolnames[!wheres.int]
		}
		else{
			xstar <- as.matrix(newdata)
		}
		
		## Define new kernel matrix
		if(!object$one.lam){ #for multiple lambdas
			H.mat <- NULL
			for(j in 1:p){
				if(is.factor(X[,j]))  H.mat[[j]] <- fn.H1(X[,j], xstar[,j]) 
				else H.mat[[j]] <- fn.H2a(X[,j], xstar[,j]) 			
			}
			if(!is.null(object$interactions) && object$parsm){ #for non-parsimonious interactions
				Tmpo <- object$interactions[[1]]
				Tmpf <- object$interactions[[2]]
				no.int <- sum(Tmpo==2)
				p1 <- p + no.int
				for(j in (p1-no.int+1):p1){
					H.mat[[j]] <- H.mat[[ Tmpf[1, j-p1+no.int] ]] * H.mat[[ Tmpf[2, j-p1+no.int] ]]
				}
				H.mat.lam <- Reduce('+', mapply('*', H.mat, object$lambda.int, SIMPLIFY=F))
			} 
			else H.mat.lam <- Reduce('+', mapply('*', H.mat, object$lambda, SIMPLIFY=F))
		}
		if(object$one.lam){ #for single lambdas
			H.mat <- 0
			for(j in 1:p) H.mat <- H.mat + fn.H2a(X[,j], xstar[,j]) 
			H.mat.lam <- object$lambda * H.mat
		}
		ystar <- as.vector(object$alpha + (H.mat.lam %*% object$w.hat))
		names(ystar) <- xrownames
	}
	ystar
}

### New feature soon
#merge.iprior <- function(x, ...) cat("hello")