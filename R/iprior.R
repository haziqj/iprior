###
### The generic function is a standard R function with a special body, usually containing 
### only a call to UseMethod
###

iprior <- function(x, y, ...) UseMethod("iprior")

## The default method
iprior.default <- function(x, y, interactions=NULL, parsm=T, one.lam=F, maxit=50000, stop.crit=0.001, report.int=100, silent=F, ...){
	ifelse(is.null(ncol(x)), Whichkernel <- is.factor(x), Whichkernel <- sapply(x, is.factor))
	x <- as.data.frame(x)
	y <- as.numeric(y)
	n <- length(y)
	
	ifelse(!one.lam, { 
		if(!is.null(interactions) && parsm){
			est <- ipriorEM3(x, y, whichkernel=Whichkernel, interactions=interactions, maxit=maxit, stop.crit=stop.crit, report.int=report.int, silent=silent)
			param <- c(est$alpha, est$lambda, est$psi)
			names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)), "psi")
			H.mat.lam <- Reduce('+', mapply('*', est$H.mat, est$lambda.int, SIMPLIFY=F))
		}
		else{
			est <- ipriorEM2(x, y, whichkernel=Whichkernel, interactions=interactions, maxit=maxit, stop.crit=stop.crit, report.int=report.int, silent=silent)
			param <- c(est$alpha, est$lambda, est$psi)
			names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)), "psi")
			H.mat.lam <- Reduce('+', mapply('*', est$H.mat, est$lambda, SIMPLIFY=F))
		} }, {
			est <- ipriorEM1(x, y, whichkernel=Whichkernel, interactions=interactions, maxit=maxit, stop.crit=stop.crit, report.int=report.int, silent=silent)
			param <- c(est$alpha, est$lambda, est$psi)
			names(param) <- c("(Intercept)", "lambda", "psi")		
			H.mat.lam <- est$lambda * est$H.mat			
		}		
	)

	##calculate fitted values
	H.mat.lamsq <- H.mat.lam %*% H.mat.lam	
	Var.Y <- est$psi*H.mat.lamsq + (1/est$psi) * diag(n)
	Var.Y.inv <- solve(Var.Y)	
	w.hat <- est$psi*H.mat.lam %*% (Var.Y.inv %*% matrix(y - est$alpha, nc=1))
	W.hat <- Var.Y.inv + tcrossprod(w.hat)
	Y.hat <- est$alpha + as.vector(crossprod(H.mat.lam, w.hat))
	
	est$call <- match.call()
	est$coefficients <- param
	est$fitted.values <- Y.hat
	est$residuals <- y-Y.hat
	est$w.hat <- w.hat
	est$yval <- y
	est$xval <- x
	est$one.lam <- one.lam
	est$parsm <- parsm
	est$interactions <- interactions
	est$sigma <- 1/sqrt(est$psi)
	est$T2 <- as.numeric(crossprod(est$w.hat) / est$psi)
	
	class(est) <- "iprior"
	est
}

print.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	kerneltypes <- c("Canonical", "Pearson", "Canonical & Pearson")
	if(all(x$kernel)) cat("\nRKHS used:", kerneltypes[2])
	else{
		if(!all(x$kernel) && !any(x$kernel)) cat("\nRKHS used:", kerneltypes[1])
		else cat("\nRKHS used:", kerneltypes[3])
	} 
	#if(x$one.lam) cat(", with single scale parameter.\n")
	#else cat(", with multiple scale parameter.\n")
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
		loglik <- dmvnorm(y-alpha, rep(0,n), Var.Y, log=T)
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
	if(!object$one.lam){ #only rename rows when using multiple lambdas
		lamnames <- paste0("lam", 1:(length(coef(object))-2))
		lamnames <- c("(Intercept)", paste(lamnames, attr(object$terms, "term.labels")[1:length(lamnames)], sep="."), "psi")
		rownames(tab) <- lamnames
	}
	#tab <- tab[-length(coef(object)),]	#removes the psi from the table
	
	res <- list(call=object$call, coefficients=tab, kernel=object$kernel, resid=object$residuals, log.lik=object$log.lik, no.iter=object$no.iter, converged=object$converged, stop.crit=object$stop.crit, one.lam=object$one.lam, T2=object$T2)
	class(res) <- "summary.iprior"
	res
}

print.summary.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	xPearson <- names(x$kern)[x$kern]
	xCanonical <- names(x$kern)[!x$kern]
	printPearson <-	paste0("Pearson (", paste(xPearson, collapse=", "), ")")
	printCanonical <- paste0("Canonical (", paste(xCanonical, collapse=", "), ")")
	cat("\n")
	cat("RKHS used:\n")
	if(!(length(xCanonical) == 0)) cat(printCanonical, "\n")
	if(!(length(xPearson) == 0)) cat(printPearson, "\n")
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
	est$call <- match.call()
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