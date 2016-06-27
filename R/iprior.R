###
### The generic function is a standard R function with a special body, usually containing 
### only a call to UseMethod
###

iprior <- function(y, ..., model=list(), control=list()) UseMethod("iprior")

### The default method
iprior.default <- function (y, ..., model=list(), control=list()) {
	
	con <- list( maxit=50000, stop.crit=1e-7, report.int=100, 
				 lambda=NULL, psi=abs(rnorm(1)),
				 progress="lite", silent=F )
	con_names <- names(con)
	con[(control_names <- names(control))] <- control
    if (length(noNms <- control_names[! control_names %in% con_names])) {
        warning("Unknown names in control options: ", paste(noNms, collapse = ", "), call.=F)
	}

	list2env(con, environment())
	silent_ <- silent
	.progress <- c("lite", "none", "full", "predloglik")
	progress <- match.arg(progress, .progress)
	if(progress == "lite"){ clean <- T; silent <- F; paramprogress <- F }
	if(progress == "none" | silent ){ clean <- T; silent <- T; paramprogress <- F }
	if(progress == "full"){ clean <- F; silent <- F; paramprogress <- T }
	if(progress == "predloglik"){ clean <- F; silent <- F; paramprogress <- F }
	cl <- match.call()
	
	### Accept objects of class 'ipriorKernel' and 'iprior'
	if (is(y, "ipriorKernel")) ipriorKernel <- y
	else if (is(y, "iprior")) {
		ipriorKernel <- y$ipriorKernel
		lambda <- y$lambda
		psi <- y$psi
		cl <- y$fullcall
	} 
	else ipriorKernel <- kernL(y, ..., model=model) #pass to kernel loader
	
	### Pass to iprior EM
	est <- ipriorEM(ipriorKernel, maxit, stop.crit, report.int, silent_, lambda, psi, clean, paramprogress)
	est$ipriorKernel <- ipriorKernel
	
	param <- c(est$alpha, est$lambda, est$psi)
	if(length(param) == 3) names(param) <- c("(Intercept)", "lambda", "psi")	
	else names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)), "psi")
	
	### Calculate fitted values
	if(maxit==0) Y.hat <- rep(est$alpha, nrow(est$H.mat.lam))
	else Y.hat <- est$alpha + as.vector(crossprod(est$H.mat.lam, est$w.hat))
	est$fitted.values <- Y.hat
	est$residuals <- ipriorKernel$Y - Y.hat
	names(est$fitted.values) <- names(ipriorKernel$Y)
	names(est$residuals) <- names(ipriorKernel$Y)
	
	### Changing the call to simply iprior
	est$fullcall <- cl
    cl[[1L]] <- as.name("iprior")
    m <- match(c("control"), names(cl), 0L)
    if(any(m > 0)) cl <- cl[-m]
	est$call <- cl
	
	### Other things to return	
	est$control <- con
	est$coefficients <- param
	est$sigma <- 1/sqrt(est$psi)
	est$T2 <- as.numeric(crossprod(est$w.hat) / est$psi)
	
	class(est) <- "iprior"
	if (is(y, "iprior")) assign(deparse(substitute(y)), est, envir=parent.frame())
	else est
}

print.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	if(x$ipriorKernel$model$kernel == "Canonical") CanOrFBM <- "Canonical" else CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ", x$ipriorKernel$model$Hurst)
	kerneltypes <- c(CanOrFBM, "Pearson", paste(CanOrFBM, "& Pearson"))
	if(all(x$ipriorKernel$whichPearson)) cat("\nRKHS used:", kerneltypes[2])
	else{
		if(!all(x$ipriorKernel$whichPearson) && !any(x$ipriorKernel$whichPearson)) cat("\nRKHS used:", kerneltypes[1])
		else cat("\nRKHS used:", kerneltypes[3])
	} 
	if(x$ipriorKernel$q == 1) cat(", with a single scale parameter.\n")
	else cat(", with multiple scale parameters.\n")
	cat("\n")
	cat("\nParameter estimates:\n")
	print(x$coefficients)
	cat("\n")
}

### The summary screen
summary.iprior <- function(object, ...){
	### Fisher information and standard errors
	se <- Fisher.fn(	alpha=object$alpha, psi=object$psi, lambda=object$lambda, 
						P.matsq=object$P.matsq, H.mat.lam=object$H.mat.lam, 
						S.mat=object$S.mat, 	Var.Y.inv=object$Var.Y.inv )
	
	### Z values to compare against (standard) Normal distribution
	zval <- coef(object) / se
	
	### Create table for summary screen
	tab <- cbind(	Estimate=round(coef(object), digits=4),
					S.E.=round(se, digits=4),
					z=round(zval, digits=3),
					"P[|Z>z|]"=round(2*pnorm(-abs(zval)), digits=3) )
	xname <- object$ipriorKernel$model$xname
	if(object$ipriorKernel$q == 1){ #only rename rows when using multiple lambdas
		lamnames <- c("(Intercept)", "lambda", "psi")
		rownames(tab) <- lamnames
	}
	else{
		lamnames <- paste0("lam", 1:(length(coef(object))-2))
		lamnames <- c("(Intercept)", paste(lamnames, xname[1:object$ipriorKernel$q], sep="."), "psi")
		rownames(tab) <- lamnames
	}
	tab <- tab[-length(coef(object)),]	#removes the psi from the table
	
	res <- list( call=object$call, coefficients=tab, whichPearson=object$ipriorKernel$whichPearson, 
				 kernel=object$ipriorKernel$model$kernel, resid=object$residuals, log.lik=object$log.lik, 
				 no.iter=object$no.iter, converged=object$converged, stop.crit=object$control$stop.crit, 
				 one.lam=object$ipriorKernel$model$one.lam, T2=object$T2, q=object$ipriorKernel$q, 
				 p=object$ipriorKernel$p, Hurst=object$ipriorKernel$model$Hurst, formula=object$formula,
				 psi.and.se=c(coef(object)[length(se)], se[length(se)]), xname=xname )
	class(res) <- "summary.iprior"
	res
}

print.summary.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
    x.names <- x$xname[1:x$q]
	xPearson <- x.names[x$whichPearson]
	xCanOrFBM <- x.names[!x$whichPearson]
	if(x$kernel == "Canonical") CanOrFBM <- "Canonical" else CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ", x$Hurst)
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
	psi.and.se <- x$psi.and.se
	sigma <- 1 / sqrt(psi.and.se[1])
	sesigma <- psi.and.se[2] * sigma^3 / 2
	# tab <- tab[-length(rownames(tab)),]
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

### Formulas
iprior.formula <- function(formula, data, model=list(), control=list()){
	
	### Pass to iprior default
	ipriorKernel <- kernL(formula, data, model=model)
	est <- iprior(ipriorKernel, control=control)	
	
	### Changing the call to simply iprior
	cl <- match.call(); est$fullcall <- cl
    cl[[1L]] <- as.name("iprior")
    m <- match(c("formula", "data"), names(cl), 0L)
    cl <- cl[c(1L, m)]
	est$call <- cl
	est$formula <- formula
	est$terms <- 
	
	class(est) <- "iprior"
	est
}

### Prediction
predict.iprior <- function (object, newdata=list(), ...) {
	list2env(object$ipriorKernel, environment())
	list2env(model, environment())
	
	if (length(newdata) == 0) ystar <- object$fitted
	else {
		if (!is.null(object$formula)) { #model has been fitted using formula interface
			mf <- model.frame(formula=object$formula, data=newdata)
			tt <- terms(mf)
			Terms <- delete.response(tt)
			xstar <- model.frame(Terms, newdata)
			xrownames <- rownames(xstar)
			xstar <- unlist(list(xstar), recursive=F)
		}
		else {
			xstar <- newdata
			xrownames <- rownames(do.call(cbind, newdata))
		}
		
		### Define new kernel matrix
		H.mat <- Hmat_list(x, kernel, whichPearson, intr, no.int, gamma, xstar)
		lambda <- object$lambda
		if (parsm && no.int > 0){ for (j in 1:no.int) lambda <- c(lambda, lambda[intr[1,j]]*lambda[intr[2,j]]) }
		H.mat.lam <- Reduce('+', mapply('*', H.mat[1:(p+no.int)], lambda[1:(p+no.int)], SIMPLIFY=F))
		
		### Calculate fitted values
		ystar <- as.vector(object$alpha + (H.mat.lam %*% object$w.hat))
		names(ystar) <- xrownames
	}
	ystar
}
