###
### The generic function is a standard R function with a special body, usually containing 
### only a call to UseMethod
###

iprior <- function(formula, data, one.lam, parsm, progress=c("lite", "none", "full", "predloglik"), ...) UseMethod("iprior")

### The default method
iprior.default <- function(x, y, interactions=NULL, parsm=T, one.lam=F, kernel=c("Canonical", "FBM"), Hurst=NULL, maxit=50000, stop.crit=1e-7, report.int=100, alpha=rnorm(1), lambda=NULL, psi=10, invmethod=c("eigen", "chol"), progress=c("lite", "none", "full", "predloglik"), silent_=F, ...){
	kernel <- match.arg(kernel)
	invmethod <- match.arg(invmethod)
	progress <- match.arg(progress)
	if(progress == "lite"){ clean <- T; silent <- F; paramprogress <- F }
	if(progress == "none" | silent_ ){ clean <- T; silent <- T; paramprogress <- F }
	if(progress == "full"){ clean <- F; silent <- F; paramprogress <- T }
	if(progress == "predloglik"){ clean <- F; silent <- F; paramprogress <- F }
	x <- as.data.frame(x)
	ifelse(is.null(ncol(x)), Whichkernel <- is.factor(x), Whichkernel <- sapply(x, is.factor))
	y <- as.numeric(y)
	n <- length(y)
	
	est <- ipriorEM(x, y, whichkernel=Whichkernel, interactions=interactions, one.lam=one.lam, parsm=parsm, kernel=kernel, gamfbm=Hurst, maxit=maxit, stop.crit=stop.crit, report.int=report.int, silent=silent, alpha.init=alpha, lambda.init=lambda, psi.init=psi, invmethod=invmethod, clean=clean, paramprogress=paramprogress)
	param <- c(est$alpha, est$lambda, est$psi)
	if(length(param) == 3) names(param) <- c("(Intercept)", "lambda", "psi")	
	else names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)), "psi")
	# colnames(est$res.param) <- names(param)
	if(is.null(Hurst)) gamma <- 0.5 else gamma <- Hurst
	
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
	if(object$q == 1){ #only rename rows when using multiple lambdas
		lamnames <- c("(Intercept)", "lambda", "psi")
		rownames(tab) <- lamnames
	}
	else{
		lamnames <- paste0("lam", 1:(length(coef(object))-2))
		lamnames <- c("(Intercept)", paste(lamnames, attr(object$terms, "term.labels")[1:length(lamnames)], sep="."), "psi")
		rownames(tab) <- lamnames
	}
	tab <- tab[-length(coef(object)),]	#removes the psi from the table
	
	res <- list(call=object$call, coefficients=tab, whichPearson=object$whichPearson, kernel=object$kernel, resid=object$residuals, log.lik=object$log.lik, no.iter=object$no.iter, converged=object$converged, stop.crit=object$stop.crit, one.lam=object$one.lam, T2=object$T2, q=object$q, p=object$p, gamma=object$gamma, formula=object$formula, psi.and.se=c(coef(object)[length(se)], se[length(se)]))
	class(res) <- "summary.iprior"
	res
}

print.summary.iprior <- function(x, ...){
	cat("\nCall:\n")
	print(x$call)
	if(!is.null(x$formula)) x.names <- names(x$whichPearson)
	else x.names <- paste0("X", 1:x$p)
	xPearson <- x.names[x$whichPearson]
	xCanOrFBM <- x.names[!x$whichPearson]
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

## Formulas
iprior.formula <- function(formula, data=list(), ...){
	mf <- model.frame(formula=formula, data=data)
	tt <- terms(mf)
	Terms <- delete.response(tt)
	#attr(attr(mf, "terms"), "intercept") <- 0
	x <- model.frame(Terms, mf)
	y <- model.response(mf)
	
	### For interactions
	tmpo <- attr(tt, "order")
	if(any(tmpo>2)) stop("iprior does not currently work with higher order interactions.")
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
	est$yval <- y
	est$xval <- x
	est$yname <- names(attr(tt, "dataClasses"))[1]
	est
}

### Prediction
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
			# wheres.int <- (colnames(xstar) == "(Intercept)")
			# xstar <- matrix(xstar[,!wheres.int], nc=p)
			# colnames(xstar) <- xcolnames[!wheres.int]
		}
		else{
			xstar <- as.matrix(newdata)
		}
		xcolnames <- colnames(xstar); xrownames <- rownames(xstar)
		
		## Define new kernel matrix
		if(!object$one.lam){ #for multiple lambdas
			H.mat <- NULL
			for(j in 1:p){
				if(is.factor(X[,j]))  H.mat[[j]] <- fn.H1(X[,j], xstar[,j])								#Pearson
				else{
					if(object$kernel=="FBM") H.mat[[j]] <- fn.H3a(X[,j], xstar[,j], gamma=object$gam)	#FBM
					else H.mat[[j]] <- fn.H2a(X[,j], xstar[,j])											#Canonical
				} 			
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
			for(j in 1:p){
				if(is.factor(X[,j]))  H.mat <- H.mat + fn.H1(X[,j], xstar[,j])							#Pearson
				else{
					if(object$kernel=="FBM") H.mat <- H.mat + fn.H3a(X[,j], xstar[,j], gamma=object$gam)	#FBM
					else H.mat <- H.mat + fn.H2a(X[,j], xstar[,j])										#Canonical
				} 			
			H.mat.lam <- object$lambda * H.mat
			}
		}
		ystar <- as.vector(object$alpha + (H.mat.lam %*% object$w.hat))
		names(ystar) <- xrownames
	}
	ystar
}

### New feature soon
#merge.iprior <- function(x, ...) cat("hello")