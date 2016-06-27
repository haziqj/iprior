###
### Kernel loader
###

kernL <- function (y, ..., model=list()) UseMethod("kernL")

kernL.default <- function (y, ..., model=list()) {
	x <- list(...)
	if (any(sapply(x, is.list))) x <- unlist(x, recursive=F)
	N <- length(y)
	p <- length(x)
	whichPearson <- unlist(lapply(x, is.factor))
	
	### Model options and checks
	mod <- list( kernel="Canonical", Hurst=0.5, interactions=NULL, parsm=T, one.lam=F, 
				 yname="y", xname=NULL)
	mod_names <- names(mod)
	mod[(model_names <- names(model))] <- model
    if (length(noNms <- model_names[!model_names %in% mod_names])) {
        warning("Unknown names in model options: ", paste(noNms, collapse = ", "), call.=F)	
	}
	.kernel <- c("Canonical", "FBM")
	mod$kernel <- match.arg(mod$kernel, .kernel)
	
	### Set up interactions, p and q
	names(mod)[3] <- "intr" #rename to something simpler
	if (!is.null(mod$intr)) { 
		if (!is.matrix(mod$intr)) { #not fitted using formula
			if (!is.character(mod$intr)) stop("Incorrect prescription of interactions.")
			mod$intr <- sapply(strsplit(mod$intr, ":"), as.numeric)
		}
		no.int <- ncol(mod$intr)
		if(mod$parsm) q <- p else q <- p + no.int
	}
	else { #no interactions
		no.int <- 0L
		q <- p
	}
	if (mod$one.lam) { #only relevant when fitted using formula
		if(q == 1) message("Option one.lam=T used with a single covariate anyway.")
		p <- q <- 1	
	}
	if (any(mod$intr > p | mod$intr < 1)) stop("Prescribed interactions out of bounds.")
	
	### Set up names
	if (is.null(mod$xname)) mod$xname <- names(x)
	else names(x) <- mod$xname[1:p]
	if (suppressWarnings(is.null(mod$xname) | any(names(x) == "") | any(is.na(names(x))))) {
		cl <- match.call()
		m <- match(c("y", "model", "control"), names(cl), 0L)
		xnamefromcall <- as.character(cl[-m])[-1]
		mod$xname <- xnamefromcall
	}
	suppressWarnings( here <- which((names(x) != "") & !is.na(names(x))) )
	mod$xname[here] <- names(x)[here]
	names(x) <- mod$xname[1:p]	
	
	
	### Set up list of H matrices
	H.mat <- Hmat_list(x, mod$kernel, whichPearson, mod$intr, no.int, mod$Hurst)
	names(H.mat) <- mod$xname[1:length(H.mat)]
	if (length(mod$xname) < length(H.mat)) {
		print(4)
		for (i in 1:ncol(mod$intr)) {
			mod$xname <- c(mod$xname, paste(mod$xname[mod$intr[1,i]], mod$xname[mod$intr[2,i]], sep=":"))
		}
		names(H.mat) <- mod$xname
	}
	
	### Set up progress bar
	# pb <- txtProgressBar(min=0, max=1, style=3)#, width=30)
	pb.count <- 0
	
	### Block B update function
	intr <- mod$intr
	environment(indx.fn) <- environment()
	H.mat2 <- H.matsq <- P.mat <- P.matsq <- S.mat <- ind <- ind1 <- ind2 <- NULL
	if (mod$one.lam | q==1) {
		P.mat <- Reduce('+', mapply('*', H.mat, 1, SIMPLIFY=F))
		P.matsq <- list(FastSquare(P.mat))
		if (mod$one.lam) mod$xname <- paste0("(", paste(names(H.mat), collapse=" + "), ")")
		H.mat <- P.mat <- list(P.mat)
		x <- list(as.data.frame(x))
		names(H.mat) <- names(x) <- mod$xname
		if (mod$one.lam) whichPearson <- F
		BlockB <- function(k){}
		S.mat <- list(matrix(0, nr=N, nc=N))
		# setTxtProgressBar(pb, 1)
	}
	else {
		### Prepare the indices (also required for indx.fn)
		z <- 1:(p+no.int)
		ind1 <- rep(z, times=(length(z)-1):0)
		ind2 <- unlist(lapply(2:length(z), function(x) c(NA,z)[-(0:x)]))
		# pb <- txtProgressBar(min=0, max=length(c(ind1,z)), style=3, title="test")#, width=30)
		
		### Cross-product terms of square kernel matrices
		for (j in 1:length(ind1)) {	#this is a list of (p+no.int)C2
			H.mat2[[j]] <- H.mat[[ ind1[j] ]] %*% H.mat[[ ind2[j] ]] + 
							H.mat[[ ind2[j] ]] %*% H.mat[[ ind1[j] ]]
			pb.count <- pb.count + 1
			# setTxtProgressBar(pb, pb.count)					
		}
		
		if (!is.null(intr) && mod$parsm) { #CASE: parsimonious interactions only
			for (k in z) { #these do not depend on lambda, so setup once
				H.matsq[[k]] <- FastSquare(H.mat[[k]])
				if(k <= p) ind[[k]] <- indx.fn(k)
				pb.count <- pb.count + 1
				# setTxtProgressBar(pb, pb.count)	 			
			}
			BlockB <- function (k) {								
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
			for (k in 1:q) {	#these do not depend on lambda, so setup once
				P.mat[[k]] <- H.mat[[k]]
				P.matsq[[k]] <- FastSquare(P.mat[[k]])	
				pb.count <- pb.count + 1
				# setTxtProgressBar(pb, pb.count)	 		
			}
			BlockB <- function (k) {	
				ind <- which(ind1==k | ind2==k)
				S.mat[[k]] <<- Reduce('+', mapply('*', H.mat2[ind], lambda[-k], SIMPLIFY=F))
			}
		}
	}
	# close(pb) 		
	BlockBstuff <- list(H.mat2=H.mat2, H.matsq=H.matsq, P.mat=P.mat, P.matsq=P.matsq, S.mat=S.mat, ind1=ind1, ind2=ind2, ind=ind, BlockB=BlockB)
	
	kernelLoaded <- list(Y=y, x=x, H.mat=H.mat, N=N, p=p, q=q, no.int=no.int, whichPearson=whichPearson, BlockBstuff=BlockBstuff, model=mod)
	class(kernelLoaded) <- "ipriorKernel"
	kernelLoaded
}

kernL.formula <- function (formula, data, model=list()) {
	mf <- model.frame(formula=formula, data=data)
	tt <- terms(mf)
	Terms <- delete.response(tt)
	x <- model.frame(Terms, mf)
	y <- model.response(mf)
	
	### For interactions
	interactions <- NULL
	tmpo <- attr(tt, "order")
	if (any(tmpo>2)) stop("iprior does not currently work with higher order interactions.")
	tmpf <- attr(tt, "factors")
	tmpf2 <- as.matrix(tmpf[-1, tmpo==2])	#this obtains 2nd order interactions
	int2 <- apply(tmpf2, 2, function(x) which(x == 1))
	if (any(tmpo==2)) interactions <- int2 

	kernelLoaded <- kernL( y=y, x, model=c( model, list(interactions=interactions,
										    yname=names(attr(tt, "dataClasses"))[1],
										    xname=attr(tt, "term.labels") ) ) )
	kernelLoaded
}

print.ipriorKernel <- function (x, ...) {
	cat("\n")
	# if (x$model$kernel == "Canonical") CanOrFBM <- "Canonical" else CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ", x$gamfbm)
	# kerneltypes <- c(CanOrFBM, "Pearson", paste(CanOrFBM, "& Pearson"))
	# if (all(x$whichPearson)) cat(kerneltypes[2], "RKHS loaded")
	# else {
		# if (!all(x$whichPearson) && !any(x$whichPearson)) cat(kerneltypes[1], "RKHS loaded")
		# else cat(kerneltypes[3], "RKHS loaded")
	# } 
	# if (x$q == 1 | x$model$one.lam) cat(", with a single scale parameter.\n")
	# else cat(", with", x$q, "scale parameters.\n")	
	cat("Sample size = ", x$N, "\n")
	cat("Number of scale parameters, p = ", x$q, "\n")
	cat("Number of interactions = ", x$no.int, "\n")
	cat("\nInfo on H matrix:\n\n")
	str(x$H.mat)
	cat("\n")
}

# ### Decide to write this in the future
# sameLambda <- function(term){
	# if(is.name(term) || !is.language(term)) return(term)	
	# if (term[[1]] == as.name("I")) return(term)
	# if (term[[1]] == as.name("^")) return(term)
	# if (term[[1]] == as.name("(")) return(term[[-1]])
	# stopifnot(is.call(term))
	# # if (length(term) == 2){
		# # cat("if4", as.character(term), "--fb--", as.character(term[[2]]), "\n")
		# # return(sameLambda(term[[2]]))
	# # }
	# c(sameLambda(term[[2]]), sameLambda(term[[3]]))
# }