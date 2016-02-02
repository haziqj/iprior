###
### iprior update
###

update.iprior <- function(mod, ...){
	newcall <- match.call()
	cl <- mod$fullcall
	
	if(is.null(mod$formula)){	#model fitted with x and y
		newarg <- list(x=mod$xval, y=mod$yval)
	}
	else{
		mdata <- charmatch(c("data"), names(cl), 0L)
		mformula <- charmatch(c("formula"), names(cl), 0L)
		newarg <- list(formula=mod$formula, data=cl[[mdata]])
		cl <- cl[-c(mdata, mformula)]
	}

	#set up the new argument for iprior
	k <- length(newarg)
	for(i in 2:length(cl)){
		newarg[[k+i-1]] <- cl[[i]]
		names(newarg)[k+i-1] <- names(cl)[i]
	}
	
	#put in the parameters
	newarg$alpha <- mod$alpha
	newarg$lambda <- mod$lambda
	newarg$psi <- mod$psi
	newarg$progress <- "lite"

	#specifically update the formula here
	mformula2 <- charmatch("formula", names(newcall), 0L)
	if(sum(mformula2) > 0){	#only do ifthere is update to formula
		newformula <- update.formula(mod$formula, newcall[[mformula2]])
		newarg$formula <- newformula
		newcall <- newcall[-mformula2]
		#if there was a new update to formula, then makes sense to specify new starting values?
		# newarg$alpha <- rnorm(1)
		# newarg$lambda <- NULL
		# newarg$psi <- 10
	}

	#here we match the user's input and update the arguments if any	
	mcommon <- charmatch(names(newcall), names(newarg), 0L)
	mcommon[1:2] <- -1
	ind1 <- which(mcommon > 0)
	ind2 <- which(mcommon == 0)
	if(sum(ind1) > 0){	#only do if there are things to match
		newcall <- newcall[c(1, ind1, ind2)]
		mcommon <- mcommon[mcommon>0]
		for(i in 1:length(mcommon)) newarg[[ mcommon[i] ]] <- newcall[[i+1]]
	}

	#following that, what remains is any other arguments the user has specified
	k <- length(newarg)
	if(sum(ind2) > 0){	#only do if there are other new user arguments
		for(i in 1:length(ind2)){
			newarg[[k+i]] <- newcall[[1+length(ind1)+i]]
			names(newarg)[k+i] <- names(newcall)[1+length(ind1)+i]
		}
	}
	
	#check if parsimonious option was touched
	if(charmatch("parsm", names(newarg)) > 0){
		if(newarg$parsm == "F") tmp.log <- F 
		else tmp.log <- T
		if(mod$parsm != tmp.log){
			newarg$lambda <- NULL
			message("Lambda was reset because parsimonious option changed.")
		} 
	}
	
	#just for aesthetics, remove progress option from call
	mprogress <- charmatch("progress", names(newarg)); print(mprogress)

	cat("Updating", deparse(substitute(mod)), "with the following call:", "\n\n")
	print(as.call(c(as.name("iprior"), newarg[-mprogress])))
	cat("\n")
	do.call("iprior", newarg)
}
