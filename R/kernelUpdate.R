###
### Update for ipriorKernel objects
###

update.ipriorKernel <- function (object, add=list(), delete=NULL, model=list()) {
	x <- object$x
	y <- object$Y
	mod <- object$model
	names(mod)[3] <- "interactions"
	mod[names(model)] <- model
	
	if (!is.null(delete)) {
		if(is.character(delete)) delete <- match(delete, names(x))
		x <- x[-delete]
		model$xname <- model$xname[-delete]
	}
	
	if (length(add) != 0) {
		if(!is.list(add)) add <- list(add)
		x <- c(x, add)
	}

	ipriorKernel <- kernL(y=y, x, model=mod)
	assign(deparse(substitute(object)), ipriorKernel, envir=parent.frame())
}

### DEPRECATED
# # update.iprior <- function(mod, ...){
	# newcall <- match.call()
	# cl <- mod$fullcall
	
	# if(is.null(mod$formula)){	#model fitted with x and y
		# newarg <- list(x=mod$xval, y=mod$yval)
		# mx <- charmatch(c("x"), names(cl), 0L)
		# my <- charmatch(c("y"), names(cl), 0L)
		# cl <- cl[-c(mx, my)]
	# }
	# else{
		# mdata <- charmatch(c("data"), names(cl), 0L)
		# mformula <- charmatch(c("formula"), names(cl), 0L)
		# newarg <- list(formula=mod$formula, data=cl[[mdata]])
		# cl <- cl[-c(mdata, mformula)]
	# }
	
	# #set up the new argument for iprior
	# k <- length(newarg)
	# if(length(cl) > 1){
		# for(i in 2:length(cl)){
			# newarg[[k+i-1]] <- cl[[i]]
			# names(newarg)[k+i-1] <- names(cl)[i]
		# }		
	# }

	# #put in the parameters
	# newarg$alpha <- mod$alpha
	# newarg$lambda <- mod$lambda
	# newarg$psi <- mod$psi
	# newarg$progress <- "lite"

	# #here we match the user's input and update the arguments if any	
	# ind.newarg <- charmatch(names(newcall), names(newarg), 0L); ind.newarg[1:2] <- -1 #picks out arguments in newarg that needs updating
	# ind.newcall1 <- which(ind.newarg > 0)		#picks out aditional arguments
	# ind.newcall2 <- which(ind.newarg == 0)	
	# ind.newarg <- ind.newarg[ind.newarg > 0]	
	# if(sum(ind.newarg) > 0){	#only do if there are things to match
		# for(i in 1:length(ind.newarg)) newarg[[ind.newarg[i]]] <- newcall[[ind.newcall1[i]]]
	# }

	# #following that, what remains is any other arguments the user has specified
	# k <- length(newarg); 
	# if(sum(ind.newcall2) > 0){	#only do if there are other new user arguments
		# for(i in 1:length(ind.newcall2)){
			# newarg[[k+i]] <- newcall[[ind.newcall2[i]]]
			# names(newarg)[k+i] <- names(newcall)[ind.newcall2[i]]
		# }
	# }


	# #check if parsimonious option was touched
	# parsm.ind <- which(!is.na(charmatch(names(newarg), "parsm")))
	# if(!all(is.na(parsm.ind))){
		# names(newarg)[parsm.ind] <- "parsm"
		# if(newarg$parsm == "F") tmp.log <- F 
		# else tmp.log <- T
		# if(mod$parsm != tmp.log){
			# newarg$lambda <- NULL
			# message("Lambda was reset because parsimonious option changed.")
		# } 
	# }
	
	# #remove progress option from call for aestheticness
	# mprogress <- charmatch("progress", names(newarg))
	
	# #aesthetic newarg
	# newarg.aesth <- newarg
	# newarg.aesth <- newarg.aesth[-mprogress]	
	
	# #specifically update the formula here
	# if(!is.null(mod$formula)){
		# if(as.formula(mod$formula) != as.formula(newarg$formula)){	#only do if formula was updated by user
			# newformula <- update.formula(mod$formula, newarg$formula, data=cl[[mdata]])
			# newarg$formula <- newformula	
			# #if there was a new update to formula, then makes sense to specify new starting values?
			# # newarg$alpha <- rnorm(1)
			# # newarg$lambda <- NULL
			# # newarg$psi <- 10
		# }	
	# }
	# else{
		# checkx <- all(dim(mod$xval) == dim(newarg$x))	#first check dimensions the same
		# if(checkx) checkx <- all(mod$xval != newarg$x)	#if yes, then check individual elements
		# checky <- all(dim(mod$yval) == dim(newarg$y))	#first check dimensions the same
		# if(checky) checky <- all(mod$yval != newarg$y)	#if yes, then check individual elements
		# if(checkx) newarg.aesth$x <- as.name("new.x") else newarg.aesth$x <- as.name("x")
		# if(checky) newarg.aesth$y <- as.name("new.y") else newarg.aesth$y <- as.name("y")
	# }
	
	# cat("Updating", deparse(substitute(mod)), "with the following call:", "\n\n")
	# print(as.call(c(as.name("iprior"), newarg.aesth)))
	# cat("\n")
	# tmp <- do.call("iprior", newarg)
	# assign(deparse(substitute(mod)), tmp, envir = .GlobalEnv)
# }