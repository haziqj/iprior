###
### The iprior EM
###

ipriorLoad <- function(x, whichkernel, one.lam, kernel, gamfbm){	
	X <- x
	N <- nrow(X)
	p <- ncol(X)
	
	### Setup the kernel matrix
	H.mat <- NULL
	for(j in 1:p){
		if(whichkernel[j]) H.mat[[j]] <- fn.H1(X[,j])					#Pearson
		else{
			if(kernel=="FBM") H.mat[[j]] <- fn.H3a(X[,j], gamma=gamfbm)	#FBM
			else H.mat[[j]] <- fn.H2a(X[,j])							#Canonical
		}
	}
	if(one.lam){ 
		H.mat <- Reduce('+', mapply('*', H.mat, 1, SIMPLIFY=F))
		H.mat <- list(H.mat)
	}
	
	list(H.mat=H.mat, p=p, whichPearson=whichkernel, one.lam=one.lam, kernel=kernel, gamfbm=gamfbm)
}