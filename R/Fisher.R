###
### The Fisher information
###

Fisher.fn <- function(alpha, psi, lambda, P.matsq, S.mat, H.mat.lam, Var.Y.inv){
	N <- nrow(Var.Y.inv)
	q <- length(lambda)
	F.mat <- NULL
	for(i in 1:q){
		F.mat[[i]] <- Var.Y.inv %*% ( psi*(2*lambda[i]*P.matsq[[i]] + S.mat[[i]]) ) 
	}
	F.mat[[q+1]] <- diag(1/psi, N)
	Fisher <- matrix(0, nr=q+2, nc=q+2)
	# Fisher[1,1] <- sum(Var.Y.inv)
	for(i in 1:(q+1)){
		for(j in 1:(q+1)){
			Fisher[i+1,j+1] <- 1/2*sum(F.mat[[i]]*F.mat[[j]])
		}
	}
	InverseFisher <- solve(Fisher[-1,-1])
	se <- sqrt(c(1/sum(Var.Y.inv), diag(InverseFisher)))
	se
}