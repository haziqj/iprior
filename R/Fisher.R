fisher <- function(alpha, psi, lambda, Psql, Sl, Hlam.mat, VarY.inv) {
	# Calculates iprior parameters' standard errors from the inverse observed
	# Fisher matrix. Used as helper function in summary().
  N <- nrow(VarY.inv)
	l <- length(lambda)
	F.mat <- NULL
	for (i in 1:l) {
		F.mat[[i]] <- VarY.inv %*% (psi * (2 * lambda[i] * Psql[[i]] +
		                                    Sl[[i]]))
	}
	F.mat[[l + 1]] <- diag(1 / psi, N) - (2 / psi ^ 2) * VarY.inv
	Fisher <- matrix(0, nrow = l + 2, ncol = l + 2)
	for (i in 1:(l + 1)) {
    for (j in 1:(l + 1)) {
			Fisher[i + 1, j + 1] <- sum(F.mat[[i]] * F.mat[[j]]) / 2
		}
	}
	InverseFisher <- solve(Fisher[-1, -1])
	se <- sqrt(c(1 / sum(VarY.inv), diag(InverseFisher)))
	se
}
