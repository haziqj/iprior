fisher <- function(alpha, psi, lambda, P.matsq, S.mat, H.mat.lam, Var.Y.inv) {
	# Calculates iprior parameters' standard errors from the inverse observed
	# Fisher matrix. Used as helper function in summary().
  N <- nrow(Var.Y.inv)
	q <- length(lambda)
	F.mat <- NULL
	for (i in 1:q) {
		F.mat[[i]] <- Var.Y.inv %*% (psi * (2 * lambda[i] * P.matsq[[i]] +
		                                    S.mat[[i]]))
	}
	F.mat[[q + 1]] <- diag(1 / psi, N) - (2 / psi ^ 2) * Var.Y.inv
	Fisher <- matrix(0, nrow = q + 2, ncol = q + 2)
	for (i in 1:(q + 1)) {
    for (j in 1:(q + 1)) {
			Fisher[i + 1, j + 1] <- sum(F.mat[[i]] * F.mat[[j]]) / 2
		}
	}
	InverseFisher <- solve(Fisher[-1, -1])
	se <- sqrt(c(1 / sum(Var.Y.inv), diag(InverseFisher)))
	se
}
