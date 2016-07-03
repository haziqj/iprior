fisher <- function(object) {
	# Calculates iprior parameters (lambda, psi) standard errors from the inverse
	# observed Fisher matrix. Used as helper function in summary().
  lambda <- object$lambda
  psi <- object$psi
  Psql <- object$Psql
  Sl <- object$Sl
  VarY.inv <- object$VarY.inv
  N <- nrow(VarY.inv)
	l <- length(lambda)

	if (is.null(Psql)) {
    # Fitted using regular EM.
    Fisher <- -optimHess(c(lambda, psi), logLik, object = object,
                        control = list(fnscale = -1))
	} else {
	  # Fitted using closed-form EM.
	  dVarY <- NULL
	  for (i in 1:l) {
	    dVarY[[i]] <- VarY.inv %*% (psi * (2 * lambda[i] * Psql[[i]] + Sl[[i]]))
	  }
	  dVarY[[l + 1]] <- diag(1 / psi, N) - (2 / psi ^ 2) * VarY.inv
	  Fisher <- matrix(0, nrow = l + 1, ncol = l + 1)
	  for (i in 1:(l + 1)) {
	    for (j in 1:(l + 1)) {
	      Fisher[i, j] <- sum(dVarY[[i]] * dVarY[[j]]) / 2
	    }
	  }
	}
	Inverse.Fisher <- solve(Fisher)
	if (any(diag(Inverse.Fisher) < 0)) {
	  warning("NaNs S.E. produced due to negative inverse Hessian\nHas the EM converged?", call. = FALSE)
	}
	suppressWarnings(se <- sqrt(c(1 / sum(VarY.inv), diag(Inverse.Fisher))))
	se
}
