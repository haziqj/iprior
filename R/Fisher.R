################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2016  Haziq Jamil
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

fisher <- function(object) {
	# Calculates iprior parameters (lambda, psi) standard errors from the inverse
	# observed Fisher matrix. Used as helper function in summary().
  lambda <- object$lambda
  psi <- object$psi
  Psql <- object$Psql
  Sl <- object$Sl
  VarY.inv <- object$VarY.inv
  force.regEM <- object$control$force.regEM
  N <- nrow(VarY.inv)
	l <- length(lambda)

	if (is.null(Sl) | force.regEM) {
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
