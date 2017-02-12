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

ipriorEM <- function(ipriorKernel, maxit = 10, stop.crit = 1e-7, report.int = 1,
                     silent = FALSE, alpha = NULL, lambda.init = NULL,
                     psi.init = NULL, clean = FALSE, paramprogress = FALSE,
                     force.regEM = FALSE, force.nlm = FALSE,
                     getHlam = FALSE, getVarY = FALSE){
  # This is the EM algorithm engine which estimates the I-prior model
  # parameters.
  #
  # Args: ipriorKernel Output from kernL() function. maxit The maximum number of
  # iterations. Defaults to 10 for debugging, but this is fed in from control
  # list. stop.crit The tolerance for the difference in log-likelihood value to
  # stop the EM. report.int The reporting interval for the EM. silent Logical,
  # if TRUE then no print report. lambda.init, psi.init Initial values for
  # lambda and psi. clean Logical, if FALSE then progress of log-likelihood
  # reported. paramprogress Logical, if TRUE then progress of parameters
  # reported. force.regEM Logical, for debugging of the regular EM routine.
  # getHlam and getVarY logical, used to save the Hlam/VarY.inv matrix and
  # return as part of Hlam()/varyinv() function

  # Declare all variables and functions to be used in this environment ---------
  ipriorEM.env <- environment()
	list2env(ipriorKernel, ipriorEM.env)
	list2env(BlockBstuff, ipriorEM.env)
	list2env(model, ipriorEM.env)
  environment(linSolvInv) <- environment(logLikEM) <- ipriorEM.env
	environment(BlockA) <- environment(BlockB) <- environment(BlockC) <- ipriorEM.env
  environment(lambdaExpand) <- environment(lambdaContract) <- ipriorEM.env
  environment(fisherNew) <- ipriorEM.env
  if (r > 0 | force.regEM | no.int.3plus > 0) {
    if (force.nlm) {
      environment(ipriorEMnlm) <- ipriorEM.env
      ipriorEMRoutine <- ipriorEMnlm
    } else if (l > 1) {
      environment(ipriorEMOptim2) <- ipriorEM.env
      ipriorEMRoutine <- ipriorEMOptim2
    } else {
      environment(ipriorEMOptim1) <- ipriorEM.env
      ipriorEMRoutine <- ipriorEMOptim1
    }
  } else {
    environment(ipriorEMClosedForm) <- ipriorEM.env
    ipriorEMRoutine <- ipriorEMClosedForm
  }

	# Initialise parameters ------------------------------------------------------
  if (is.null(alpha)) alpha <- as.numeric(mean(Y))
	if (is.null(psi.init)) psi <- abs(rnorm(1)) else psi <- psi.init
	if (is.null(lambda.init)) lambda <- abs(rnorm(l, sd = 0.1))
	else{
		if (length(lambda.init) != l) {
		  stop(paste("Incorrect dimension of lambda initial values. vector of
		             length", l, "required."), call. = FALSE)
		} else {
		  lambda <- lambda.init
		}
	}

	# Results storage, and for use in progress() ---------------------------------
	res.loglik <- matrix(NA, nrow = maxit + 1, ncol = 3)	# loglik, predlik, delta
	res.param <- matrix(NA, nrow = maxit + 1, ncol = 2 + l)
	rownames(res.loglik) <- paste0("Iteration ", 0:maxit, ":")
	colnames(res.loglik) <- c("Log-lik.", "Pred.log-l.", "Delta_i,i-1")
	rownames(res.param) <- paste0("Iteration ", 0:maxit, ":")
	if (l == 1) colnames(res.param) <- c("(Intercept)", "lambda", "psi")
	else colnames(res.param) <- c("(Intercept)", paste0("lambda", 1:l), "psi")

	# Function to calculate Hlam.mat ---------------------------------------------
	if (q == 1) {
		hlamFn <- function(x = lambda, env = ipriorEM.env) {
		  assign("Hlam.mat", x[1] * Pl[[1]], envir = env)
		}
	}
	else {
	  hlamFn <- function(x = lambda, env = ipriorEM.env) {
			assign("Hlam.mat", Reduce("+", mapply("*", Hl[1:q], x[1:q],
			                                      SIMPLIFY = FALSE)), envir = env)
		}
	}

	# Checks and begin iterations ------------------------------------------------
	if (report.int == 0)	report.int <- maxit
	i <- 0
	check.naught <- 0
	Hlam.mat <- is.VarYneg <- s <- u <- V <- VarY.inv <- w.hat <- W.hat <- 0
	BlockA()  # lambda expanded here
	log.lik0 <- logLikEM()
	log.lik1 <- log.lik0 + 2 * stop.crit
	res.loglik[1,1] <- log.lik0
	lambdaContract()  # for printing
	res.param[1,] <- c(alpha, lambda, psi)
	if (!silent) {
		if (clean) {
		  cat(format(paste0("Iteration " , 0, ":"), width = 16, just = "left"),
		      "Log-likelihood = ", ipriorEMprettyLoglik(log.lik0), " ", sep = "" )
		} else {
			head.tab <- format(" ", width = 16, just = "right")
			if (paramprogress) {
			  head.tab <- c(head.tab, format(c(colnames(res.loglik),
			                                   colnames(res.param)[-1]),
			                                 width = 11, just = "right"))
			}
			else {
			  head.tab <- c(head.tab, format(c(colnames(res.loglik)), width = 11,
			                                 just = "right"))
			}
			cat(head.tab, "\n")		#prints the table headers
			if (paramprogress) {
			  ipriorEMprettyIter(c(res.loglik[1,], res.param[1,-1]), 0)
			}
			else {
			  ipriorEMprettyIter(res.loglik[1,], 0)
			}
		}
	}
	if (!silent) {
	  pb <- txtProgressBar(min = 0, max = report.int * 10, style = 1, char = ".")
	}
	if (is.VarYneg) {
	  warning(paste("Variance of Y is not positive definite at iteration", i),
	          call. = FALSE)
	}

	# The EM algorithm routine ---------------------------------------------------
	while ((i != maxit) && (abs(log.lik0 - log.lik1) > stop.crit)) {
	  i <- i + 1
		log.lik0 <- log.lik1

    # Update for parameters lambda and psi -------------------------------------
		BlockC()  # obtains Var.Y.inv and updates w.hat and W.hat
		ipriorEMRoutine()

		# New value of log-likelihood ----------------------------------------------
		BlockA()  # performs Hlam.mat update and eigendecomposition
		log.lik1 <- logLikEM()

		# Storage ------------------------------------------------------------------
		dloglik <- log.lik1 - log.lik0
		dloglikold <- res.loglik[i,3]

		# Calculate predicted log-likelihood at this iteration ---------------------
		a <- ifelse((0 < dloglik) & (dloglik < dloglikold), dloglik/dloglikold, 0)
		predloglik <- log.lik1 - dloglik + dloglik/(1 - a)
		res.loglik[i + 1, ] <- c(log.lik1, predloglik, dloglik)
		lambdaContract()  # for printing
		res.param[i + 1, ] <- c(alpha, lambda, psi)

		# Report and conclusion ----------------------------------------------------
		check.naught <- max(0, i %% report.int)
		if (log.lik1 < log.lik0) {
		  warning(paste("Log-likelihood decreased at iteration", i), call. = FALSE)
		}
		if (!is.na(check.naught) && check.naught == 0 && !silent) {
			if (clean) {
			  cat("\n", format(paste0("Iteration " , i, ":"), width = 16,
			                   just = "left"), "Log-likelihood = ",
			      ipriorEMprettyLoglik(log.lik1), " ", sep = "")
			} else{
				cat("\n")
				if (paramprogress) {
				  ipriorEMprettyIter(c(res.loglik[i + 1, ], res.param[i + 1, -1]), i)
				}
				else ipriorEMprettyIter(res.loglik[i + 1,], i)
			}
		}
		if (!silent) setTxtProgressBar(pb, i)
		if (i %% (report.int * 10) == 0 && !silent) {
		  # Reset progress bar.
		  pb <- txtProgressBar(min = i, max = (report.int * 10) + i, style = 1,
		                       char = ".")
		}
		if (is.VarYneg) {
		  warning(paste("Variance of Y is not positive definite at iteration", i),
		          call. = FALSE)
		}
	}

	# Need the w.hat if maxit is zero --------------------------------------------
	if (maxit == 0L) {
	  VarY.inv <- linSolvInv()
	  tmp1 <-  psi * Hlam.mat
	  tmp2 <- VarY.inv %*% matrix(Y - alpha, ncol = 1)
	  w.hat <- psi * Hlam.mat %*% (VarY.inv %*% matrix(Y - alpha, ncol = 1))
	}

	# Do you just need Hlam or VarY.inv? -----------------------------------------
	if (getHlam) return(Hlam.mat)
	if (getVarY) return(VarY.inv)

	# # Calculate standard errors --------------------------------------------------
	# se <- NULL
	# if (!not.finalEM) {
	#   se <- fisherNew()  # makes use of the logLik() on the ipriorKernel
	# }

	# Calculate fitted value -----------------------------------------------------
	Y.hat <- alpha + as.vector(crossprod(Hlam.mat, w.hat))

	# Final report ---------------------------------------------------------------
	if (!silent && check.naught != 0L) {
		if (clean) {
		  cat("\n", format(paste0("Iteration " , i, ":"), width = 16, just = "left"),
		      "Log-likelihood = ", ipriorEMprettyLoglik(log.lik1), " ", sep = "")
		} else{
			cat("\n")
			if (paramprogress) {
			  ipriorEMprettyIter(c(res.loglik[i + 1, ], res.param[i + 1,-1]), i)
			}
			else ipriorEMprettyIter(res.loglik[i + 1, ], i)
		}
	}
	res.loglik <- res.loglik[1:(i + 1), ]
	res.param <- res.param[1:(i + 1), ]
	if (!silent) close(pb)
	converged <- !(abs(log.lik0 - log.lik1) > stop.crit)
	if (!silent && converged) {
	  cat("EM complete.\n")#, "\nNumber of iterations =", i, "\n")
	} else if (!silent) {
	  cat("EM NOT CONVERGED!\n")#, "\nNumber of iterations =", i, "\n")
	}


	list(alpha = alpha, lambda = lambda, psi = psi, log.lik = log.lik1,
	     no.iter = i, fitted.values = Y.hat,
	     # Psql = Psql, Sl = Sl, VarY.inv = VarY.inv, Hlam.mat = Hlam.mat,
	     w.hat = w.hat, converged = converged,
	     res.loglik = res.loglik, res.param = res.param)
}
