################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2018  Haziq Jamil
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

iprior_em_closed <- function(mod, maxit = 500, stop.crit = 1e-5, silent = FALSE,
                             theta0 = NULL, lambda0 = NULL, psi0 = NULL,
                             mixed = FALSE, omega = 0) {
  # The closed-form EM algorithm method for estimating I-prior models. This only
  # works when there are no other hyperparameters to estimate other than lambda
  # and psi.
  #
  # Args: An ipriorKernel object (mod), number of maximum EM iterations (maxit),
  # the stopping criterion to determine convergence (stop.crit), an option to
  # suppress reporting (silent), the initial values (theta0) and a logical
  # argument used internally to determine whether this EM routine is part of the
  # iprior_mixed() routine or not.
  #
  # Returns: A list containing the optimised theta and parameters, y.hat and w,
  # loglik values, standard errors, number of iterations, time taken, and
  # convergence information.

  # Declare all variables and functions to be used into environment ------------
  iprior.env <- environment()
  list2env(mod, iprior.env)
  list2env(BlockBStuff, iprior.env)
  environment(BlockB) <- iprior.env
  environment(em_loop_logical) <- iprior.env
  maxit <- max(1, maxit)  # cannot have maxit <= 0

  # Initialise -----------------------------------------------------------------
  # if (is.null(lambda0)) {
  #   lambda0 <- rnorm(p)
  #   if (p == 1) lambda0 <- abs(lambda0)
  # }
  # if (is.null(psi0)) psi0 <- abs(rnorm(1))
  # lambda <- lambda0
  # psi <- psi0
  if (is.null(theta0)) theta0 <- rnorm(mod$thetal$n.theta)
  psi <- theta_to_psi(theta0, mod)
  lambda.new <- lambda <- theta_to_collapsed_param(theta0, mod)[seq_len(mod$p)]
  niter <- 0
  loglik <- rep(NA, maxit)
  Hl <- expand_Hl_and_lambda(Hl, rep(1, p), intr, intr.3plus)$Hl  # expand Hl

  # The EM loop ----------------------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (em_loop_logical()) {
    # Block A ------------------------------------------------------------------
    Hlam <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
    eigen_Hlam(Hlam, iprior.env)  # assign u and V to environment
    z <- psi * u ^ 2 + 1 / psi  # eigenvalues of Vy

    # Block C ------------------------------------------------------------------
    zinv.Vt <- t(V) / z
    Vy.inv.y <- as.numeric(crossprod(y, V) %*% zinv.Vt)
    w <- psi * Hlam %*% Vy.inv.y
    W <- V %*% zinv.Vt + tcrossprod(w)

    # Update lambda ------------------------------------------------------------
    for (k in seq_len(p)) {
      lambda <- expand_Hl_and_lambda(lambda[1:p], lambda[1:p], intr, NULL)$lambda
      BlockB(k)  # Updates Pl, Psql, and Sl
      T1 <- sum(Psql[[k]] * W)
      T2 <- crossprod(y, Pl[[k]]) %*% w - sum(Sl[[k]] * W) / 2
      lambda.new[k] <- as.numeric(T2 / T1)
    }
    lambda <- (1 + omega) * lambda.new - omega * lambda[1:p]

    # Update psi ---------------------------------------------------------------
    Hlamsq <- V %*% (t(V) * u ^ 2)
    T3 <- crossprod(y) + sum(Hlamsq * W) - 2 * crossprod(y, Hlam %*% w)
    psi.new <- sqrt(max(0, as.numeric(sum(diag(W)) / T3)))
    psi <- (1 + omega) * psi.new - omega * psi

    # Calculate log-likelihood ---------------------------------------------------
    logdet <- sum(log(z))
    loglik[niter + 1] <- -n / 2 * log(2 * pi) - logdet / 2 -
      crossprod(y, Vy.inv.y) / 2

    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, niter)
  }

  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate fitted values  ---------------------------------------------------
  y.hat <- get_y.hat(u, V, w)

  # Calculate standard errors --------------------------------------------------
  if (!isTRUE(mixed)) {
    Vy.inv <- V %*% zinv.Vt
    dVy <- NULL
    for (i in seq_len(p)) {
      dVy[[i]] <- Vy.inv %*% (psi * (2 * lambda[i] * Psql[[i]] + Sl[[i]]))
    }
    dVy[[p + 1]] <- diag(1 / psi, n) - (2 / psi ^ 2) * Vy.inv
    Fi <- matrix(0, nrow = p + 1, ncol = p + 1)
    for (i in seq_len(p + 1)) {
      for (j in seq_len(p + 1)) {
        Fi[i, j] <- sum(dVy[[i]] * dVy[[j]]) / 2
      }
    }
    se <- sqrt(diag(solve(Fi)))
  } else {
    se <- NULL
  }

  # Clean up and close ---------------------------------------------------------
  convergence <- niter == maxit
  param <- kernel_to_param(kernels, lambda[1:p])
  theta <- param_to_theta(param, estl, logpsi = log(psi))$theta
  param.full <- theta_to_collapsed_param(theta, mod)

  if (!silent) {
    close(pb)
    if (isTRUE(mixed)) cat("")
    else if (convergence) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  list(theta = theta, param.full = param.full, se = se,
       loglik = as.numeric(na.omit(loglik)), w = as.numeric(w), y.hat = y.hat,
       niter = niter,  start.time = start.time, end.time = end.time,
       time = time.taken, convergence = convergence, message = NULL)
}
