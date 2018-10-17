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

iprior_em_reg <- function(mod, maxit = 500, stop.crit = 1e-5, silent = FALSE,
                          theta0 = NULL, psi.reg = FALSE, mixed = FALSE,
                          omega = 0) {
  # The regular version of the EM algorithm for estimating I-prior models. This
  # maximises the M-step using the quasi-Newton L-BFGS algorithm.
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
  y <- mod$y
  n <- mod$n
  environment(em_loop_logical) <- environment()
  maxit <- max(1, maxit)  # cannot have maxit <= 0

  # Initialise -----------------------------------------------------------------
  if (is.null(theta0)) theta <- rnorm(mod$thetal$n.theta)
  else theta <- theta0
  psi <- theta_to_psi(theta, mod)
  niter <- 0
  loglik <- rep(NA, maxit)

  # The EM loop ----------------------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (em_loop_logical()) {
    # Block A ------------------------------------------------------------------
    Hlam <- get_Hlam(mod, theta)
    list2env(eigen_Hlam(Hlam), environment())
    z <- psi * u ^ 2 + 1 / psi  # eigenvalues of Vy

    # Block C ------------------------------------------------------------------
    zinv.Vt <- t(V) / z
    Vy.inv.y <- as.numeric(crossprod(y, V) %*% zinv.Vt)
    w <- psi * Hlam %*% Vy.inv.y
    W <- V %*% zinv.Vt + tcrossprod(w)

    # Update parameters other than psi -----------------------------------------
    if (!isTRUE(psi.reg)) psi.optim <- NULL
    else psi.optim <- psi
    res.optim <- optim(theta, QEstep, psi = psi.optim, object = mod, w = w, W = W,
                       method = "L-BFGS", control = list(trace = FALSE))
    theta.new <- res.optim$par
    theta <- (1 + omega) * theta.new - omega * theta

    # Update psi ---------------------------------------------------------------
    if (isTRUE(psi.reg)) {
      psi <- theta_to_psi(theta, mod)
    } else {
      Hlamsq <- V %*% (t(V) * u ^ 2)
      T3 <- crossprod(y) + sum(Hlamsq * W) - 2 * crossprod(y, Hlam %*% w)
      psi.new <- sqrt(max(0, as.numeric(sum(diag(W)) / T3)))
      psi <- (1 + omega) * psi.new - omega * psi
      theta[grep("psi", names(mod$thetal$theta))] <- log(psi)
    }

    # Calculate log-likelihood -------------------------------------------------
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
    tmp <- optimHess(theta, loglik_iprior, object = mod)
    tmp <- eigenCpp(-tmp)
    u <- tmp$val + 1e-9
    V <- tmp$vec
    Fi.inv <- V %*% t(V) / u
    se <- sqrt(diag(Fi.inv))
    se <- convert_se(se, theta, mod)  # delta method to convert to parameter s.e.
  } else {
    se <- NULL
  }

  # Clean up and close ---------------------------------------------------------
  convergence <- niter == maxit
  param.full <- theta_to_collapsed_param(theta, mod)
  names(theta) <- names(mod$thetal$theta)

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

QEstep <- function(theta, psi = NULL, object, w, W) {
  # This is the Q function (E-step) for the EM algorithm. It is defined as
  # Q(theta) = psi  * sum(y ^ 2) + tr(Vy %*% W) - 2 * psi * crossprod(y, Hlam
  # %*% w)
  #
  # Args: parameters theta and psi (optional), the ipriorKernel object, the
  # vector w and matrix W.
  #
  # Returns: Numeric.
  if (is.null(psi)) psi <- theta_to_psi(theta, object)
  Hlam <- get_Hlam(object, theta)
  Vy <- psi * fastSquare(Hlam) + diag(1 / psi, object$n)
  res <- psi * sum(object$y ^ 2) + sum(Vy * W) -
    2 * psi * crossprod(object$y, Hlam %*% w)
  as.numeric(res)
}
