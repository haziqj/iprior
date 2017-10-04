iprior_em_closed <- function(mod, maxit = 500, stop.crit = 1e-5, silent = FALSE,
                             lambda0 = NULL, psi0 = NULL) {
  # Declare all variables and functions to be used into environment ------------
  iprior.env <- environment()
  list2env(mod, iprior.env)
  list2env(BlockBStuff, iprior.env)
  environment(BlockB) <- iprior.env
  environment(em_loop_logical) <- iprior.env
  maxit <- max(1, maxit)  # cannot have maxit <= 0

  # Initialise -----------------------------------------------------------------
  if (is.null(lambda0)) {
    lambda0 <- rnorm(p)
    if (p == 1) lambda0 <- abs(lambda0)
  }
  if (is.null(psi0)) psi0 <- abs(rnorm(1))
  lambda <- lambda0
  psi <- psi0
  niter <- 0
  loglik <- rep(NA, maxit)
  Hl <- expand_Hl_and_lambda(Hl, lambda, intr, intr.3plus)$Hl

  # The EM loop ----------------------------------------------------------------
  if (!silent) pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()

  while (em_loop_logical()) {
    # Block A --------------------------------------------------------------------
    Hlam <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
    list2env(eigen_Hlam(Hlam), environment())
    z <- psi * u ^ 2 + 1 / psi  # eigenvalues of Vy

    # Block C --------------------------------------------------------------------
    zinv.Vt <- t(V) / z
    Vy.inv.y <- as.numeric(crossprod(y, V) %*% zinv.Vt)
    w <- psi * Hlam %*% Vy.inv.y
    W <- V %*% zinv.Vt + tcrossprod(w)

    # Update lambda ------------------------------------------------------------
    for (k in seq_len(p)) {
      lambda <- expand_Hl_and_lambda(lambda, lambda, intr, NULL)$lambda
      BlockB(k)  # Updates Pl, Psql, and Sl
      T1 <- sum(Psql[[k]] * W)
      T2 <- crossprod(y, Pl[[k]]) %*% w - sum(Sl[[k]] * W) / 2
      lambda[k] <- as.numeric(T2 / T1)
    }

    # Update psi ---------------------------------------------------------------
    Hlamsq <- V %*% (t(V) * u ^ 2)
    T3 <- crossprod(y) + sum(Hlamsq * W) - 2 * crossprod(y, Hlam %*% w)
    psi <- sqrt(max(0, as.numeric(sum(diag(W)) / T3)))

    # Calculate log-likelihood ---------------------------------------------------
    logdet <- sum(log(z))
    loglik[niter + 1] <- -n / 2 * log(2 * pi) - logdet / 2 - crossprod(y, Vy.inv.y) / 2

    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, niter)
  }

  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate standard errors --------------------------------------------------
  # >>> needs adjusting depending on est.lambda / est.psi = T/F <<<
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

  # Clean up and close ---------------------------------------------------------
  convergence <- niter == maxit
  param <- kernel_to_param(kernels, lambda[1:p])
  theta <- param_to_theta(param, estl, logpsi = log(psi))$theta
  param.full <- theta_to_collapsed_param(theta, mod)

  if (!silent) {
    close(pb)
    if (convergence) cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }

  list(theta = theta, param.full = param.full,
       loglik = as.numeric(na.omit(loglik)),
       se = se, niter = niter, w = as.numeric(w), start.time = start.time,
       end.time = end.time, time = time.taken, convergence = convergence,
       message = NULL)

}

em_loop_logical <- function() {
  # Helper function to determine when to stop the while loop for the EM
  # algorithm.
  # If niter == 0 return TRUE because must complete 1 iteration.
  # If niter == 1 then stop if maxit == 1, otherwise continue (nothing to compare).
  # If niter > 1 then just check whether maxit reached or stop.crit reached.
  ll.diff <- loglik[niter] - loglik[niter - 1]
  crit1 <- (niter != maxit)
  crit2 <- (abs(ll.diff) > stop.crit)
  if (niter == 0) {
    return(TRUE)
  } else if (niter == 1) {
    return(crit1)
  } else {
    if (ll.diff < 0) {
      warning(paste0("Log-likelihood decreased at iteration ", niter),
              call. = FALSE)
    }
    return(crit1 & crit2)
  }
}
