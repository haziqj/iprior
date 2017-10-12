iprior_em_reg <- function(mod, maxit = 500, stop.crit = 1e-5, silent = FALSE,
                          theta0 = NULL, psi.reg = FALSE) {
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
    if (is.ipriorKernel_nys(mod)) {
      Hlam <- get_Hlam(mod, theta, get_Xl.nys(mod))
      list2env(eigen_Hlam_nys(Hlam), environment())
    } else {
      Hlam <- get_Hlam(mod, theta)
      list2env(eigen_Hlam(Hlam), environment())
    }
    z <- psi * u ^ 2 + 1 / psi  # eigenvalues of Vy

    # Block C ------------------------------------------------------------------
    zinv.Vt <- t(V) / z
    Vy.inv.y <- as.numeric(crossprod(y, V) %*% zinv.Vt)
    w <- psi * (V %*% crossprod(V * u, Vy.inv.y))
    W <- V %*% zinv.Vt + tcrossprod(w)

    # Update parameters other than psi -----------------------------------------
    if (!isTRUE(psi.reg)) psi <- NULL
    res.optim <- optim(theta, QEstep, psi = psi, object = mod, w = w, W = W,
                       method = "L-BFGS", control = list(trace = FALSE))
    theta <- res.optim$par

    # Update psi ---------------------------------------------------------------
    if (isTRUE(psi.reg)) {
      psi <- theta_to_psi(theta, mod)
    } else {
      Hlamsq <- V %*% (t(V) * u ^ 2)
      Hlam.w <- V %*% crossprod(V * u, w)
      T3 <- crossprod(y) + sum(Hlamsq * W) - 2 * crossprod(y, Hlam.w)
      psi <- sqrt(max(0, as.numeric(sum(diag(W)) / T3)))
      theta[ grep("psi", names(mod$thetal$theta))] <- log(psi)
    }

    # Calculate log-likelihood ---------------------------------------------------
    logdet <- sum(log(z))
    loglik[niter + 1] <- -n / 2 * log(2 * pi) - logdet / 2 - crossprod(y, Vy.inv.y) / 2

    niter <- niter + 1
    if (!silent) setTxtProgressBar(pb, niter)
  }

  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate standard errors --------------------------------------------------
  tmp <- optimHess(theta, loglik_iprior, object = mod)
  tmp <- eigenCpp(-tmp)
  u <- tmp$val + 1e-9
  V <- tmp$vec
  Fi.inv <- V %*% t(V) / u
  se <- sqrt(diag(Fi.inv))
  se <- convert_se(se, theta, mod)  # delta method to convert to parameter s.e.

  # Clean up and close ---------------------------------------------------------
  convergence <- niter != maxit
  param.full <- theta_to_collapsed_param(theta, mod)

  if (!silent) {
    close(pb)
    if (convergence) cat("Converged after", niter, "iterations.\n")
    else cat("Convergence criterion not met.\n")
  }

  list(theta = theta, param.full = param.full,
       loglik = as.numeric(na.omit(loglik)),
       se = se, niter = niter, w = as.numeric(w), start.time = start.time,
       end.time = end.time, time = time.taken,
       convergence = as.numeric(!convergence), message = NULL)
}

QEstep <- function(theta, psi = NULL, object, w, W) {
  # Q(theta) = psi  * sum(y ^ 2) + tr(Vy %*% W) - 2 * psi * crossprod(y,
  # Hlam %*% w)
  if (is.null(psi)) psi <- theta_to_psi(theta, object)
  if (is.ipriorKernel_nys(mod)) {
    Hlam <- get_Hlam(object, theta, get_Xl.nys(mod))
    list2env(eigen_Hlam_nys(Hlam), environment())
  } else {
    Hlam <- get_Hlam(object, theta)
    list2env(eigen_Hlam(Hlam), environment())
  }
  Hlamsq <- V %*% (t(V) * u ^ 2)
  Hlam.w <- V %*% crossprod(V * u, w)
  Vy <- psi * Hlamsq + diag(1 / psi, object$n)
  res <- psi * sum(object$y ^ 2) + sum(Vy * W) -
    2 * psi * crossprod(object$y, Hlam.w)
  as.numeric(res)
}
