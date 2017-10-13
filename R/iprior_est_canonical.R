iprior_canonical <- function(mod, theta0 = NULL, control) {
  iprior.env <- environment()
  w <- loglik <- NULL
  start.time <- Sys.time()
  res <- optim(theta0, loglik_canonical, object = mod, env = iprior.env,
               trace = TRUE, get.w = TRUE, method = "L-BFGS", control = control,
               hessian = TRUE)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  tmp <- eigenCpp(-res$hessian)
  u <- tmp$val + 1e-9
  V <- tmp$vec
  Fi.inv <- V %*% t(V) / u
  se <- sqrt(diag(Fi.inv))
  se <- convert_se(se, res$par, mod)  # delta method to convert to parameter s.e.
  loglik <- as.numeric(na.omit(loglik))
  param.full <- theta_to_collapsed_param(res$par, mod)

  list(theta = res$par, param.full = param.full, loglik = loglik,
       se = se, niter = res$count[1], w = as.numeric(w), start.time = start.time,
       end.time = end.time, time = time.taken, convergence = res$convergence,
       message = res$message)

}

loglik_canonical <- function(theta, object, trace = FALSE, env = NULL,
                             get.w = FALSE) {
  psi <- theta_to_psi(theta, object)
  lambda <- theta_to_collapsed_param(theta, mod)[seq_len(object$p)]
  X <- matrix(unlist(object$Xl), nrow = object$n)
  XtX <- crossprod(X)
  Lambda <- diag(lambda)
  Vy.inv <- psi * (
    diag(1, object$n) - X %*% solve(solve(psi ^ 2 * Lambda %*% XtX %*% Lambda) + XtX, t(X))
  )
  Vy.inv.y <- Vy.inv %*% object$y
  logdet <- determinant(Vy.inv)$mod
  res <- as.numeric(
    -object$n / 2 * log(2 * pi) + logdet / 2 - crossprod(object$y, Vy.inv.y) / 2
  )

  if (isTRUE(trace)) {
    loglik <- get("loglik", envir = env)
    loglik <- c(loglik, res)
    assign("loglik", loglik, envir = env)
  }

  if (isTRUE(get.w)) {
    w <- psi * (X %*% Lambda) %*% crossprod(X, Vy.inv.y)
    assign("w", w, envir = env)
  }

  res
}
