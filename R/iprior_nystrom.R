loglik_nystrom <- function(theta, object, trace = FALSE, env = NULL,
                           get.w = FALSE) {
  # Args: theta (hyperparameters), object (an ipriorKernel object), and options
  # trace (logical) and env (the environment of optim) to be used with optim to
  # get the log-likelihood values and w (if get.w == TRUE)
  #
  # Output: The log-likelihood value given theta of the I-prior model.
  psi <- theta_to_psi(theta, object)
  AB <- get_Hlam(object, theta, get_Xl.nys(object))
  list2env(eigen_Hlam_nys(AB), environment())

  y <- object$y  # y has already been standardised!
  n <- object$n
  z <- psi * u ^ 2 + 1 / psi
  logdet <- sum(log(z))
  Vy.inv.y <- vy_inv_a(1 / z, V, y)

  res <- -n / 2 * log(2 * pi) - logdet / 2 - crossprod(y, Vy.inv.y) / 2

  if (isTRUE(trace)) {
    loglik <- get("loglik", envir = env)
    loglik <- c(loglik, res)
    assign("loglik", loglik, envir = env)
  }

  if (isTRUE(get.w)) {
    w <- psi * (V %*% ((t(V) * u) %*% Vy.inv.y))
    assign("w", w, envir = env)
  }

  as.numeric(res)
}

eigen_Hlam_nys <- function(Hlam) {
  # Args: The (truncated) scaled kernel matrix Hlam, which is m x n in size,
  # where m = object$nys.size
  #
  # Output: A list of the eigen values and vectors of Hlam.
  m <- nrow(Hlam)
  A <- Hlam[, seq_len(m)]
  B <- Hlam[, -seq_len(m)]
  tmp1 <- eigenCpp(A)
  # print(c(lambda, psi))
  # if (any(tmp1$val + 1e-7 < 0)) print(tmp1$val + 1e-7)
  U <- tmp1$vectors
  C.tmp <- U * rep(1 / sqrt(tmp1$val + 1e-7), each = nrow(U))
  C <- C.tmp %*% crossprod(U, B)
  Q <- A + tcrossprod(C)
  tmp2 <- eigenCpp(Q)
  # if (any(tmp2$val + 1e-7 < 0)) print(tmp2$val + 1e-7)
  u <- tmp2$values + 1e-7
  R <- tmp2$vectors
  V <- rbind(A, t(B)) %*% tcrossprod(C.tmp, U) %*%
    (R * rep(1 / sqrt(u), each = nrow(R)))

  list(u = u, V = V)
}
