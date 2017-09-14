#' @export
logLik.ipriorKernel_Nystrom <- function(object, theta = NULL, Nys.adj = FALSE, ...) {
  # Initialise parameters ------------------------------------------------------
  alpha <- as.numeric(mean(object$Y))
  beta <- theta[-length(theta)]  # log(lambda)
  delta <- theta[length(theta)]  # log(psi)
  if (is.null(beta)) beta <- abs(rnorm(object$l, sd = 0.1))
  else {
    if (length(beta) != object$l) {
      stop(paste("Incorrect dimension of lambda initial values. vector of
                 length", object$l, "required."), call. = FALSE)
    }
  }
  if (is.null(delta)) delta <- abs(rnorm(1))

  # Calculate log-likelihood value ---------------------------------------------
  list2env(Nystrom_eigen(object, exp(beta), exp(delta)), envir = environment())
  if (isTRUE(Nys.adj)) z <- c(z, rep(1 / exp(delta), object$n - object$model$Nys.kern))
  logdet <- sum(log(z))
  res <- (-object$n / 2) * log(2 * pi) - logdet / 2 - crossprod(object$Y - alpha, a) / 2
  # print(c(exp(theta), res))
  as.numeric(res)
}

#' @export
logLik.ipriorMod_Nystrom <- function(object, theta = NULL, Nys.adj = FALSE, ...) {
  if (is.null(theta)) {
    if (isTRUE(Nys.adj)) object$loglik.adj
    else object$loglik
  }
  else logLik(object$ipriorKernel, theta, Nys.adj)
}

#' @export
deviance.ipriorKernel_Nystrom <- function(object, theta = NULL, Nys.adj = FALSE, ...) {
  -2 * logLik(object, theta)
}

#' @export
deviance.ipriorMod_Nystrom <- function(object, theta = NULL, Nys.adj = FALSE, ...) {
  if (is.null(theta)) {
    if (isTRUE(Nys.adj)) -2 * logLik(object)
    else -2 * logLik(object)
  }
  else deviance(object$ipriorKernel, theta, Nys.adj)
}

Nystrom_eigen <- function(object, lambda, psi) {
  this.env <- environment()
  list2env(object, this.env)
  list2env(model, this.env)
  environment(.lambdaExpand) <- environment(.lambdaContract) <- this.env
  Nys.m <- Nys.kern

  # Calculate Hlam.mat ---------------------------------------------------------
  if (q == 1) {
    hlamFn <- function(x = lambda, env = this.env) {
      assign("Hlam.mat", x[1] * Hl[[1]], envir = env)
    }
  }
  else {
    hlamFn <- function(x = lambda, env = this.env) {
      assign("Hlam.mat", Reduce("+", mapply("*", Hl[1:q], x[1:q],
                                            SIMPLIFY = FALSE)), envir = env)
    }
  }
  .lambdaExpand(x = lambda, env = this.env)
  hlamFn()

  # Eigendecomposition and Nystrom approximation -------------------------------
  A <- Hlam.mat[, 1:Nys.m]
  B <- Hlam.mat[, -(1:Nys.m)]
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
  z <- psi * u ^ 2 + 1 / psi
  a <- (V * rep(1 / z, each = nrow(V))) %*% (crossprod(V, Y - mean(Y)))

  list(u = u, V = V, z = z, a = a, A = A, B = B)
}
