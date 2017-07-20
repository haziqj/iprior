ipriorNystrom <- function(y, ..., size, seed = NULL, lambda.init = NULL,
                          psi.init = NULL) {
  UseMethod("ipriorNystrom")
}

ipriorNystrom.default <- function(y, ..., size, seed = NULL, lambda.init = NULL,
                                  psi.init = NULL, model = list()) {
  x.names <- as.character(as.list(match.call(expand.dots = FALSE))$...)
  ipriorNystrom.env <- environment()

  # Take subsamples and prepare kernel -----------------------------------------
  if (!is.null(seed)) set.seed(seed)
  Nys.samp <- sample(seq_along(y))
  model$Nys.kern <- size
  model$Nys.samp <- Nys.samp
  if (is.null(model$xname)) model$xname <- x.names
  model$kernel <- "FBM"
  mod <- kernL(y, ..., model = model)

  # Initialise parameters ------------------------------------------------------
  lambda <- lambda.init
  psi <- psi.init
  if (is.null(lambda)) lambda <- abs(rnorm(mod$l, sd = 0.1))
  else {
    if (length(lambda) != mod$l) {
      stop(paste("Incorrect dimension of lambda initial values. vector of
                 length", mod$l, "required."), call. = FALSE)
    }
  }
  if (is.null(psi)) psi <- abs(rnorm(1))
  theta <- c(lambda, psi)

  # Minimise deviance = -2 * log-likelihood ------------------------------------
  res.optim <- optim(theta, deviance, object = mod, method = "L-BFGS",
                     control = list(trace = 1), hessian = TRUE)

  # Results --------------------------------------------------------------------
  theta <- res.optim$par
  res <- list(alpha = mean(mod$Y), lambda = exp(theta[-length(theta)]),
              psi = exp(theta[length(theta)]),
              se = diag(solve(res.optim$hessian)), ipriorKernel = mod,
              loglik = -2 * res.optim$value, optim = res.optim)
  class(res) <- "ipriorMod_Nystrom"
  res
}

ipriorNystrom.formula <- function(formula, data, size, seed = NULL, lambda.init = NULL,
                                  psi.init = NULL) {

}

print.ipriorMod_Nystrom <- function(x, ...) {
  print(c(x$lambda, x$psi))
}

print.ipriorKernel_Nystrom <- function(x, ...) {
  class(x) <- "ipriorKernel"
  print(x)
}

logLik.ipriorKernel_Nystrom <- function(object, theta = NULL, ...) {
  # Initialise parameters ------------------------------------------------------
  alpha <- as.numeric(mean(object$Y))
  beta <- theta[-length(theta)]
  delta <- theta[length(theta)]
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
  # if(...) z <- c(z, rep(0, n - Nys.m))
  logdet <- sum(log(z))
  res <- (-n / 2) * log(2 * pi) - logdet / 2 - crossprod(object$Y - alpha, a) / 2
  print(c(exp(theta), res))
  as.numeric(res)
}

deviance.ipriorKernel_Nystrom <- function(object, theta = NULL, ...) {
  return(-2 * logLik(object, theta))
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
  if (any(tmp1$val + 1e-7 < 0)) print(tmp1$val + 1e-7)
  U <- tmp1$vectors
  C.tmp <- U * rep(1 / sqrt(tmp1$val + 1e-7), each = nrow(U))
  C <- C.tmp %*% crossprod(U, B)
  Q <- A + tcrossprod(C)
  tmp2 <- eigenCpp(Q)
  if (any(tmp2$val + 1e-7 < 0)) print(tmp2$val + 1e-7)
  u <- tmp2$values + 1e-7
  R <- tmp2$vectors
  V <- rbind(A, t(B)) %*% tcrossprod(C.tmp, U) %*%
    (R * rep(1 / sqrt(u), each = nrow(R)))
  z <- psi * u ^ 2 + 1 / psi
  a <- (V * rep(1 / z, each = nrow(V))) %*% (crossprod(V, Y - mean(Y)))

  list(u = u, V = V, z = z, a = a)
}

fitted.ipriorMod_Nystrom <- function(object, ...) {
  this.env <- environment()
  list2env(object, this.env)
  list2env(Nystrom_eigen(ipriorKernel, lambda, psi), envir = this.env)
  y.hat <- alpha + psi * (V * rep(u ^ 2, each = nrow(V))) %*% (t(V) %*% a)
  # Effectively this is SR estimates...
  as.numeric(y.hat)[order(object$ipriorKernel$model$Nys.samp)]
}

plot.ipriorMod_Nystrom <- function() {
  1
}
