get_Hlam <- function(object, theta = NULL) {
  which.poly <- is.kern_poly(object$kernels)
  if (!is.null(theta)) {
    tmp <- theta_to_param(theta, object)
    kernels <- tmp$kernels
    lambda <- tmp$lambda
  } else {
    kernels <- object$kernels
    lambda <- rep(1, object$p)
  }
  Hl <- get_Hl(object$Xl, kernels = kernels, lambda = lambda)

  # # if polynomial kernels are involved, then need to recalculate Hl
  # if (any(which.poly)) {
  #   Hl <- get_Hl(object$Xl, kernels = object$kernels, lambda = lambda)
  # } else {
  #   Hl <- object$Hl
  # }

  # Expand lambda here

  # Expand Hl here

  # Obtain Hlam
  res <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
  attr(res ,"kernel") <- "Hlam"
  res
}

eigen_Hlam <- function(Hlam) {
  tmp <- eigenCpp(Hlam)
  list(u = tmp$val, V = tmp$vec)
}

vy_inv_a <- function(u, V, a) {
  (V * rep(u, each = nrow(V))) %*% crossprod(V, a)
}

loglik_iprior <- function(theta, object, debug = FALSE) {
  psi <- theta_to_psi(theta)
  list2env(eigen_Hlam(get_Hlam(object, theta)), environment())

  y <- object$y  # y has already been standardised!
  n <- object$n
  z <- psi * u ^ 2 + 1 / psi
  logdet <- sum(log(z))
  Vy.inv.y <- vy_inv_a(1 / z, V, y)

  # for debug
 if (isTRUE(debug)) {
   tmp <- theta_to_param(theta, object)
   param <- as.matrix(tmp[, 1:4])
   print(c(as.numeric(na.omit(c(param))), psi))
 }

  res <- -n / 2 * log(2 * pi) - logdet / 2 - crossprod(y, Vy.inv.y) / 2
  as.numeric(res)
}


iprior2 <- function(y, ..., kernel = "linear", interactions = NULL,
                    est.lambda = TRUE, est.hurst = TRUE,
                    est.lengthscale = TRUE, est.offset = TRUE, control = list()) {
  mod <- kernL2(y = y, ..., kernel = kernel, est.lambda = est.lambda,
                est.hurst = est.hurst, est.lengthscale = est.lengthscale,
                est.offset = est.offset)

  control_ <- list(
    # fnscale   = -2,  # minimise the deviance
    # trace     = 1,    # trace of the optim
    theta0 = NULL
  )
  control.names <- names(control_)
  control_[(control.names <- names(control))] <- control
  control <- control_




  if (is.null(control$theta0)) {
    theta0 <- rnorm(mod$nt)
  } else {
    if (length(theta.start) != mod$nt) {
      stop(paste("Incorrect number of parameters specified. Should be", nt))
    }
  }

  iprior_direct(mod, loglik_iprior, theta0)
}

iprior_direct <- function(mod, estimation.method, theta.init) {
  res <- optim(theta.init, estimation.method, object = mod, debug = FALSE,
               method = "L-BFGS", control = list(fnscale = -2, trace = 1))
  param <- theta_to_param(res$par, mod)
  c(collapse_param(param)$param, psi = theta_to_psi(res$par))
}

# mod <- kernL2(stack.loss, stackloss$Air.Flow, stackloss$Water.Temp,
#               stackloss$Acid.Conc., kernel = "fbm")
