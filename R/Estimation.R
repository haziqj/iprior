get_Hlam <- function(object, theta = NULL) {
  which.poly <- is.kern_poly(object$kernels)
  if (!is.null(theta)) {
    tmp <- theta_to_param(theta, object$theta.start$na, object$which.pearson,
                          object$poly.degree)
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
  V %*% diag(u, nrow(V)) %*% t(V) %*% a
}

loglik_iprior <- function(theta = NULL, object, debug = FALSE) {
  logpsi <- 0
  if (!is.null(theta)) {
    logpsi <- theta[length(theta)]
    theta <- theta[-length(theta)]
  }
  psi <- exp(logpsi)
  list2env(eigen_Hlam(get_Hlam(object, theta)), environment())
  y <- object$y
  n <- object$n
  alpha <- attr(y, "scaled:center")
  y.cen <- as.numeric(y) - alpha
  z <- psi * u ^ 2 + 1 / psi
  logdet <- sum(log(z))
  Vy.inv.y <- vy_inv_a(1 / z, V, y.cen)

  # for debug
 if (isTRUE(debug)) {
   tmp <- theta_to_param(theta, object$theta.start$na, object$which.pearson,
                         object$poly.degree)
   param <- as.matrix(tmp[, 1:4])
   print(c(as.numeric(na.omit(c(param))), psi))
 }

  res <- -n / 2 * log(2 * pi) - logdet / 2 - crossprod(y.cen, Vy.inv.y) / 2
  as.numeric(res)
}

iprior2 <- function(object) {
  res <- optim(rnorm(length(object$theta.start$theta) + 1),
               loglik_iprior, object = object, debug = FALSE,
               method = "L-BFGS", control = list(fnscale = -1))
  res
}
