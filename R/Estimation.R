get_Hlam <- function(object, theta, xstar = list(NULL)) {
  # which.poly <- is.kern_poly(object$kernels)
  tmp <- theta_to_param(theta, object)
  kernels <- tmp$kernels
  lambda <- tmp$lambda

  # do we need to calculate Hl entries each time? maybe so for when estimating
  # hurst, offset, lengthscale, and polynomial kernels.
  Hl <- get_Hl(object$Xl, xstar, kernels = kernels, lambda = lambda)
  expand_Hl_and_lambda(Hl, lambda, object$intr, object$intr.3plus, environment())

  res <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
  kernels.to.add <- sapply(Hl, attr, "kernel")
  attr(res ,"kernel") <- paste(kernels.to.add, collapse = " + ")
  res
}

eigen_Hlam <- function(Hlam) {
  # Args: The scaled kernel matrix Hlam.
  #
  # Output: A list of the eigen values and vectors of Hlam.
  tmp <- eigenCpp(Hlam)
  list(u = tmp$val, V = tmp$vec)
}

vy_inv_a <- function(u, V, a) {
  # Args: Eigenvalues u, eigenvectors V (of Vy), and a vector a.
  #
  # Output: A vector Vy^{-1} %*% a.
  (V * rep(u, each = nrow(V))) %*% crossprod(V, a)
}

loglik_iprior <- function(theta, object, trace = FALSE, env = NULL,
                          get.w = FALSE) {
  # Args: theta (hyperparameters), object (an ipriorKernel object), and options
  # trace (logical) and env (the environment of optim) to be used with optim to
  # get the log-likelihood values and w (if get.w == TRUE)
  #
  # Output: The log-likelihood value given theta of the I-prior model.
  psi <- theta_to_psi(theta, object)
  list2env(eigen_Hlam(get_Hlam(object, theta)), environment())

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




iprior_direct <- function(mod, estimation.method, theta.init, control) {
  # Args: An ipriorKernel object (mod), one of the direct estimation methods
  # e.g. loglik_direct or loglik_nystrom, the initial values (theta.init) and a
  # list of control options to pass to optim.
  #
  # Output: A list containing theta, loglik, se, niter, and time.
  this.env <- environment()
  w <- loglik <- NULL
  start.time <- Sys.time()
  res <- optim(theta.init, estimation.method, object = mod, env = this.env,
               trace = TRUE, get.w = TRUE, method = "L-BFGS", control = control,
               hessian = TRUE)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  tmp <- eigenCpp(-res$hessian)
  u <- tmp$val + 1e-9
  V <- tmp$vec
  Fi <- V %*% t(V) / u
  se <- sqrt(diag(Fi))
  se <- convert_se(se, res$par, mod)  # delta method to convert to parameter s.e.
  loglik <- as.numeric(na.omit(loglik))
  param.full <- theta_to_collapsed_param(res$par, mod)
  list(theta = res$par, param.full = param.full, loglik = loglik,
       se = se, niter = res$count[1], w = as.numeric(w), start.time = start.time,
       end.time = end.time, time = time.taken, convergence = res$convergence,
       message = res$message)
}

convert_se <- function(se, theta, object) {
  theta.names <- names(object$theta)
  res <- se

  types <- c("lambda", "hurst", "offset", "lengthscale", "psi")
  for (i in seq_along(types)) {
    type <- types[i]
    ind <- grep(type, theta.names)
    if (type == "lambda" & length(ind) == 1)
      res[ind] <- se[ind] * exp(theta[ind])
    else if (type == "hurst")
      res[ind] <- se[ind] * dnorm(theta[ind])
    else
      res[ind] <- se[ind] * exp(theta[ind])
  }

  res
}

iprior_fixed <- function(mod) {
  w <- loglik <- NULL
  start.time <- Sys.time()
  loglik_iprior(mod$thetal$theta, mod, trace = TRUE, get.w = TRUE,
                env = environment())
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  param.full <- theta_to_collapsed_param(mod$thetal$theta, mod)
  list(theta = NULL, param.full = param.full, loglik = loglik,
       se = NA, niter = NA, w = as.numeric(w), start.time = start.time,
       end.time = end.time, time = time.taken, convergence = NA, message = NA,
       niter = NA)
}
