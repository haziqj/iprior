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

iprior_canonical <- function(mod, theta0 = NULL, control = list()) {
  # The efficient optimisation method for estimating I-prior models with
  # canonical kernels ONLY (no interactions). Like iprior_direct(), this
  # minimises the marginal deviance (-2 * logLik) using a quasi-Newton algorithm
  # (L-BFGS).
  #
  # Args: An ipriorKernel object (mod), one of the direct estimation methods
  # e.g. loglik_direct or loglik_nystrom, the initial values (theta0) and a list
  # of control options to pass to optim.
  #
  # Returns: A list containing the optimised theta and parameters, y.hat and w,
  # loglik values, standard errors, number of iterations, time taken, and
  # convergence information.

  # Default optim control list -------------------------------------------------
  control_ <- list(
    fnscale = -2,
    trace   = TRUE,
    maxit   = 100,
    REPORT  = 10,
    factr   = 1e7
  )
  control <- update_control(control, control_)

  iprior.env <- environment()
  w <- loglik <- NULL
  start.time <- Sys.time()
  res <- optim(theta0, loglik_canonical, object = mod, env = iprior.env,
               trace = TRUE, method = "L-BFGS", control = control,
               hessian = TRUE)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate fitted values  ---------------------------------------------------
  w <- get_w(lambda, X, Vy.inv.y, psi)  # NOTE: psi, X, lambda and Vy.inv.y have
  y.hat <- get_y.hat(lambda, X, w)      # been assigned to environment by loglik

  # Calculate standard errors --------------------------------------------------
  tmp <- eigenCpp(-res$hessian)
  u <- tmp$val + 1e-9
  V <- tmp$vec
  Fi.inv <- V %*% (t(V) / u)
  se <- sqrt(diag(Fi.inv))
  se <- convert_se(se, res$par, mod)  # delta method to convert to parameter s.e.

  # Clean up and close ---------------------------------------------------------
  loglik <- as.numeric(na.omit(loglik))
  param.full <- theta_to_collapsed_param(res$par, mod)

  list(theta = res$par, param.full = param.full, se = se, loglik = loglik,
       w = w, y.hat = y.hat, niter = res$count[1], start.time = start.time,
       end.time = end.time, time = time.taken, convergence = res$convergence,
       message = res$message)
}

loglik_canonical <- function(theta, object, trace = FALSE, env = NULL) {
  # The log-likelihood function for I-prior models with linear kernels only (no
  # interanctions either).
  #
  # Args: theta (hyperparameters), object (an ipriorKernel object), and optional
  # trace (logical) and env (the environment of optim) to be used with optim in
  # order to append the log-likelihood value to the existing vector, and also
  # obtain the stuff required to calculate w and y.hat in the env environment.
  #
  # Returns: The log-likelihood value given theta of the I-prior model. Also, if
  # trace = TRUE then a vector of loglik values is written to the env
  # environment, along with Vy.invy, lambda, X and psi.
  psi <- theta_to_psi(theta, object)
  lambda <- theta_to_collapsed_param(theta, object)[seq_len(object$p)]
  X <- matrix(unlist(object$Xl), nrow = object$n)
  X <- scale(X, center = TRUE, scale = FALSE)
  XtX <- crossprod(X)
  Lambda <- diag(lambda, ncol(X))
  Vy.inv <- psi * (
    diag(1, object$n) - X %*% solve(solve(psi ^ 2 * Lambda %*% XtX %*% Lambda) + XtX, t(X))
  )
  Vy.inv.y <- Vy.inv %*% object$y
  logdet <- determinant(Vy.inv)$mod
  res <- as.numeric(
    (-object$n / 2) * log(2 * pi) + logdet / 2 - crossprod(object$y, Vy.inv.y) / 2
  )

  if (isTRUE(trace)) {
    loglik <- get("loglik", envir = env)
    loglik <- c(loglik, res)
    assign("loglik", loglik, envir = env)
    assign("Vy.inv.y", Vy.inv.y, envir = env)
    assign("X", X, envir = env)
    assign("lambda", lambda, envir = env)
    assign("psi", psi, envir = env)
  }

  res
}
