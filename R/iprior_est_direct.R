################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2017  Haziq Jamil
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

iprior_direct <- function(mod, estimation.method, theta0, control) {
  # The direct optimisation method for estimating I-prior models. This minimises
  # the marginal deviance (-2 * logLik) using a quasi-Newton algorithm (L-BFGS).
  #
  # Args: An ipriorKernel object (mod), one of the direct estimation methods
  # e.g. loglik_direct or loglik_nystrom, the initial values (theta0) and a list
  # of control options to pass to optim.
  #
  # Returns: A list containing the optimised theta and parameters, loglik values,
  # standard errors, number of iterations, time taken, and convergence
  # information.
  iprior.env <- environment()
  w <- loglik <- NULL

  start.time <- Sys.time()
  res <- optim(theta0, estimation.method, object = mod, env = iprior.env,
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

loglik_iprior <- function(theta, object, trace = FALSE, env = NULL,
                          get.w = FALSE) {
  # The log-likelihood function for I-prior models.
  #
  # Args: theta (hyperparameters), object (an ipriorKernel object), and options
  # trace (logical) and env (the environment of optim) to be used with optim to
  # get the log-likelihood values and w (if get.w == TRUE)
  #
  # Returns: The log-likelihood value given theta of the I-prior model. Also, if
  # trace = TRUE then a vector of loglik values is written to the env
  # environment. If get.w = TRUE then a vector of w (posterior mean of the
  # I-prior random effects) are written to env as well.
  psi <- theta_to_psi(theta, object)
  eigen_Hlam(get_Hlam(object, theta), environment())

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

