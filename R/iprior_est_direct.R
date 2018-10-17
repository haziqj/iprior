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

iprior_direct <- function(mod, lik.fn, theta0, control = list(),
                          method = "L-BFGS") {
  # The direct optimisation method for estimating I-prior models. This minimises
  # the marginal deviance (-2 * logLik) using a quasi-Newton algorithm (L-BFGS).
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
    fnscale = -2,  # minimising log-likelihood on deviance scale
    trace   = 1,
    maxit   = 100,
    REPORT  = 10,
    factr   = 1e7
  )
  control <- update_control(control, control_)

  iprior.env <- environment()
  w <- loglik <- NULL

  start.time <- Sys.time()
  res <- optim(theta0, lik.fn, object = mod, env = iprior.env,
               trace = TRUE, method = method, control = control,
               hessian = TRUE)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Calculate fitted values  ---------------------------------------------------
  w <- get_w(u, V, Vy.inv.y, psi)  # NOTE: psi, V, u and Vy.inv.y have been
  y.hat <- get_y.hat(u, V, w)      # assigned to environment by loglik

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

loglik_iprior <- function(theta, object, trace = FALSE, env = NULL) {
  # The log-likelihood function for I-prior models.
  #
  # Args: theta (hyperparameters), object (an ipriorKernel object), and optional
  # trace (logical) and env (the environment of optim) to be used with optim in
  # order to append the log-likelihood value to the existing vector, and also
  # obtain the stuff required to calculate w and y.hat in the env environment.
  #
  # Returns: The log-likelihood value given theta of the I-prior model. Also, if
  # trace = TRUE then a vector of loglik values is written to the env
  # environment, along with Vy.invy, u, V and psi.
  psi <- theta_to_psi(theta, object)
  eigen_Hlam(get_Hlam(object, theta), environment())

  y <- object$y  # y has already been standardised!
  n <- object$n

  z <- psi * u ^ 2 + 1 / psi
  logdet <- sum(log(z))
  Vy.inv.y <- A_times_a(1 / z, V, y)

  res <- -n / 2 * log(2 * pi) - logdet / 2 - crossprod(y, Vy.inv.y) / 2

  if (isTRUE(trace)) {
    loglik <- get("loglik", envir = env)
    loglik <- c(loglik, res)
    assign("loglik", loglik, envir = env)
    assign("Vy.inv.y", Vy.inv.y, envir = env)
    assign("u", u, envir = env)
    assign("V", V, envir = env)
    assign("psi", psi, envir = env)
  }

  as.numeric(res)
}

