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

loglik_nystrom <- function(theta, object, trace = FALSE, env = NULL) {
  # The log-likelihood function for I-prior model using Nystrom approximation of
  # the kernel matrices.
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
  AB <- get_Hlam(object, theta)
  list2env(eigen_Hlam_nys(AB), environment())

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
