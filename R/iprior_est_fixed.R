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

iprior_fixed <- function(mod) {
  # When all that is required is to obtain the posterior estimate of the
  # regression function (i.e. no estimation of any hyperparameters), then this
  # estimatation method does just that.
  #
  # Args: An ipriorKernel object (mod)
  #
  # Returns: A list containing the optimised theta and parameters, loglik
  # values, standard errors, number of iterations, time taken, and convergence
  # information.
  w <- loglik <- NULL

  start.time <- Sys.time()
  if (is.ipriorKernel_nys(mod)) {
    loglik_nystrom(mod$thetal$theta, mod, trace = TRUE, get.w = TRUE,
                   env = environment())
  } else {
    loglik_iprior(mod$thetal$theta, mod, trace = TRUE, get.w = TRUE,
                  env = environment())
  }
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  param.full <- theta_to_collapsed_param(mod$thetal$theta, mod)

  list(theta = NULL, param.full = param.full, loglik = loglik,
       se = NA, niter = NA, w = as.numeric(w), start.time = start.time,
       end.time = end.time, time = time.taken, convergence = NA, message = NA,
       niter = NA)
}
