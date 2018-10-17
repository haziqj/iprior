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

iprior_mixed <- function(mod, theta0 = NULL, em.maxit = 5, stop.crit = 1e-5,
                         silent = FALSE, control.optim = list(),
                         optim.method = "L-BFGS") {
  # This method does several EM steps before switching to the direct
  # optimisation method. It combines both the stability of the EM with the fast
  # convergence of the direct method.
  #
  # Args: An ipriorKernel object (mod), number of maximum EM iterations
  # (em.maxit), the stopping criterion to determine convergence (stop.crit), an
  # option to suppress reporting (silent), the initial values (theta0) and an
  # optional list of controls for optim.
  #
  # Returns: A list containing the optimised theta and parameters, loglik
  # values, standard errors, number of iterations, time taken, and convergence
  # information.

  # Default optim control list -------------------------------------------------
  control.optim_ <- list(
    fnscale = -2,
    trace   = ifelse(isTRUE(silent), 0, 1),
    maxit   = 100,
    REPORT  = 10,
    factr   = 1e7
  )
  control.optim <- update_control(control.optim, control.optim_)

  # First pass to EM routine ---------------------------------------------------
  if (!isTRUE(silent))
    cat(paste0("Running ", em.maxit, " initial EM iterations\n"))
  start.time <- Sys.time()
  em.method <- iprior_method_checker(mod, "em")
  if (em.method["em.closed"]) {
    tmp <- iprior_em_closed(mod, em.maxit, stop.crit, silent, theta0,
                            mixed = TRUE)
  }
  if (em.method["em.reg"]) {
    tmp <- iprior_em_reg(mod, em.maxit, stop.crit, silent, theta0, mixed = TRUE)
  }

  # Then pass to direct maximisation routine -----------------------------------
  if (!isTRUE(silent)) cat("Now switching to direct optimisation\n")
  res <- iprior_direct(mod, loglik_iprior, tmp$theta, control.optim,
                       optim.method)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Update time, call, maxit, niter, lb, error, brier --------------------------
  res$time <- time.taken
  res$start.time <- start.time
  res$end.time <- end.time
  res$loglik <- c(tmp$loglik, res$loglik)

  res
}
