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

iprior_parallel <- function(mod, method = "direct",
                            control = list(silent = FALSE, restarts = TRUE,
                                           no.cores = parallel::detectCores(),
                                           par.maxit = 100,
                                           optim.method = "L-BFGS")) {
  # This estimation method is not really any specific method per se, but rather
  # a  function that enables the selected method to be run multiple times from
  # different random starting points. The starting point that gives the highest
  # likelihood value after three iterations is selected.
  #
  # Args: An ipriorKernel object (mod), the I-prior esitmation method, and a
  # list of controls (the some one in the main iprior() function).
  #
  # Returns: A list containing the optimised theta and parameters, loglik
  # values, standard errors, number of iterations, time taken, and convergence
  # information.

  # Set up controls ------------------------------------------------------------
  if (control$restarts == 1) control$restarts <- parallel::detectCores()
  control$no.cores <- min(parallel::detectCores(), control$restarts)
  if (!is.null(control$theta0)) {
    message("Ignoring theta0 control options with random restarts.")
  }
  control$theta0 <- NULL
  if (!isTRUE(control$silent)) {
    cat("Performing", control$restarts, "random restarts on", control$no.cores,
        "cores\n")
    snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
    pb <- txtProgressBar(min = 0, max = control$restarts, style = 1)
  } else {
    snow.options.list <- list()
  }

  # The multithreading bit -----------------------------------------------------
  start.time <- Sys.time()
  cl <- parallel::makeCluster(control$no.cores)
  doSNOW::registerDoSNOW(cl)
  res <- foreach::`%dopar%`(
    foreach::foreach(
      i = seq_len(control$restarts),
      .packages = "iprior",
      .options.snow = snow.options.list,
      .errorhandling = "remove"
    ), {
      new.control          <- control
      new.control$restarts <- 0
      new.control$maxit    <- control$par.maxit
      new.control$silent   <- TRUE
      tmp <- iprior(mod, control = new.control, method = method)
      list(
        theta  = tmp$theta,
        loglik = tmp$loglik,
        niter  = tmp$niter
      )
    }
  )
  if (!isTRUE(control$silent)) close(pb)
  parallel::stopCluster(cl)

  # Find best starting value ---------------------------------------------------
  tmp <- sapply(res, function(x) x$loglik[length(x$loglik)])
  names(tmp) <- paste("Run", seq_along(tmp))
  best.run <- which(tmp == max(tmp))
  best.niter <- res[[best.run]]$niter
  best.loglik <- res[[best.run]]$loglik
  if (!isTRUE(control$silent)) {
    cat("Log-likelihood from random starts:\n")
    print(tmp)
    cat("Continuing on Run", best.run, "\n")
  }

  # Continue updating the best model -------------------------------------------
  control$restarts <- 0
  control$theta0   <- res[[best.run]]$theta
  control$maxit    <- control$maxit - 3
  res <- iprior(mod, method = method, control = control)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)

  # Update time taken ----------------------------------------------------------
  res$time <- time.taken
  res$start.time <- start.time
  res$end.time <- end.time
  res$niter <- res$niter + best.niter
  res$loglik <- c(best.loglik, res$loglik)

  res
}
