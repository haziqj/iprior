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

#' @export
iprior_cv <- function(...) UseMethod("iprior_cv")

#' Perform a cross-validation experiment with the iprior function
#'
#' A convenience function to perform a k-fold cross-validation experiment and
#' obtain mean squared error of prediction. Most of the arguments are similar to
#' \code{iprior()} and \code{kernL()}.
#'
#' Uses a multicore loop to fit the folds by default, set \code{par.cv = FALSE}
#' to not use multithreading.
#'
#' @inheritParams iprior
#' @param folds The number of cross-validation folds. Set equal to sample size
#'   or \code{Inf} to perform leave-one-out cross-validation.
#' @param par.cv Logical. Multithreading to fit the models? Defaults to
#'   \code{TRUE}.
#'
#' @return An \code{iprior_xv} object containing a data frame of the
#'   cross-validated values such as the log-likelihood, training MSE and test
#'   MSE.
#'
#' @examples
#' \dontrun{
#'
#' # 5-fold CV experiment
#' (mod.cv <- iprior_cv(y ~ X, gen_smooth(100), kernel = "se", folds = 5))
#'
#' # LOOCV experiment
#' (mod.cv <- iprior_cv(y ~ X, gen_smooth(100), kernel = "se", folds = Inf))
#'
#' # Can also get root MSE
#' print(mod.cv, "RMSE")
#' }
#'
#' @name iprior_cv
#' @export
iprior_cv.default <- function(y, ..., folds = 2, par.cv = TRUE,
                              kernel = "linear", method = "direct",
                              control = list(), interactions = NULL,
                              est.lambda = TRUE, est.hurst = FALSE,
                              est.lengthscale = FALSE, est.offset = FALSE,
                              est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                              psi = 1, nystrom = FALSE, nys.seed = NULL) {
  n <- length(y)
  if (folds > n) folds <- n
  if (folds <= 1) {
    warning("Number of folds needs to be 2 or greater. Defaulting to 2.")
    folds <- 2
  }

  # Get the test samples for cross-validation ----------------------------------
  samp <- seq_len(n)
  n.cv <- rep(n %/% folds, folds)
  n.cv[1] <- n.cv[1] + n %% folds
  n.cv.cum <- c(0, cumsum(n.cv))
  cv.samp <- list(NULL)
  for (i in seq_len(folds)) {
    cv.samp[[i]] <- samp[(n.cv.cum[i] + 1):n.cv.cum[i + 1]]
  }

  # Checks and settings --------------------------------------------------------
  original.silent <- control$silent
  control$silent <- TRUE
  if (!isTRUE(original.silent)) {
    pb <- txtProgressBar(min = 0, max = folds, style = 1)
    snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
  } else {
    snow.options.list <- list()
  }

  if (isTRUE(control$restarts)) par.cv <- FALSE

  # Fit full model to get proper starting values -------------------------------
  if (is.null(control$theta0)) {
    control.full <- control
    control.full$silent <- original.silent
    control.full$restarts <- TRUE
    if (!isTRUE(original.silent)) cat("Fitting full model\n")
    full.mod <- iprior(y = y, ..., kernel = kernel, interactions = interactions,
                       est.lambda = est.lambda, est.hurst = est.hurst,
                       est.lengthscale = est.lengthscale, est.offset = est.offset,
                       est.psi = est.psi, fixed.hyp = fixed.hyp, lambda = lambda,
                       psi = psi, nystrom = nystrom, nys.seed = nys.seed,
                       method = method, control = control.full)
    control$theta0 <- get_theta(full.mod)
  }

  # The cross-validation routine -----------------------------------------------
  if (!isTRUE(original.silent)) cat("Performing cross-validation routine\n")
  if (!isTRUE(par.cv)) {
    # The non-multithreading version
    res <- as.data.frame(matrix(NA, ncol = 3, nrow = folds))
    colnames(res) <- c("Log-lik", "Training", "Test")
    for (k in seq_along(cv.samp)) {
      the.train.samp <- seq_along(y)[-cv.samp[[k]]]
      mod <- iprior.ipriorKernel(
        kernL(y = y, ..., train.samp = the.train.samp, kernel = kernel,
              interactions = interactions, est.lambda = est.lambda,
              est.hurst = est.hurst, est.lengthscale = est.lengthscale,
              est.offset = est.offset, est.psi = est.psi,
              fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
              nystrom = nystrom, nys.seed = nys.seed),
        method = method, control = control
      )
      res[k, ] <- c(logLik(mod), get_rmse(mod))
      if (!isTRUE(original.silent)) setTxtProgressBar(pb, k)
    }
  } else {
    # The Multithreading version
    cl <- parallel::makeCluster(parallel::detectCores())
    doSNOW::registerDoSNOW(cl)
    res <- foreach::`%dopar%`(
      foreach::foreach(
        k = seq_along(cv.samp),
        .combine = rbind,
        .packages = "iprior",
        .options.snow = snow.options.list,
        .errorhandling = "remove"
      ), {
        the.train.samp <- seq_along(y)[-cv.samp[[k]]]
        mod <- iprior.ipriorKernel(
          kernL(y = y, ..., train.samp = the.train.samp, kernel = kernel,
                interactions = interactions, est.lambda = est.lambda,
                est.hurst = est.hurst, est.lengthscale = est.lengthscale,
                est.offset = est.offset, est.psi = est.psi,
                fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
                nystrom = nystrom, nys.seed = nys.seed),
          method = method, control = control
        )
        c("Log-lik" = logLik(mod), get_rmse(mod))
      }
    )
    parallel::stopCluster(cl)
  }
  if (!isTRUE(original.silent)) close(pb)

  structure(list(res = res, folds = folds, n = n), class = "iprior_xv")
}

#' @rdname iprior_cv
#' @export
iprior_cv.formula <- function(formula, data, folds = 2, one.lam = FALSE,
                              par.cv = TRUE, kernel = "linear",
                              method = "direct", control = list(),
                              est.lambda = TRUE, est.hurst = FALSE,
                              est.lengthscale = FALSE, est.offset = FALSE,
                              est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                              psi = 1, nystrom = FALSE, nys.seed = NULL, ...) {
  list2env(formula_to_xy(formula = formula, data = data, one.lam = one.lam),
           envir = environment())
  iprior_cv.default(y, Xl.formula = Xl, interactions = interactions,
                    folds = folds, par.cv = par.cv, kernel = kernel,
                    method = method, control = control, est.lambda = est.lambda,
                    est.hurst = est.hurst, est.lengthscale = est.lengthscale,
                    est.offset = est.offset, est.psi = est.psi,
                    fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
                    nystrom = nystrom, nys.seed = nys.seed, ...)
}

#' @export
print.iprior_xv <- function(x, result = c("RMSE", "MSE"), ...) {
  cv.method <- ifelse(x$folds == x$n, "Leave-one-out Cross Validation",
                      paste0(x$folds, "-fold Cross Validation"))
  result <- match.arg(toupper(result), c("RMSE", "MSE"))
  res <- x$res
  if (result == "MSE") res[, -1] <- res[, -1] ^ 2
  res <- apply(res, 2, mean)

  cat("Results from", cv.method, "\n")
  cat(paste0("Training ", result, ":"), res[2], "\n")
  cat(paste0("Test ", result, "    :"), res[3], "\n")
}




