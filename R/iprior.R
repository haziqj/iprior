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

#' @export
iprior <- function(...) UseMethod("iprior")
# iprior <- function(y, ..., kernel = "linear", method = c("direct", "em",
#                                                          "fixed", "canonical",
#                                                          "mixed"),
#                    control = list(), interactions = NULL, est.lambda = TRUE,
#                    est.hurst = FALSE, est.lengthscale = FALSE,
#                    est.offset = FALSE, est.psi = TRUE, fixed.hyp = NULL,
#                    lambda = 1, psi = 1, nystrom = FALSE, nys.seed = NULL,
#                    object, iter.update, model = list()) {
#   UseMethod("iprior")
# }

#' Fit an I-prior regression model
#'
#' A function to perform regression using I-priors. The I-prior model parameters
#' may be estimated in a number of ways: direct minimisation of the marginal
#' deviance, EM algorithm, fixed hyperparameters, or using a Nystrom kernel
#' approximation.
#'
#' The \code{iprior()} function is able to take formula based input and
#' non-formula. When not using formula, the syntax is as per the default S3
#' method. That is, the response variable is the vector \code{y}, and any
#' explanatory variables should follow this, and separated by commas.
#'
#' As described \link[=kernL2]{here}, the model can be loaded first into an
#' \code{ipriorKernel} object, and then passed to the \code{iprior()} function
#' to perform the estimation.
#'
#' @inheritParams kernL2
#' @param object An \code{ipriorKernel2} or \code{ipriorMod} object.
#' @param method The estimation method. One of: \itemize{ \item{\code{"direct"}
#'   - for the direct minimisation of the marginal deviance using
#'   \code{optim()}'s L-BFGS method} \item{\code{"em"} - for the EM algorithm}
#'   \item{\code{"mixed"} - combination of the direct and EM methods}
#'   \item{\code{"fixed"} - for just obtaining the posterior regression function
#'   with fixed hyperparameters (default method when setting \code{fixed.hyp =
#'   TRUE})} \item{\code{"canonical"} - an efficient estimation method which
#'   takes advantage of the structure of the linear kernel} }
#' @param control 1
#' @param one.lam 1
#' @param iter.update 1
#'
#' @return 1
#'
#' @examples 1
#'
#' @name iprior
#' @export
iprior.default <- function(y, ..., kernel = "linear", method = c("direct", "em",
                                                                 "fixed",
                                                                 "canonical",
                                                                 "mixed"),
                           control = list(), interactions = NULL,
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                           psi = 1, nystrom = FALSE, nys.seed = NULL,
                           model = list()) {
  # Load the I-prior model -----------------------------------------------------
  if (is.ipriorKernel2(y)) {
    mod <- y
  } else {
    mod <- kernL2(y = y, ..., kernel = kernel, interactions = interactions,
                  est.lambda = est.lambda, est.hurst = est.hurst,
                  est.lengthscale = est.lengthscale, est.offset = est.offset,
                  est.psi = est.psi, fixed.hyp = fixed.hyp, lambda = lambda,
                  psi = psi, nystrom = nystrom, nys.seed = nys.seed,
                  model = model)
  }

  # Set up controls ------------------------------------------------------------
  method <- tolower(method)
  method <- match.arg(method, c("direct", "em", "fixed", "canonical", "mixed"))

  control_ <- list(
    maxit     = 100,
    em.maxit  = 5,  # for mixed method
    stop.crit = 1e-8,  # sqrt(.Machine$double.eps), roughly 1e-8
    theta0    = NULL,
    # lambda0   = NULL,
    # psi0      = NULL,
    silent    = FALSE,
    report    = 10,
    psi.reg   = FALSE,  # option for iprior_em_reg()
    restarts  = 0,
    no.cores  = parallel::detectCores()
  )
  control <- update_control(control, control_)  # see iprior_helper.R
  control.optim <- list(
    fnscale = -2,
    trace   = ifelse(isTRUE(control$silent), 0, 1),
    maxit   = max(0, control$maxit - 1),
    REPORT  = control$report
  )

  # Starting values ------------------------------------------------------------
  if (is.null(control$theta0)) {
    theta0 <- rnorm(mod$thetal$n.theta)  # rep(0, mod$thetal$n.theta)
  } else {
    theta0 <- control$theta0
    if (length(theta0) != mod$thetal$n.theta) {
      stop(paste("Incorrect number of parameters specified. Should be",
                 mod$thetal$n.theta))
    }
  }

  # Send to correct estimation method ------------------------------------------
  if (as.numeric(control$restarts) >= 1) {
    res <- iprior_parallel(mod, method, control)
    res$est.method <- paste0(
      gsub("\\.", "", res$est.method), " with random restarts."
    )
  } else {
    est.method <- iprior_method_checker(mod, method)
    if (est.method["fixed"]) {
      res <- iprior_fixed(mod)
      res$est.method <- "I-prior fixed."
      res$est.conv <- ""
    } else if (est.method["canonical"]) {
      res <- iprior_canonical(mod, theta0, control.optim)
      res$est.method <- "Efficient canonical method."
    } else {
      if (est.method["em.closed"]) {
        res <- iprior_em_closed(mod, control$maxit, control$stop.crit,
                                control$silent, theta0)
        res$est.method <- "Closed-form EM algorithm."
      }
      if (est.method["em.reg"]) {
        res <- iprior_em_reg(mod, control$maxit, control$stop.crit,
                             control$silent, theta0)
        res$est.method <- "Regular EM algorithm."
      }
      if (est.method["direct"]) {
        res <- iprior_direct(mod, loglik_iprior, theta0, control.optim)
        res$est.method <- "Direct optimisation method."
      }
      if (est.method["nystrom"]) {
        res <- iprior_direct(mod, loglik_nystrom, theta0, control.optim)
        res$est.method <- "Nystrom approximated optimisation."
      }
      if (est.method["mixed"]) {
        res <- iprior_mixed(mod, theta0, control$em.maxit, control$stop.crit,
                            control$silent, control.optim)
        res$est.method <- paste0("EM algorithm (", control$em.maxit,
                                 " steps) + direct minimisation.")
      }
      if (res$conv == 0)
        res$est.conv <- paste0("Converged to within ", control$stop.crit,
                               " tolerance.")
      else
        res$est.conv <- "Convergence criterion not met."
    }
  }

  intercept <- attr(mod$y, "scaled:center")
  if (is.null(intercept)) intercept <- mean(mod$y)
  res$intercept <- intercept
  res$coefficients <- reduce_theta(res$param.full, mod$estl)$theta.reduced
  tmp <- predict_iprior(mod$y, get_Hlam(mod, res$theta), res$w, res$intercept)
  res$fitted.values <- tmp$y
  names(res$fitted.values) <- attr(mod$y, "dimnames")[[1]]
  res$residuals <- tmp$resid
  res$train.error <- tmp$train.error
  res$ipriorKernel <- mod
  res$method <- method
  res$control <- control

  res$fullcall <- match.call()
  cl <- mod$call
  cl[[1L]] <- as.name("iprior")
  # names(cl)[2] <- ""  # get rid of "y ="
  res$call <- cl

  class(res) <- "ipriorMod"
  res
}

#' @rdname iprior
#' @export
iprior.formula <- function(formula, data, kernel = "linear", one.lam = FALSE,
                           method = c("direct", "em", "fixed", "canonical",
                                      "mixed"), control = list(),
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                           psi = 1, nystrom = FALSE, nys.seed = NULL,
                           model = list(), ...) {
  # Simply load the kernel and pass to iprior.default() ------------------------
  mod <- kernL2.formula(formula, data, kernel = kernel, one.lam = one.lam,
                        est.lambda = est.lambda, est.hurst = est.hurst,
                        est.lengthscale = est.lengthscale,
                        est.offset = est.offset, est.psi = est.psi,
                        fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
                        nystrom = nystrom, nys.seed = nys.seed, model = model)
  res <- iprior.default(y = mod, method = method, control = control)
  res
}

#' @describeIn iprior Re-run or continue running the EM algorithm from last
#'   attained parameter values in object \code{ipriorMod}.
iprior.ipriorMod <- function(object, method = NULL, control = list(),
                             iter.update = 100, ...) {
  mod           <- object$ipriorKernel
  con           <- object$control
  con$theta0    <- object$theta
  con$maxit     <- iter.update
  control.names <- names(con)
  con[(control.names <- names(control))] <- control
  control <- con

  if (is.null(method)) method <- object$method
  if (method == "em") method.msg <- "EM"
  else method.msg <- method

  message(paste0("Updating iprior model with ", iter.update, " iterations using ",
                 method.msg, " method."))

  # Pass to iprior.default ----------------------------------------------------
  res <- iprior.default(y = mod, method = method, control = con)

  # Update time, call, maxit, niter, lb, error, brier --------------------------
  new.time.diff <- res$end.time - res$start.time
  old.time.diff <- object$end.time - object$start.time
  res$time <- as.time(new.time.diff + old.time.diff)
  res$start.time <- object$start.time
  res$end.time <- object$end.time + new.time.diff
  res$call <- object$call
  res$control$maxit <- iter.update + object$control$maxit
  res$niter <- res$niter + object$niter
  res$loglik <- c(object$loglik, res$loglik)

  res
}
