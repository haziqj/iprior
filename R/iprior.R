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
iprior <- function(...) UseMethod("iprior")
# iprior <- function(y, ..., kernel = "linear", method = "direct",
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
#' As described \link[=kernL]{here}, the model can be loaded first into an
#' \code{ipriorKernel} object, and then passed to the \code{iprior()} function
#' to perform the estimation.
#'
#' @inheritParams kernL
#' @param object An \code{ipriorKernel} or \code{ipriorMod} object.
#' @param method The estimation method. One of: \itemize{ \item{\code{"direct"}
#'   - for the direct minimisation of the marginal deviance using
#'   \code{optim()}'s L-BFGS method} \item{\code{"em"} - for the EM algorithm}
#'   \item{\code{"mixed"} - combination of the direct and EM methods}
#'   \item{\code{"fixed"} - for just obtaining the posterior regression function
#'   with fixed hyperparameters (default method when setting \code{fixed.hyp =
#'   TRUE})} \item{\code{"canonical"} - an efficient estimation method which
#'   takes advantage of the structure of the linear kernel} }
#' @param control (Optional) A list of control options for the estimation
#'   procedure: \describe{ \item{\code{maxit}}{The maximum number of iterations
#'   for the quasi-Newton optimisation or the EM algorithm. Defaults to
#'   \code{100}.} \item{\code{em.maxit}}{For \code{method = "mixed"}, the number
#'   of EM steps before switching to direct optimisation. Defaults to \code{5}.}
#'   \item{\code{stop.crit}}{The stopping criterion for the EM and L-BFGS
#'   algorithm, which is the difference in successive log-likelihood values.
#'   Defaults to \code{1e-8}.} \item{\code{theta0}}{The initial values for the
#'   hyperparameters. Defaults to random starting values.}
#'   \item{\code{report}}{The interval of reporting for the \code{optim()}
#'   function.} \item{\code{restarts}}{The number of random restarts to perform.
#'   Defaults to \code{0}. It's also possible to set it to \code{TRUE}, in which
#'   case the number of random restarts is set to the total number of available
#'   cores.} \item{\code{no.cores}}{The number of cores in which to do random
#'   restarts. Defaults to the total number of available cores.}
#'   \item{\code{omega}}{The overrelaxation parameter for the EM algorithm - a
#'   value between 0 and 1.}}
#' @param iter.update The number of iterations to perform when calling the
#'   function on an \code{ipriorMod} object. Defaults to \code{100}.
#'
#' @return An \code{ipriorMod} object. Several accessor functions have been
#'   written to obtain pertinent things from the \code{ipriorMod} object. The
#'   \code{print()} and \code{summary()} methods display the relevant model
#'   information.
#'
#' @seealso \link[=optim]{optim}, \link[=update.ipriorMod]{update},
#'   \link[=check_theta]{check_theta}, print, summary, plot, coef, sigma,
#'   fitted, predict, logLik, deviance.
#'
#' @examples
#'
#' # Formula based input
#' (mod.stackf <- iprior(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
#'                       data = stackloss))
#' mod.toothf <- iprior(len ~ supp * dose, data = ToothGrowth)
#' summary(mod.toothf)
#'
#' # Non-formula based input
#' mod.stacknf <- iprior(y = stackloss$stack.loss,
#'                       Air.Flow = stackloss$Air.Flow,
#'                       Water.Temp = stackloss$Water.Temp,
#'                       Acid.Conc. = stackloss$Acid.Conc.)
#' mod.toothnf <- iprior(y = ToothGrowth$len, ToothGrowth$supp, ToothGrowth$dose,
#'                       interactions = "1:2")
#'
#' # Formula based model option one.lam = TRUE
#' # Sets a single scale parameter for all variables
#' modf <- iprior(stack.loss ~ ., data = stackloss, one.lam = TRUE)
#' modnf <- iprior(y = stackloss$stack.loss, X = stackloss[1:3])
#' all.equal(coef(modnf), coef(modnf))  # both models are equivalent
#'
#' # Fit models using different kernels
#' dat <- gen_smooth(n = 100)
#' mod <- iprior(y ~ X, dat, kernel = "fbm")  # Hurst = 0.5 (default)
#' mod <- iprior(y ~ X, dat, kernel = "poly3")  # polynomial degree 3
#'
#' # Fit models using various estimation methods
#' mod1 <- iprior(y ~ X, dat)
#' mod2 <- iprior(y ~ X, dat, method = "em")
#' mod3 <- iprior(y ~ X, dat, method = "canonical")
#' mod4 <- iprior(y ~ X, dat, method = "mixed")
#' mod5 <- iprior(y ~ X, dat, method = "fixed", lambda = coef(mod1)[1],
#'                psi = coef(mod1)[2])
#' c(logLik(mod1), logLik(mod2), logLik(mod3), logLik(mod4),
#'   logLik(mod5))
#'
#' \dontrun{
#'
#' # For large data sets, it is worth trying the Nystrom method
#' mod <- iprior(y ~ X, gen_smooth(5000), kernel = "se", nystrom = 50,
#'               est.lengthscale = TRUE)  # a bit slow
#' plot_fitted(mod, ci = FALSE)
#' }
#'
#' @name iprior
#' @export
iprior.default <- function(y, ..., kernel = "linear", method = "direct",
                           control = list(), interactions = NULL,
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                           psi = 1, nystrom = FALSE, nys.seed = NULL,
                           model = list(), train.samp, test.samp) {
  # Load the I-prior model -----------------------------------------------------
  if (is.ipriorKernel(y)) {
    mod <- y
  } else {
    mod <- kernL(y = y, ..., kernel = kernel, interactions = interactions,
                  est.lambda = est.lambda, est.hurst = est.hurst,
                  est.lengthscale = est.lengthscale, est.offset = est.offset,
                  est.psi = est.psi, fixed.hyp = fixed.hyp, lambda = lambda,
                  psi = psi, nystrom = nystrom, nys.seed = nys.seed,
                  model = model, train.samp = train.samp, test.samp = test.samp)
  }
  if (is.categorical(mod)) {
    warning("Categorical responses loaded. Consider using iprobit package.",
            call. = FALSE)
  }

  # Set up controls ------------------------------------------------------------
  method <- tolower(method)
  method <- match.arg(method, c("direct", "em", "fixed", "canonical", "mixed"))

  control_ <- list(
    maxit     = 100,
    em.maxit  = 5,  # for mixed method
    par.maxit = 5,  # for parallel method
    stop.crit = 1e-8,
    theta0    = NULL,
    # lambda0   = NULL,
    # psi0      = NULL,
    silent    = FALSE,
    report    = 10,
    psi.reg   = FALSE,  # option for iprior_em_reg()
    restarts  = 0,
    no.cores  = parallel::detectCores(),
    optim.method = "L-BFGS",
    omega     = 0
  )
  control <- update_control(control, control_)  # see iprior_helper.R
  control.optim <- list(
    fnscale = -2,
    trace   = ifelse(isTRUE(control$silent), 0, 1),
    maxit   = max(0, control$maxit - 1),
    REPORT  = control$report,
    factr   = control$stop.crit / .Machine$double.eps
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
    if (est.method["fixed"]) {  # FIXED METHOD
      res <- iprior_fixed(mod)
      res$est.method <- "Estimation with fixed hyperparameters."
      res$est.conv <- "Convergence not assessed."
    } else if (est.method["canonical"]) {
      res <- iprior_canonical(mod, theta0, control.optim)
      res$est.method <- "Efficient canonical method."
    } else {
      if (est.method["em.closed"]) {  # EM CLOSED-FORM
        res <- iprior_em_closed(mod, control$maxit, control$stop.crit,
                                control$silent, theta0, omega = control$omega)
        res$est.method <- "Closed-form EM algorithm."
      }
      if (est.method["em.reg"]) {  # EM REGULAR
        res <- iprior_em_reg(mod, control$maxit, control$stop.crit,
                             control$silent, theta0)
        res$est.method <- "Regular EM algorithm."
      }
      if (est.method["direct"]) {  # DIRECT OPTIMISATION
        res <- iprior_direct(mod, loglik_iprior, theta0, control.optim,
                             control$optim.method)
        res$est.method <- "Direct optimisation method."
      }
      if (est.method["nystrom"]) {  # NYSTROM DIRECT OPTIMISATION
        res <- iprior_direct(mod, loglik_nystrom, theta0, control.optim,
                             control$optim.method)
        res$est.method <- "Nystrom approximated optimisation."
      }
      if (est.method["mixed"]) {  # MIXED METHOD
        res <- iprior_mixed(mod, theta0, control$em.maxit, control$stop.crit,
                            control$silent, control.optim, control$optim.method)
        res$est.method <- paste0("EM algorithm (", control$em.maxit,
                                 " steps) + direct minimisation.")
      }
      if (res$conv == 0)
        res$est.conv <- paste0("Converged to within ", control$stop.crit,
                               " tolerance.")
      else if (res$conv == 1)
        res$est.conv <- "Convergence criterion not met."
      else
        res$est.conv <- res$message
    }
  }
  res$ipriorKernel <- mod
  res$coefficients <- reduce_theta(res$param.full, mod$estl)$theta.reduced

  # Fitted and predicted values ------------------------------------------------
  tmp <- mse_iprior(mod$y, y.hat = res$y.hat)  # no intercept added
  res$fitted.values <- tmp$y + get_intercept(mod)
  names(res$fitted.values) <- attr(mod$y, "dimnames")[[1]]
  res$residuals <- tmp$resid
  res$train.error <- tmp$train.error
  if (!is.null(mod$y.test) & !is.null(mod$Xl.test)) {
    res$test <- structure(predict_iprior(res, mod$Xl.test, mod$y.test),
                          class = "ipriorPredict")
  }

  # Other stuff ----------------------------------------------------------------
  res$method <- method
  res$control <- control
  # res$fullcall <- match.call()
  res$call <- fix_call_default(match.call(), "iprior")
  res$ipriorKernel$call <- fix_call_default(match.call(), "kernL")

  class(res) <- "ipriorMod"
  res
}

#' @rdname iprior
#' @export
iprior.formula <- function(formula, data, kernel = "linear", one.lam = FALSE,
                           method = "direct", control = list(),
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                           psi = 1, nystrom = FALSE, nys.seed = NULL,
                           model = list(), train.samp, test.samp, ...) {
  # Simply load the kernel and pass to iprior.default() ------------------------
  mod <- kernL.formula(formula, data, kernel = kernel, one.lam = one.lam,
                       est.lambda = est.lambda, est.hurst = est.hurst,
                       est.lengthscale = est.lengthscale,
                       est.offset = est.offset, est.psi = est.psi,
                       fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
                       nystrom = nystrom, nys.seed = nys.seed, model = model,
                       train.samp = train.samp, test.samp = test.samp, ...)
  res <- iprior.default(y = mod, method = method, control = control)
  res$call <- fix_call_formula(match.call(), "iprior")
  res$ipriorKernel$call <- fix_call_formula(match.call(), "kernL")
  res
}

#' @describeIn iprior Takes in object of type \code{ipriorKernel}, a loaded and
#'   prepared I-prior model, and proceeds to estimate it.
#' @export
iprior.ipriorKernel <- function(object, method = "direct",
                                 control = list(), ...) {
  res <- iprior.default(y = object, method = method, control = control)

  # Fix call -------------------------------------------------------------------
  res$object$call <- ipriorKernel.call <- object$call
  if (is.null(object$formula)) {
    res$call <- fix_call_default(ipriorKernel.call, "iprior")
  } else {
    res$call <- fix_call_formula(ipriorKernel.call, "iprior")
  }

  res
}

#' @describeIn iprior Re-run or continue running the EM algorithm from last
#'   attained parameter values in object \code{ipriorMod}.
#' @export
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
