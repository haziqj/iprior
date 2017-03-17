################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2016  Haziq Jamil
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
iprior <- function(...) {
  # The S3 generic function for objects of class "iprior"
  UseMethod("iprior")
}

#' Fit an I-prior regression model
#'
#' A function to perform linear regression using I-priors. The I-prior model is
#' fitted via maximum likelihood using an EM algorithm.
#'
#' The \code{iprior()} function is able to take formula based input and
#' non-formula. When not using formula, the syntax is as per the default S3
#' method. That is, the response variable is the vector \code{y}, and any
#' explanatory variables should follow this, and separated by commas.
#'
#' As described \link[=kernL]{here}, the model can be loaded first into an
#' \code{ipriorKernel} object, and then passed to the \code{iprior()} function
#' to perform the EM algorithm.
#'
#' If an \code{ipriorMod} object is input, then the EM starts from the last
#' obtained parameter values. This is particularly useful for very heavy models,
#' or models which have not yet converged after reaching the maximum number of
#' iterations. Running \code{iprior()} just continues the EM algorithm.
#'
#' There are several model options available which primarily controls the number
#' and placement of scale parameters \code{lambda} in the model, although these
#' are not applicable when running the function on \code{ipriorMod} or
#' \code{ipriorKernel} objects.
#'
#' @inheritParams kernL
#' @param object This is either an object of class formula (when fitting using
#'   formula interface), \code{ipriorKernel} or \code{ipriorModel}. This is used
#'   when not using formula or \code{"y, x"} input.
#' @param control (optional) A list of control options for the EM algorithm and
#'   output: \describe{\item{\code{maxit}}{The maximum number of iterations
#'   until the EM stops. Defaults to \code{50000}.} \item{\code{stop.crit}.}{The
#'   EM stopping criteria, which is the difference in succesive log-likelihood
#'   values. Defaults to \code{1e-7}.} \item{\code{progress}}{Option for the
#'   reporting of the EM while the function is running. Choose from one of
#'   \code{"lite"} (default), \code{"full"} (log-likelihood and parameters
#'   trace), \code{"predloglik"} (log-likelihood trace only) or \code{"none"}.
#'   Visit the
#'   \href{https://github.com/haziqjamil/iprior/wiki/The-predicted-log-likelihood-feature}{Wiki}
#'   page for more information.} \item{\code{report}}{The EM reports every
#'   \code{report} iterations. Defaults to \code{100}.}
#'   \item{\code{silent}}{(logical) Should the EM report should be printed or
#'   not? This is the same as setting \code{progress = "none"}.}
#'   \item{\code{lambda, psi, sigma}}{These are options to set the initial
#'   values of the parameters. For convenience, the user may choose to input one
#'   of \code{psi} or \code{sigma}, but not both, since \code{psi = 1 / sigma ^
#'   2}.} \item{\code{intercept}}{It is possible to set a fixed value for the
#'   intercept (not recommended).}}
#'
#' @return An object of class \code{ipriorMod} which is a list of 24 items. The
#'   more important items are described below. \describe{ \item{\code{alpha,
#'   lambda, psi, coefficients, sigma}}{The last attained parameter values after
#'   running the EM algorithm. This can also be extracted via \code{coef()}}
#'   \item{\code{log.lik}}{The last attained log-likelihood value. This can also
#'   be extracted via \code{logLik()}.} \item{\code{no.iter}}{The number of
#'   iterations the EM algorithm ran for.} \item{\code{Hlam.mat}}{This is the
#'   scaled kernel matrix of dimension \code{n} by \code{n}.}
#'   \item{\code{VarY.inv}}{The variance-covariance matrix of the marginal
#'   distribution of \code{y}.} \item{\code{w.hat}}{The vector of posterior mean
#'   estimates of the I-prior random effects.} \item{\code{fitted.values}}{These
#'   are posterior estimates of \code{y}, i.e. the fitted values. This can also
#'   be extracted via \code{fitted()} or \code{\link{predict}()}}
#'   \item{\code{residuals}}{The vector of residuals. This can also be extracted
#'   via \code{resid()}.} }
#'
#' @examples
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
#' mod.toothnf <- iprior(y = ToothGrowth$len,
#'                       supp = ToothGrowth$supp,
#'                       dose = ToothGrowth$dose,
#'                       model = list(interactions = "1:2"))
#'
#' # Formula based model option one.lam = TRUE
#' # Sets a single scale parameter for all variables
#' modf <- iprior(stack.loss ~ ., data = stackloss, model = list(one.lam = TRUE))
#' modnf <- iprior(y = stackloss$stack.loss, x = stackloss[1:3])
#'
#' # Example of using the FBM kernel for smoothing models
#' mod <- kernL(y ~ x, datfbm, model = list(kernel = "FBM"))  # Hurst = 0.5 (default)
#' mod <- kernL(y ~ x, datfbm, model = list(kernel = "FBM,0.75"))  # custom Hurst
#'
#' # Fit the model using EM starting at a specific parameter value
#' mod.fit <- iprior(mod, control = list(lambda = 8.41, psi = 0.33))
#'
#' @name iprior
#' @export
iprior.default <- function(y, ..., model = list(), control = list()) {
  # Set up the controls for the EM algorithm -----------------------------------
  con <- list(maxit = 50000, stop.crit = 1e-07, report = 100, intercept = NULL,
              lambda = NULL, psi = NULL, sigma = NULL, theta = NULL,
              progress = "lite", silent = FALSE, force.regEM = FALSE,
              force.nlm = FALSE)
  con_names <- names(con)
  con[(control_names <- names(control))] <- control
  if (length(noNms <- control_names[!control_names %in% con_names])) {
    warning("Unknown names in control options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  list2env(con, environment())
  silent_ <- silent
  progress <- match.arg(progress, c("lite", "none", "full", "predloglik"))
  if (progress == "lite") {
    clean         <- TRUE
    silent        <- FALSE
    paramprogress <- FALSE
  }
  if (progress == "none" | silent) {
    clean         <- TRUE
    silent        <- TRUE
    paramprogress <- FALSE
  }
  if (progress == "full") {
    clean         <- FALSE
    silent        <- FALSE
    paramprogress <- TRUE
  }
  if (progress == "predloglik") {
    clean         <- FALSE
    silent        <- FALSE
    paramprogress <- FALSE
  }
  if (silent_) silent <- silent_

  # Check initial values for parameters ----------------------------------------
  par.check1 <- !is.null(theta) & any(c(!is.null(lambda), !is.null(psi),
                                        !is.null(sigma)) )
  if (par.check1) {
    stop("Starting values stated for both theta and lambda/psi/sigma.",
         call. = FALSE)
  }
  par.check2 <- !is.null(psi) & !is.null(sigma)
  if (par.check2) {
    stop("Only one of psi or sigma can be initiated.", call. = FALSE)
  }
  if (!is.null(theta)) {
    psi <- theta[length(theta)]
    lambda <- theta[-length(theta)]
  }
  if (!is.null(sigma)) psi <- 1 / sqrt(sigma)

  # Set yname ------------------------------------------------------------------
  cl <- match.call()
  ynamefromcall <- as.character(cl[2])
  check.yname <- is.null(model$yname)
  if (check.yname) model$yname <- ynamefromcall

  # Pass to kernel loader and then EM routine ----------------------------------
  if (is.ipriorKernel(y)) {
    ipriorKernel <- y
  } else {
    # When using y and x to call iprior
    ipriorKernel <- kernL(y, ..., model = model)  # pass to kernel loader
  }
  est <- ipriorEM(ipriorKernel, maxit, stop.crit, report, silent, intercept,
                  lambda, psi, clean, paramprogress, force.regEM, force.nlm)
  est$ipriorKernel <- ipriorKernel
  est$sigma <- 1/sqrt(est$psi)
  if (ipriorKernel$model$rootkern) {
    # Do if GPR
    param <- c(est$alpha, est$lambda ^ 2 * est$psi, est$sigma)
    psi.or.sigma <- "sigma"
  } else {
    param <- c(est$alpha, est$lambda, est$psi)
    psi.or.sigma <- "psi"
  }
  if (length(param) == 3) {
    names(param) <- c("(Intercept)", "lambda", psi.or.sigma)
  } else {
    names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)),
                      psi.or.sigma)
  }

  # Fix xname ------------------------------------------------------------------
  mx <- match(c("y", "model", "control"), names(cl), 0L)
  xnamefromcall <- as.character(cl[-mx])[-1]
  check.xname <- grepl("\\.\\.", est$ipriorKernel$model$xname)
  if (any(check.xname)) {
    est$ipriorKernel$model$xname[check.xname] <- xnamefromcall[check.xname]
    est$ipriorKernel$model$lamnamesx <- est$ipriorKernel$model$xname[whereOrd(est$ipriorKernel$model$order)]
  }

  # Calculate fitted values and residuals --------------------------------------
  est$residuals     <- ipriorKernel$Y - est$fitted.values
  names(est$fitted.values) <- names(est$residuals) <- names(ipriorKernel$Y)

  # Changing the call to simply iprior -----------------------------------------
  est$fullcall <- cl
  cl[[1L]] <- as.name("iprior")
  m <- match(c("control"), names(cl), 0L)
  if (any(m > 0)) cl <- cl[-m]
  est$call <- cl

  # Other things to return -----------------------------------------------------
  est$control      <- con
  est$coefficients <- param
  est$T2           <- as.numeric(crossprod(est$w.hat)/est$psi)

  # Gaussian Process Regression estimates --------------------------------------


  class(est) <- "ipriorMod"
  est
}

#' @rdname iprior
#' @export
iprior.formula <- function(formula, data = parent.frame(), model = list(),
                           control = list(), ...) {
  # Pass to iprior default -----------------------------------------------------
  ipriorKernel <- kernL(formula, data, model = model)
  est <- iprior.default(y = ipriorKernel, control = control)

  # Changing the call to simply iprior -----------------------------------------
  cl <- match.call()
  est$fullcall <- cl
  cl[[1L]] <- as.name("iprior")
  m <- match(c("formula", "data"), names(cl), 0L)
  cl <- cl[c(1L, m)]
  est$call <- cl
  names(est$call)[2] <- "formula"
  est$formula <- formula
  est$terms <- class(est) <- "ipriorMod"
  est
}

#' @describeIn iprior Takes in object of type \code{ipriorKernel} and estimates
#'   the parameters of the model via the EM algorithm.
#' @export
iprior.ipriorKernel <- function(object, control = list(), ...) {
  est <- iprior.default(y = object, control = control)

  # Fix the call ---------------------------------------------------------------
  cl <- est$ipriorKernel$call
  cl[[1L]] <- as.name("iprior")
  est$call <- cl

  est
}

#' @describeIn iprior Re-run or continue running the EM algorithm from last
#'   attained parameter values in object \code{ipriorMod}.
#' @export
iprior.ipriorMod <- function(object, control = list(), ...) {
  ipriorKernel <- object$ipriorKernel
  lambda       <- object$lambda
  psi          <- object$psi

  con <- list(lambda = lambda, psi = psi)
  con <- c(control, con)

  est <- iprior.default(y = ipriorKernel, control = con)
  est$fullcall <- object$fullcall
  assign(deparse(substitute(object)), est, envir = parent.frame())
}

#' @export
print.ipriorMod <- function(x, ...) {
  kernel <- x$ipriorKernel$model$kernel
  Hurst <- x$ipriorKernel$model$Hurst
  cat("\nCall:\n")
  print(x$call)
  if (length(unique(Hurst)) > 1) {
    FBM <- "Fractional Brownian Motion with multiple Hurst coef."
  } else {
    FBM <- paste0("Fractional Brownian Motion with Hurst coef. ",
                  unique(Hurst))
  }

  kerneltypes <- c("Pearson", "Canonical", FBM)
  which.kern <- c("Pearson", "Canonical", "FBM") %in% unique(kernel)
  #                  "Pearson & Canonical", paste("Pearson &", FBM),
  #                  paste("Canonical &", FBM),
  #                  paste("Pearson, Canonical, &", FBM))
  if (x$ipriorKernel$model$rootkern) {
    if (sum(which.kern) == 1) {
      cat("\nGPR with", kerneltypes[which.kern])
    } else if (sum(which.kern) == 2) {
      cat("\nGPR with", paste(kerneltypes[which.kern], collapse = " & "))
    } else {
      cat("\nGPR with", paste("Pearson, Canonical, &", FBM))
    }
    cat(" covariance kernel.\n")
  } else {
    if (sum(which.kern) == 1) {
      cat("\nRKHS used:", kerneltypes[which.kern])
    } else if (sum(which.kern) == 2) {
      cat("\nRKHS used:", paste(kerneltypes[which.kern], collapse = " & "))
    } else {
      cat("\nRKHS used:", paste("Pearson, Canonical, &", FBM))
    }
    if (x$ipriorKernel$l == 1) {
      cat(", with a single scale parameter.\n")
    } else {
      cat(", with multiple scale parameters.\n")
    }
  }
  # cat("\n")
  cat("\nParameter estimates:\n")
  print(x$coefficients)
  cat("\n")
}

#' @export
summary.ipriorMod <- function(object, ...) {
  # Standard errors from inverse observed Fisher matrix ------------------------
  se <- fisher(object)
  if (object$ipriorKernel$model$rootkern) {
    se[c(-1, -length(se))] <- se[c(-1, -length(se))] * 2 * object$lambda * object$psi
  }

  # Z values to compare against (standard) Normal distribution -----------------
  zval <- coef(object)/se

  # Create table for summary output --------------------------------------------
  tab <- cbind(Estimate   = round(coef(object), digits = 4),
               S.E.       = round(se, digits = 4),
               z          = round(zval, digits = 3),
               `P[|Z>z|]` = round(2 * pnorm(-abs(zval)), digits = 3))
  xname <- object$ipriorKernel$model$xname
  if (object$ipriorKernel$l == 1) {
    # only rename rows when using multiple lambdas
    lamnames <- c("(Intercept)", "lambda", "psi")
    rownames(tab) <- lamnames
  } else {
    lamnames <- paste0("lam", 1:(length(coef(object)) - 2))
    lamnames <- c("(Intercept)", paste(lamnames,
                                       object$ipriorKernel$model$lamnamesx,
                                       sep = "."), "psi")
    rownames(tab) <- lamnames
  }
  tab <- tab[-length(coef(object)), ]  # removes the psi from the table

  res <- list(call = object$call, coefficients = tab,
              whichPearson = object$ipriorKernel$whichPearson,
              kernel = object$ipriorKernel$model$kernel,
              resid = object$residuals, log.lik = object$log.lik,
              no.iter = object$no.iter, converged = object$converged,
              stop.crit = object$control$stop.crit,
              one.lam = object$ipriorKernel$model$one.lam, T2 = object$T2,
              l = object$ipriorKernel$l, p = object$ipriorKernel$p,
              Hurst = object$ipriorKernel$model$Hurst,formula = object$formula,
              psi.and.se = c(coef(object)[length(se)], se[length(se)]),
              xname = xname, no.int = object$ipriorKernel$no.int,
              optim.converged = object$optim.converged,
              rootkern = object$ipriorKernel$model$rootkern)
  class(res) <- "ipriorSummary"
  res
}

#' @export
print.ipriorSummary <- function(x, ...) {
  # The print out of the S3 summary method for iprior objects.
  cat("\nCall:\n")
  print(x$call)
  x.names <- x$xname[1:x$p]
  x.pea <- x.names[isPea(x$kernel)]
  x.can <- x.names[isCan(x$kernel)]
  x.fbm <- x.names[isFBM(x$kernel)]
  printPea <- paste0("Pearson (", paste(x.pea, collapse = ", "), ")")
  printCan <- paste0("Canonical (", paste(x.can, collapse = ", "), ")")

  cat("\n")
  if (x$rootkern) cat("GPR covariance kernel:\n")
  else cat("RKHS used:\n")
  if (!(length(x.pea) == 0)) cat(printPea, "\n")
  if (!(length(x.can) == 0)) cat(printCan, "\n")
  if (!(length(x.fbm) == 0)) {
    Hurst <- x$Hurst[isFBM(x$kernel)]
    un.Hurst <- unique(Hurst)
    uH <- length(un.Hurst)
    printFBM <- list(NULL)
    for (i in 1:uH) {
      printFBM[[i]] <- paste0("Fractional Brownian Motion with Hurst coef. ",
                              un.Hurst[i], " (", paste(x.fbm[Hurst == un.Hurst[i]],
                                                       collapse = ", "), ")")
    }
    printFBM <- paste(printFBM, collapse = "\n")
    cat(printFBM, "\n")
  }
  if (!x$rootkern) {
    if (x$l == 1) {
      cat("with a single scale parameter.\n")
    } else {
      cat("with multiple scale parameters.\n")
    }
  }
  cat("\n")
  cat("Residuals:\n")
  print(summary(x$resid)[-4])
  cat("\n")
  tab <- x$coefficients
  psi.and.se <- x$psi.and.se
  sigma <- 1/sqrt(psi.and.se[1])
  sesigma <- psi.and.se[2] * sigma ^ 3 / 2
  printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  if (!is.null(x$optim.converged)) {
    cat("Routine converged via EM and direct optimisation.")
  } else {
    if (x$converged) {
      cat("EM converged to within", x$stop.crit, "tolerance.")
    } else {
      cat("EM failed to converge.")
    }
    cat(" No. of iterations:", x$no.iter)
  }
  cat("\nStandard deviation of errors:", signif(sigma, digits = 4),
      "with S.E.:", round(sesigma, digits = 4))
  # cat("\nT2 statistic:", signif(x$T2, digits = 4), "on ??? degrees of freedom.")
  cat("\nLog-likelihood value:", x$log.lik, "\n")
  cat("\n")
}
