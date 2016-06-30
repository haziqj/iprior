#' Fit an I-prior regression model
#'
#' Description
#'
#' Details
#'
#' @param object This is either an object of class formula (when fitting using
#'   formula interface), \code{ipriorKernel} or \code{ipriorModel}. When left
#'   \code{NULL}, then \code{y} must be used, along with the
#'   covariates/independent variables.
#' @param y Vector of response variables.
#' @param ... Only for when fitting using non-formula, enter the variables
#'   (vectors or matrices) separated by commas. No other options applicable
#'   here.
#' @param model (optional) List of model options. Not used for
#'   \code{ipriorKernel} or \code{ipriorModel} objects.
#' @param control (optional) list of control options for EM algorithm and
#'   output.
#' @param data Data frame containing variables when using formula interface.
#'
#' @return An object of class \code{ipriorMod}.
#'
#' @examples (mod.iprior <- iprior(stack.loss ~ ., data = stackloss))
#'
#' @name iprior
#' @export
iprior <- function(object, ...) {
  # The S3 generic function for objects of class "iprior"
  UseMethod("iprior")
}

#' @rdname iprior
#' @export
iprior.default <- function(object = NULL, model = list(),
                           control = list(), y, ...) {
  # Set up the controls for the EM algorithm -----------------------------------
  con <- list(maxit = 50000, stop.crit = 1e-07, report.int = 100, lambda = NULL,
              psi = abs(rnorm(1)), progress = "lite", silent = FALSE)
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
  cl <- match.call()

  # Pass to kernel loader and then EM routine ----------------------------------
  if (!is.null(y)) {
    # When using y and x to call iprior
    ipriorKernel <- kernL(y, ..., model = model)  # pass to kernel loader
  } else {
    # Otherwise, it is guaranteed to be an ipriorKernel object via the methods
    ipriorKernel <- object
  }
  est <- ipriorEM(ipriorKernel, maxit, stop.crit, report.int, silent_, lambda,
                  psi, clean, paramprogress)
  est$ipriorKernel <- ipriorKernel
  param <- c(est$alpha, est$lambda, est$psi)
  if (length(param) == 3) {
    names(param) <- c("(Intercept)", "lambda", "psi")
  } else {
    names(param) <- c("(Intercept)", paste0("lambda", 1:length(est$lambda)),
                      "psi")
  }

  # Calculate fitted values and residuals --------------------------------------
  if (maxit == 0) {
    Y.hat <- rep(est$alpha, nrow(est$H.mat.lam))
  } else {
    Y.hat <- est$alpha + as.vector(crossprod(est$H.mat.lam, est$w.hat))
  }
  est$fitted.values <- Y.hat
  est$residuals     <- ipriorKernel$Y - Y.hat
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
  est$sigma        <- 1/sqrt(est$psi)
  est$T2           <- as.numeric(crossprod(est$w.hat)/est$psi)

  class(est) <- "ipriorMod"
  est
}

#' @rdname iprior
#' @export
iprior.formula <- function(object, data, model = list(), control = list(), ...) {
  # Formula based S3 constructor function for iprior.

  # Pass to iprior default -----------------------------------------------------
  ipriorKernel <- kernL(object, data, model = model)
  est <- iprior.default(object = ipriorKernel, control = control, y = NULL)

  # Changing the call to simply iprior -----------------------------------------
  cl <- match.call()
  est$fullcall <- cl
  cl[[1L]] <- as.name("iprior")
  m <- match(c("object", "data"), names(cl), 0L)
  cl <- cl[c(1L, m)]
  est$call <- cl
  names(est$call)[2] <- "formula"
  est$formula <- object
  est$terms <- class(est) <- "ipriorMod"
  est
}

#' @describeIn iprior Takes in object of type \code{ipriorKernel} and estimates
#'   the parameters of the model via the EM algorithm.
#' @export
iprior.ipriorKernel <- function(object, control = list(), ...) {
  est <- iprior.default(object, control = control, y = NULL)
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

  est <- iprior.default(object = ipriorKernel, control = con, y = NULL)
  est$fullcall <- object$fullcall
  assign(deparse(substitute(object)), est, envir = parent.frame())
}

#' @export
print.ipriorMod <- function(x, ...) {
  whichPearson <- x$ipriorKernel$whichPearson
  cat("\nCall:\n")
  print(x$call)
  if (x$ipriorKernel$model$kernel == "Canonical") {
    CanOrFBM <- "Canonical"
  } else {
    CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ",
                       x$ipriorKernel$model$Hurst)
  }
  kerneltypes <- c(CanOrFBM, "Pearson", paste(CanOrFBM, "& Pearson"))
  if (all(x$ipriorKernel$whichPearson)) {
    cat("\nRKHS used:", kerneltypes[2])
  } else {
    if (!all(whichPearson) && !any(whichPearson)) {
      cat("\nRKHS used:", kerneltypes[1])
    } else {
      cat("\nRKHS used:", kerneltypes[3])
    }
  }
  if (x$ipriorKernel$q == 1) {
    cat(", with a single scale parameter.\n")
  } else {
    cat(", with multiple scale parameters.\n")
  }
  cat("\n")
  cat("\nParameter estimates:\n")
  print(x$coefficients)
  cat("\n")
}

#' Summary screen for iprior models.
#'
#' Summary screen for iprior models.
#'
#' Summary screen for iprior models.
#'
#' @param object Objects of class \code{iprior}.
#'
#' @examples
#' \donttest{mod <- iprior(Hwt ~ ., data=MASS::cats)}
#' \donttest{summary(mod)}
#'
#'
#' @export
summary.ipriorMod <- function(object, ...) {
  # Standard errors from inverse observed Fisher matrix ------------------------
  se <- fisher(alpha = object$alpha, psi = object$psi, lambda = object$lambda,
                  P.matsq = object$P.matsq, H.mat.lam = object$H.mat.lam,
                  S.mat = object$S.mat, Var.Y.inv = object$Var.Y.inv)

  # Z values to compare against (standard) Normal distribution -----------------
  zval <- coef(object)/se

  # Create table for summary output --------------------------------------------
  tab <- cbind(Estimate   = round(coef(object), digits = 4),
               S.E.       = round(se, digits = 4),
               z          = round(zval, digits = 3),
               `P[|Z>z|]` = round(2 * pnorm(-abs(zval)), digits = 3))
  xname <- object$ipriorKernel$model$xname
  if (object$ipriorKernel$q == 1) {
    # only rename rows when using multiple lambdas
    lamnames <- c("(Intercept)", "lambda", "psi")
    rownames(tab) <- lamnames
  } else {
    lamnames <- paste0("lam", 1:(length(coef(object)) - 2))
    lamnames <- c("(Intercept)", paste(lamnames, xname[1:object$ipriorKernel$q],
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
              q = object$ipriorKernel$q, p = object$ipriorKernel$p,
              Hurst = object$ipriorKernel$model$Hurst,formula = object$formula,
              psi.and.se = c(coef(object)[length(se)], se[length(se)]),
              xname = xname)
  class(res) <- "ipriorSummary"
  res
}

#' @export
print.ipriorSummary <- function(x, ...) {
  # The print out of the S3 summary method for iprior objects.
  cat("\nCall:\n")
  print(x$call)
  x.names <- x$xname[1:x$q]
  xPearson <- x.names[x$whichPearson]
  xCanOrFBM <- x.names[!x$whichPearson]
  if (x$kernel == "Canonical") {
    CanOrFBM <- "Canonical"
  } else {
    CanOrFBM <- paste0("Fractional Brownian Motion with Hurst coef. ", x$Hurst)
  }
  printPearson <- paste0("Pearson (", paste(xPearson, collapse = ", "), ")")
  printCanOrFBM <- paste0(CanOrFBM, " (", paste(xCanOrFBM, collapse = ", "), ")")
  cat("\n")
  cat("RKHS used:\n")
  if (!(length(xCanOrFBM) == 0)) cat(printCanOrFBM, "\n")
  if (!(length(xPearson) == 0)) cat(printPearson, "\n")
  if (x$q == 1) {
    cat("with a single scale parameter.\n")
  } else {
    cat("with multiple scale parameters.\n")
  }
  cat("\n")
  cat("Residuals:\n")
  print(summary(x$resid)[-4])
  cat("\n")
  tab <- x$coefficients
  psi.and.se <- x$psi.and.se
  sigma <- 1/sqrt(psi.and.se[1])
  sesigma <- psi.and.se[2] * sigma ^ 3 / 2
  printCoefmat(tab, P.value = T, has.Pvalue = T)
  cat("\n")
  if (x$converged) {
    cat("EM converged to within", x$stop.crit, "tolerance.")
  } else {
    cat("EM failed to converge.")
  }
  cat(" No. of iterations:", x$no.iter)
  cat("\nStandard deviation of errors:", signif(sigma, digits = 4),
      "with S.E.:", round(sesigma, digits = 4))
  cat("\nT2 statistic:", signif(x$T2, digits = 4), "on ??? degrees of freedom.")
  cat("\nLog-likelihood value:", x$log.lik, "\n")
  cat("\n")
}
