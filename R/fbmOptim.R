#' Find the Hurst coefficient of a FBM I-prior model
#'
#' From an \code{ipriorKernel} object, a golden section search of (0,1) is
#' performed using \code{optimize()}.
#'
#' Currently, \code{fbmOptim()} will only search for a single optimum Hurst
#' coefficient. If the model specified has multiple Hurst coefficients, this
#' will be overridden with the single optimal Hurst value.
#'
#' @param object An object of class \code{ipriorKernel}.
#' @param method One of \code{c("ipriorOptim", "iprior")} for model fitting of
#'   the final I-prior model.
#' @param silent (logical) Run the optimisation silently or not.
#'
#' @return An \code{ipriorMod} object.
#' @export
#'
#' @examples
#' mod <- kernL(y ~ ., datfbm, model = list(kernel = "FBM"))
#' (mod.iprior <- fbmOptim(mod))
fbmOptim <- function(object, method = c("ipriorOptim", "iprior"),
                     silent = FALSE) {
  if (!is.ipriorKernel(object)) {
    stop("Input objects of class ipriorKernel only.", call. = FALSE)
  }
  if (!any(isFBM(object$model$kernel))) {
    stop("This only works if at least one of the kernels is FBM.",
         call. = FALSE)
  }
  Hurst <- object$model$Hurst[isFBM(object$model$kernel)]
  un.Hurst <- unique(Hurst)
  if (length(un.Hurst) > 1) {
    warning("Currently, fbmOptim() will search for a single optimum Hurst coefficient.",
         call. = FALSE)
  }

  method <- match.arg(method,  c("ipriorOptim", "iprior"))

  res <- stats::optimise(fbmOptimDeviance, c(0, 1),
                         object = object, silent = silent)
  update.ipriorKernel(object, round(res$min, 5))

  if (!silent) cat("Optimum Hurst coefficient found.\n")
  if (!silent) cat("\nPreparing iprior output... ")
  if (method == "iprior") {
    mod.fit <- iprior(object, control = list(silent = TRUE))
  }
  if (method == "ipriorOptim") {
    mod.fit <- ipriorOptim(object, control = list(silent = TRUE))
  }

  if (!silent) cat("DONE.\n")

  mod.fit
}

fbmOptimDeviance <- function(gamma, object, silent = FALSE) {
  # Returns deviance of an ipriorKernel object for a particular Hurst coefficient
  # gamma. Used to find optimum value of gamma in fbmOptim().
  if (!silent) cat("Hurst = ", gamma, "\n")
  update.ipriorKernel(object, gamma)
  mod.fit <- ipriorOptim(object, control = list(silent = TRUE))
  return(deviance.ipriorMod(mod.fit))
}

# mod <- kernL(y ~., datfbm, model = list(kernel = "FBM"))
