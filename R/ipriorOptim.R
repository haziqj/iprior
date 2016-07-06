#' Estimate an I-prior model using a combination of EM algorithm and direct
#' optimisation
#'
#' This is a wrapper function for \code{iprior()} and \code{optim()} which
#' estimates an I-prior model that has been stored in an \code{ipriorKernel}
#' object.
#'
#' The EM algorithm is slow to converge at times, but every iteration is
#' guaranteed to increase the likelihood value. On the other hand a direct
#' maximisation of the I-prior likelihood may sometimes result in
#' ill-conditioned variance parameter due to the nature of the parameterisation
#' of the I-prior model. Thus, an ideal implementation is a combination of EM
#' and direct optimisation.
#'
#' First, the EM algorithm is performed for a maximum of five iterations. The
#' parameters are then passed to \code{optim()} and the negative log-likelihood
#' is minimised. The method used for optim is \code{"L-BFGS-B"}, as the
#' \code{psi} parameter of the I-prior model needs to be contrained to be
#' greater than zero.
#'
#' @param object An object of class \code{ipriorKernel}.
#' @param control A list of controls for the initial EM algorithm fit. Refer to
#'   \code{\link{iprior}} for a full list of available controls.
#'
#' @return An object of class \code{ipriorMod}.
#'
#' @seealso \code{\link{iprior}} and \code{\link{kernL}}.
#'
#' @examples
#' (mod <- kernL(stack.loss ~ ., stackloss))
#' mod.iprior <- ipriorOptim(mod)
#' summary(mod.iprior)
#'
#' @export
ipriorOptim <- function(object, control = list()) {
  if (!is.ipriorKernel(object)) {
    stop("Input objects of class ipriorKernel only.", call. = FALSE)
  }

  control$maxit <- 5
  control$report <- 1
  mod.iprior <- iprior(object, control = control)
  silent <- mod.iprior$control$silent

  if (!silent) cat("\nNow switching to optim...\n\n")
  mod.optim <- optim(par = mod.iprior$coef[-1], fn = logLik, object = object,
                     method = "L-BFGS-B", lower = c(rep(-Inf, object$l), 1e-9),
                     control = list(trace = !silent, fnscale = -1))

  if (!silent) cat("\nPreparing iprior output... ")
  theta <- mod.optim$par

  suppressWarnings(mod.iprior <- iprior(object, control = list(silent = TRUE,
                                                               maxit = 1,
                                                               theta = theta)))
  if (!silent) cat("DONE.\n")
  mod.iprior
}
