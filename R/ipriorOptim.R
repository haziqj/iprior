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
#' First, the EM algorithm is performed for a maximum of three iterations. This
#' can be changed by passing a different \code{maxit} value to the list of
#' control options. The parameters are then passed to \code{optim()} and the
#' negative log-likelihood is minimised. The method used for optim is
#' \code{"L-BFGS-B"}, as the \code{psi} parameter of the I-prior model needs to
#' be constrained to be greater than zero.
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
ipriorOptim <- function(object, control = list(maxit = 3, report = 1)) {
  if (!is.ipriorKernel(object)) {
    stop("Input objects of class ipriorKernel only.", call. = FALSE)
  }

  if (is.null(control$maxit)) control$maxit <- 3
  if (is.null(control$report)) control$report <- 1
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
  mod.iprior$optim.converged <- TRUE
  mod.iprior
}
