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

#' Extract the variance of the responses
#'
#' Extract the variance of the responses, \eqn{\mathbf V_y}, of a fitted I-prior
#' model.
#'
#' For an I-prior model, the variance is given by \deqn{\mathbf V_y = \psi\mathbf H_\lambda^2 + \psi^{-1}\mathbf I_n}.
#'
#' @param object An object of class \code{ipriorMod} or \code{ipriorKernel}.
#'
#' @examples
#' mod.fit <- iprior(stack.loss ~ ., stackloss)
#' Var.Y <- vary(mod.fit)
#' str(Var.Y)
#'
#' @export
vary <- function(object) {
  if (!is.ipriorMod(object)) {
    stop("Input objects of class ipriorMod only.", call. = FALSE)
  }

  Hlam.mat <- Hlam(object)
  Hlam.matsq <- fastSquare(Hlam.mat)
  Var.Y <- object$psi * Hlam.matsq + diag(1 / object$psi, nrow(Hlam.mat))
  Var.Y
}

# varyinv <- function(object, theta = NULL) {
#   UseMethod("varyi")
# }
#
# varyinv.ipriorKernel <- function(object, theta = NULL) {
#   lambda <- theta[-length(theta)]
#   psi <- theta[length(theta)]
#   tmp <- ipriorEM(object, maxit = 0, silent = TRUE, lambda.init = lambda,
#                   psi.init = psi, clean = TRUE, getVarY = TRUE)
#   tmp
# }
