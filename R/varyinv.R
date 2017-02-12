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
#' Extract the variance of the responses of an I-prior model.
#'
#' For \code{ipriorKernel} objects, random values for \code{theta} are used.
#'
#' @param object An object of class \code{ipriorMod} or \code{ipriorKernel}.
#' @param theta (optional) Evaluates the log-likelihood at \code{theta} which is
#'   of the form \code{theta = c(lambda, psi)}.
#'
#' @export
varyinv <- function(object, theta = NULL) {
  UseMethod("varyinv")
}

#' @name varyinv
#' @export
varyinv.ipriorMod <- function(object, theta = NULL) {
  tmp <- with(object, {
    if (!is.null(theta)) {
      lambda <- theta[-length(theta)]
      psi <- theta[length(theta)]
    }
    ipriorEM(ipriorKernel, maxit = 0, silent = TRUE, lambda.init = lambda,
             psi.init = psi, clean = TRUE, getVarY = TRUE)
  })
  tmp
}

#' @name varyinv
#' @export
varyinv.ipriorKernel <- function(object, theta = NULL) {
  lambda <- theta[-length(theta)]
  psi <- theta[length(theta)]
  tmp <- ipriorEM(object, maxit = 0, silent = TRUE, lambda.init = lambda,
                  psi.init = psi, clean = TRUE, getVarY = TRUE)
  tmp
}











