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

#' Extract the scaled kernel matrix
#'
#' Extract the scaled kernel matrix \eqn{\mathbf H_\lambda} of an
#' \code{ipriorMod} or \code{ipriorKernel} object.
#'
#' The maximum likelihood values for \code{lambda} are used by default for
#' fitted I-prior models. For \code{ipriorKernel} objects, random values for
#' \code{lambda} are used.
#'
#' @param object An object of class \code{ipriorMod} or \code{ipriorKernel}.
#' @param lambda (optional) Values of the scale parameters.
#'
#' @examples
#' # Extracting from an ipriorMod object
#' mod.fit <- iprior(stack.loss ~ ., stackloss)
#' H.mat1 <- Hlam(mod.fit)
#' str(H.mat1)
#'
#' # Extracting from an ipriorKernel object at a specified lambda value
#' mod <- mod.fit$ipriorKernel
#' H.mat2 <- Hlam(mod, mod.fit$lambda)
#'
#' # They are both equal
#' all.equal(H.mat1, H.mat2)
#'
#' @export
Hlam <- function(object, lambda = NULL) {
  UseMethod("Hlam")
}

#' @name Hlam
#' @export
Hlam.ipriorMod <- function(object, lambda = NULL) {
  tmp <- with(
    object,
    ipriorEM(ipriorKernel, maxit = 0, silent = TRUE, lambda.init = lambda,
             psi.init = psi, clean = TRUE, getHlam = TRUE)
  )
  class(tmp) <- NULL
  tmp
}

#' @name Hlam
#' @export
Hlam.ipriorKernel <- function(object, lambda = NULL) {
  tmp <- ipriorEM(object, maxit = 0, silent = TRUE, lambda.init = lambda,
                  psi.init = NULL, clean = TRUE, getHlam = TRUE)
  class(tmp) <- NULL
  tmp
}
