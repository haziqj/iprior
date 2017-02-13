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

varyinv <- function(object, theta = NULL) {
  # Calculate the inverse variance of Y from ipriorMod and ipriorKernel objects.
  # Used in fisher() function.
  UseMethod("varyinv")
}

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

varyinv.ipriorKernel <- function(object, theta = NULL) {
  lambda <- theta[-length(theta)]
  psi <- theta[length(theta)]
  tmp <- ipriorEM(object, maxit = 0, silent = TRUE, lambda.init = lambda,
                  psi.init = psi, clean = TRUE, getVarY = TRUE)
  tmp
}
