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

.reorder_ipriorKernel <- function(object, Nys.samp) {
  # y and X
  object$Y <- object$Y[Nys.samp]
  tmp <- lapply(object$x, rwa_1, smp = Nys.samp)
  mostattributes(tmp) <- attributes(object$x)
  object$x <- tmp

  # Hl
  object$Hl <- lapply(object$Hl, rwa_2, smp = Nys.samp)

  # In BlockBstuff
  if (!is.null(object$BlockBstuff$H2l))
    object$BlockBstuff$H2l <- lapply(object$BlockBstuff$H2l, rwa_2, smp = Nys.samp)
  if (!is.null(object$BlockBstuff$Hsql))
    object$BlockBstuff$Hsql <- lapply(object$BlockBstuff$Hsql, rwa_2, smp = Nys.samp)
  if (!is.null(object$BlockBstuff$Pl))
    object$BlockBstuff$Pl <- lapply(object$BlockBstuff$Pl, rwa_2, smp = Nys.samp)
  if (!is.null(object$BlockBstuff$Psql))
    object$BlockBstuff$Psql <- lapply(object$BlockBstuff$Psql, rwa_2, smp = Nys.samp)
  if (!is.null(object$BlockBstuff$Sl))
    object$BlockBstuff$Sl <- lapply(object$BlockBstuff$Sl, rwa_2, smp = Nys.samp)

  object
}

rwa_1 <- function(z, smp) {
  # Reorder with attributes
  res <- z[smp, ]
  mostattributes(res) <- attributes(z)
  res
}

rwa_2 <- function(z, smp) {
  # Reorder with attributes
  res <- z[smp, smp]
  mostattributes(res) <- attributes(z)
  res
}
