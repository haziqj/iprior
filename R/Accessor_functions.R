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

get_intercept <- function(object) {
  check_and_get_ipriorKernel(object)
  object$intercept
}

get_y <- function(object) {
  check_and_get_ipriorKernel(object)
  if (is.iprobit(object)) {
    warning("Numerised categorical variables.", call. = FALSE)
  }
  res <- as.numeric(object$y + get_intercept(object))
  names(res) <- rownames(object$y)
  res
}
