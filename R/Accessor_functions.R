################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2017  Haziq Jamil
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

#' Accessor functions for \code{ipriorMod} objects.
#'
#' @param object An \code{ipriorMod} object.
#' @param units Units for object size.
#' @param standard Standard for object size.
#' @param theta (Optional) Value of hyperparameters to evaluate the kernel
#'   matrix.
#' @param xstar (Optional) If not supplied, then a square, symmetric kernel
#'   matrix is returned using the data as input points. Otherwise, the kernel
#'   matrix is evaluated with respect to this set of data as well. It must be a
#'   list of vectors/matrices with similar dimensions to the original data.
#'
#' @name Accessors
NULL

#' @rdname Accessors
#' @export
get_intercept <- function(object) {
  check_and_get_ipriorKernel(object)
  object$intercept
}

#' @rdname Accessors
#' @export
get_y <- function(object) {
  check_and_get_ipriorKernel(object)
  if (is.iprobit(object)) {
    warning("Numerised categorical variables.", call. = FALSE)
  }
  res <- as.numeric(object$y + get_intercept(object))
  names(res) <- rownames(object$y)
  res
}

#' @rdname Accessors
#' @export
get_size <- function(object, units = "kB", standard = "SI") {
  check_and_get_ipriorKernel(object)
  print(object.size(object), units = units, standard = standard)
}

#' @rdname Accessors
#' @export
get_hyp <- function(object) {
  check_and_get_ipriorMod(object)
  object$param.full
}

#' @rdname Accessors
#' @export
get_lambda <- function(object) {
  tmp <- get_hyp(object)
  tmp[grep("lambda", names(tmp))]
}

#' @rdname Accessors
#' @export
get_psi <- function(object) {
  tmp <- get_hyp(object)
  tmp[grep("psi", names(tmp))]
}

#' @rdname Accessors
#' @export
get_se <- function(object) {
  check_and_get_ipriorMod(object)
  expand_theta(object$se, object$ipriorKernel$thetal$theta.drop, NA)
}

#' @rdname Accessors
#' @export
get_kernels <- function(object) {
  if (is.ipriorMod(object)) theta <- object$theta
  if (is.ipriorKernel2(object)) theta <- object$thetal$theta
  check_and_get_ipriorKernel(object)
  param.tab <- theta_to_param(theta, object)
  res <- param.tab$kernel
  names(res) <- object$xname
  res
}

#' @rdname Accessors
#' @export
get_kern_matrix <- function(object, theta = NULL, xstar = list(NULL)) {
  if (is.ipriorMod(object)) {
    # estl <- object$ipriorKernel$estl
    # til.cond <- (
    #   !isTRUE(estl$est.hurst) & !isTRUE(estl$est.lengt) & !isTRUE(estl$est.offs)
    # )
    res <- get_Hlam(object$ipriorKernel, object$theta, FALSE)
    return(res)
  } else if (is.ipriorKernel2(object)) {
    # estl <- object$estl
    # til.cond <- (
    #   !isTRUE(estl$est.hurst) & !isTRUE(estl$est.lengt) & !isTRUE(estl$est.offs)
    # )
    res <- get_Hlam(object, object$theta, FALSE)
    return(res)
  }
}

#' @rdname Accessors
#' @export
get_mse <- function(object) {
  check_and_get_ipriorMod(object)
  object$train.error
}

#' @rdname Accessors
#' @export
get_estl <- function(object) {
  check_and_get_ipriorKernel(object)
  unlist(object$estl)
}

#' @rdname Accessors
#' @export
get_method <- function(object) {
  check_and_get_ipriorMod(object)
  cat(object$est.method)
}

#' @rdname Accessors
#' @export
get_convergence <- function(object) {
  check_and_get_ipriorMod(object)
  cat(object$est.conv)
}

#' @rdname Accessors
#' @export
get_niter <- function(object) {
  check_and_get_ipriorMod(object)
  niter <- object$niter
  maxit <- object$control$maxit
  cat("Iterations:", paste0(niter, "/", maxit, "."))
}

#' @rdname Accessors
#' @export
get_time <- function(object) {
  check_and_get_ipriorMod(object)
  object$time
}
