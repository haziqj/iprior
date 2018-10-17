################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2018  Haziq Jamil
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

#' Reproducing kernels for the I-prior package
#'
#' The kernel functions used in this package are: \itemize{ \item{The
#' (canonical) linear kernel} \item{The fractional Brownian motion (fBm) kernel
#' with Hurst index \eqn{\gamma}} \item{The Pearson kernel} \item{The (scaled)
#' \eqn{d}-degree polynomial kernel with offset \eqn{c}} \item{The squared
#' exponential (SE) kernel with lengthscale \eqn{l}} }
#'
#' The Pearson kernel is used for nominal-type variables, and thus
#' \code{\link{factor}}-type variables are treated with the Pearson kernel
#' automatically when fitting I-prior models. The other kernels are for
#' continuous variables, and each emits different properties of functions.
#'
#' The linear kernel is used for "straight-line" functions. In addition, if
#' squared, cubic, or higher order terms are to be modelled, then the polynomial
#' kernel is suitable for this purpose. For smoothing models, the fBm kernel is
#' preferred, although the SE kernel may be used as well.
#'
#' @param x A vector, matrix or data frame.
#' @param y (Optional) vector, matrix or data frame. \code{x} and \code{y} must
#'   have identical column sizes.
#' @param gamma The Hurst coefficient for the fBm kernel.
#' @param l The lengthscale for the SE kernel.
#' @param c The offset for the polynomial kernel. This is a value greater than
#'   zero.
#' @param d The degree for the polynomial kernel. This is an integer value
#'   greater than oe equal to two.
#' @param lam.poly The scale parameter for the polynomial kernel.
#' @param centre Logical. Whether to centre the data (default) or not.
#'
#' @return A matrix whose \code{[i, j]} entries are given by \eqn{h(\code{x[i]},
#'   \code{y[j]})}, with \code{h} being the appropriate kernel function. The
#'   matrix has dimensions \code{m} by \code{n} according to the lengths of
#'   \code{y} and \code{x} respectively. When a single argument \code{x} is
#'   supplied, then \code{y} is taken to be equal to \code{x}, and a symmetric
#'   \code{n} by \code{n} matrix is returned.
#'
#'   The matrix has a \code{"kernel"} attribute indicating which type of kernel
#'   function was called.
#'
#' @examples
#' kern_linear(1:3)
#' kern_fbm(1:5, 1:3, gamma = 0.7)
#'
#' @references \url{http://phd.haziqj.ml/intro/}
#'
#' @name kernel
#' @aliases kernels
NULL

#' @rdname kernel
#' @export
kern_canonical <- function(x, y = NULL, centre = TRUE) {
  list2env(kern_check_xy(x, y, centre), environment())

  if (is.null(y)) res <- tcrossprod(x)
  else res <- tcrossprod(y, x)

  attributes(res)$kernel <- "linear"
  res
}

#' @rdname kernel
#' @export
kern_linear <- kern_canonical

#' @rdname kernel
#' @export
kern_pearson <- function(x, y = NULL) {
  # vectors only
  ytmp <- y
  if (is.null(ytmp)) y <- x
  if (any(!is.factor(x), !is.factor(y))) {
    warning("Non-factor type vector used with Pearson kernel.", call. = FALSE)
  }

  # Combine x and y, unfactorise them and work with numbers --------------------
  x <- factor(x); y <- factor(y)
  z <- unlist(list(x, y))  # simply doing c(x, y) messes with the factors
  z <- as.numeric(z)
  x <- z[seq_along(x)]; y <- z[-seq_along(x)]
  if (any(is.na(match(y, x)))) {
    stop("The vector y contains elements not belonging to x.")
  }
  prop <- table(x) / length(x)

  unqy <- sort(unique(y))
  unqx <- sort(unique(x))
  tmpx <- lapply(unqx, function(k) which(x == k))
  tmpy <- lapply(unqy, function(k) which(y == k))
  tmp <- lapply(seq_along(unqy),
                function(k) expand.grid(tmpy[[k]], tmpx[[unqy[k]]]))

  # Side note: can avoid for loop below by combining the list tmp using
  # do.call(rbind, tmp) or the faster data.table option
  # as.matrix(data.table::rbindlist(tmp)) but tests found that this is actually
  # slower.
  res <- matrix(-1, nrow = length(y), ncol = length(x))
  for (i in seq_along(unqy)) {
    res[as.matrix(tmp[[i]])] <- 1 / prop[unqy[i]] - 1
  }

  attributes(res)$kernel <- "pearson"
  res
}

#' @rdname kernel
#' @export
kern_fbm <- function(x, y = NULL, gamma = 0.5, centre = TRUE) {
  list2env(kern_check_xy(x, y, FALSE), environment())

  A <- matrix(0, n, n)
  index.mat <- upper.tri(A)
  index <- which(index.mat, arr.ind = TRUE)
  xcrossprod <- tcrossprod(x)
  tmp1 <- diag(xcrossprod)[index[, 1]]
  tmp2 <- diag(xcrossprod)[index[, 2]]
  tmp3 <- xcrossprod[index]
  A[index.mat] <- tmp1 + tmp2 - 2 * tmp3
  A <- A + t(A)
  A <- abs(A) ^ gamma
  rvec <- apply(A, 1, sum)
  s <- sum(rvec)

  if (is.null(y)) {
    if (isTRUE(centre)) {
      rvec1 <- tcrossprod(rvec, rep(1, n))
      res <- (A - rvec1 / n - t(rvec1) / n + s / (n ^ 2)) / -2
    } else {
      a <- matrix(diag(xcrossprod), n, n) ^ gamma
      res <- (A - a - t(a)) / -2
    }
  } else {
    rvec1 <- tcrossprod(rep(1, m), rvec)
    B <- matrix(0, m, n)
    indexy <- expand.grid(1:m, 1:n)
    ynorm <- apply(y, 1, function(x) sum(x ^ 2))
    xycrossprod <- tcrossprod(y, x)
    tmp1 <- ynorm[indexy[, 1]]
    tmp2 <- diag(xcrossprod)[indexy[, 2]]
    tmp3 <- as.numeric(xycrossprod)
    B[, ] <- tmp1 + tmp2 - 2 * tmp3
    # neg.B <- B[B < 0]
    # if (length(neg.B) > 0) {
    #   warning(c("These numbers are negative:", neg.B,
    #             ". Positive values taken for root."))
    # }
    B <- abs(B) ^ gamma
    if (isTRUE(centre)) {
      qvec <- apply(B, 1, sum)
      qvec1 <- tcrossprod(qvec, rep(1, n))
      res <- (B - qvec1 / n - rvec1 / n + s / (n ^ 2)) / (-2)
    } else {
      bx <- matrix(diag(xcrossprod), nrow = m, ncol = n, byrow = TRUE) ^ gamma
      by <- matrix(ynorm, nrow = m, ncol = n) ^ gamma
      res <- (B - bx - by) / -2
    }
  }

  attributes(res)$kernel <- paste0("fbm,", gamma)
  res
}

#' @rdname kernel
#' @export
kern_se <- function(x, y = NULL, l = 1, centre = TRUE) {
  list2env(kern_check_xy(x, y, centre), environment())
  xcrossprod <- tcrossprod(x)

  if (is.null(y)) {
    A <- matrix(0, n, n)
    index.mat <- upper.tri(A)
    index <- which(index.mat, arr.ind = TRUE)
    tmp1 <- diag(xcrossprod)[index[, 1]]
    tmp2 <- diag(xcrossprod)[index[, 2]]
    tmp3 <- xcrossprod[index]
    A[index.mat] <- tmp1 + tmp2 - 2 * tmp3
    A <- A + t(A)
    xmxp.norm <- A
  } else {
    B <- matrix(0, m, n)
    indexy <- expand.grid(1:m, 1:n)
    ynorm <- apply(y, 1, function(x) sum(x ^ 2))
    xycrossprod <- tcrossprod(y, x)
    tmp1 <- ynorm[indexy[, 1]]
    tmp2 <- diag(xcrossprod)[indexy[, 2]]
    tmp3 <- as.numeric(xycrossprod)
    B[, ] <- tmp1 + tmp2 - 2 * tmp3
    xmxp.norm <- B
  }

  res <- exp(-xmxp.norm / (2 * l ^ 2))
  attributes(res)$kernel <- paste0("se,", l)
  if (isTRUE(centre)) return(kern_centre(res))
  else return(res)
}

#' @rdname kernel
#' @export
kern_poly <- function(x, y = NULL, c = 0, d = 2, lam.poly = 1, centre = TRUE) {
  if (!(as.integer(d) == d)) {
    stop("Non-integer value for polynomial degree d.", call. = FALSE)
  }
  if (d <= 1) {
    stop("Polynomial degree must be greater than 1.", call. = FALSE)
  }
  if (c < 0) {
    stop("Polynomial offset must be greater than 0.", call. = FALSE)
  }

  x.ip <- kern_canonical(x, y, centre = centre)
  res <- (lam.poly * x.ip + c) ^ d
  attributes(res)$kernel <- paste0("poly", d, ",", c)
  res
}

kern_check_xy <- function(x, y, centre.xy) {
  # Helper function to determine whether inputs for kernel functions are
  # structurally similar (same number of columns), and centres the data
  # appropriately.
  #
  # Args: Data x and y, and logical centre.xy.
  #
  # Returns: A list of processed data x and y, and the number of data points n
  # and m for the data respectively.
  if (is.vector(x)) x <- matrix(x)
  else x <- as.matrix(x)
  n <- nrow(x)

  m <- NULL
  if (!is.null(y)) {
    if (is.vector(y)) y <- matrix(y)
    else y <- as.matrix(y)
    if (ncol(y) != ncol(x)) stop("New data y is structurally unsimilar to x.")
    m <- nrow(y)
  }

  if (isTRUE(centre.xy)) {
    x <- scale(x, center = TRUE, scale = FALSE)
    if (!is.null(y)) {
      x.centre <- attr(x ,"scaled:center")
      y <- scale(y, center = x.centre, scale = FALSE)
    }
  }

  list(x = x, y = y, n = n, m = m)
}

kern_centre <- function(mat) {
  # Helper function to centre the kernels.
  #
  # Args: An uncentred kernel matrix.
  #
  # Returns: A centred kernel matrix.
  rvec <- apply(mat, 1, sum)
  cvec <- apply(mat, 2, sum)
  svec <- sum(mat)
  n <- ncol(mat)
  mat - rvec / n - rep(cvec, each = nrow(mat)) / n + svec / (n ^ 2)
}
