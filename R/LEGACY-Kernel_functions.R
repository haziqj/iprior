#' ################################################################################
#' #
#' #   iprior: Linear Regression using I-priors
#' #   Copyright (C) 2017  Haziq Jamil
#' #
#' #   This program is free software: you can redistribute it and/or modify
#' #   it under the terms of the GNU General Public License as published by
#' #   the Free Software Foundation, either version 3 of the License, or
#' #   (at your option) any later version.
#' #
#' #   This program is distributed in the hope that it will be useful,
#' #   but WITHOUT ANY WARRANTY; without even the implied warranty of
#' #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' #   GNU General Public License for more details.
#' #
#' #   You should have received a copy of the GNU General Public License
#' #   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#' #
#' ################################################################################
#'
#' fnH1 <- kern_pearson
#'
#' fnH2 <- kern_linear
#'
#' fnH3 <- kern_fbm
#'
#' #' fn.H1 <- function(x, y = NULL) {
#' #'   # Pearson kernel.
#' #'   # Args: x,y vector of type "factor"
#' #' 	ytmp <- y
#' #' 	if (is.null(ytmp)) y <- x
#' #' 	if (any(is.numeric(x), is.numeric(y))) {
#' #' 	  warning("Non-factor type vector used with Pearson kernel.", call. = FALSE)
#' #' 	}
#' #' 	# Combine x and y, unfactorise them and work with numbers --------------------
#' #' 	x <- factor(x); y <- factor(y)
#' #' 	z <- unlist(list(x,y))  # simply doing c(x,y) messes with the factors
#' #' 	z <- as.numeric(z)
#' #' 	x <- z[1:length(x)]; y <- z[(length(x) + 1):length(z)]
#' #' 	if (any(is.na(match(y,x)))) {
#' #' 	  stop("The vector y contains elements not belonging to x.")
#' #' 	}
#' #' 	prop <- table(x)/length(x)
#' #'
#' #' 	unqy <- sort(unique(y))
#' #' 	unqx <- sort(unique(x))
#' #' 	tmpx <- lapply(unqx, function(k) which(x == k))
#' #' 	tmpy <- lapply(unqy, function(k) which(y == k))
#' #' 	tmp <- lapply(1:length(unqy),
#' #' 	              function(k) expand.grid(tmpy[[k]], tmpx[[unqy[k]]]))
#' #' 	# Side note: can avoid for loop below by combining the list tmp using
#' #' 	# do.call(rbind, tmp) or the faster data.table option
#' #' 	# as.matrix(data.table::rbindlist(tmp)) but tests found that this is actually
#' #' 	# slower.
#' #'
#' #' 	mat <- matrix(-1, nrow = length(y), ncol = length(x))
#' #' 	for (i in 1:length(unqy)) {
#' #' 		mat[as.matrix(tmp[[i]])] <- 1/prop[unqy[i]] - 1
#' #' 	}
#' #' 	class(mat) <- "Pearson"
#' #' 	mat
#' #' }
#' #'
#' #' # NOT USED YET
#' #' fn.H1a <- function(x, y = NULL) {
#' #'   # Efficient Pearson kernel
#' #'   # Args: x,y vector of type "factor"
#' #'   ytmp <- y
#' #'   if (is.null(ytmp)) y <- x
#' #'   if (any(is.numeric(x), is.numeric(y))) {
#' #'     warning("Non-factor type vector used with Pearson kernel.", call. = FALSE)
#' #'   }
#' #'   # Combine x and y, unfactorise them and work with numbers --------------------
#' #'   x <- factor(x); y <- factor(y)
#' #'   z <- unlist(list(x,y))  # simply doing c(x,y) messes with the factors
#' #'   z <- as.numeric(z)
#' #'   x <- z[1:length(x)]; y <- z[(length(x) + 1):length(z)]
#' #'   if (any(is.na(match(y,x)))) {
#' #'     stop("The vector y contains elements not belonging to x.")
#' #'   }
#' #'   prop <- table(x)/length(x)
#' #'
#' #'   unqy <- sort(unique(y))
#' #'   unqx <- sort(unique(x))
#' #'   tmp <- match(unqy, unqx)
#' #'
#' #'   mat <- matrix(-1, nrow = length(unqy), ncol = length(unqx))
#' #'   for (i in 1:length(tmp)) {
#' #'     mat[i, i] <- 1 / prop[unqy[i]] - 1
#' #'   }
#' #'   class(mat) <- "EfficientPearson"
#' #'   mat
#' #' }
#' #'
#' #'
#' #'
#' #' # The following three functions are able to take matrices instead of vector
#' #' # inputs. This is required for the one.lam = TRUE option. These functions are
#' #' # exported
#' #'
#' #' #' @rdname kernel
#' #' #' @export
#' #' fnH1 <- function(x, y = NULL){
#' #' 	res <- 0
#' #' 	if ((ncol(x) > 1) && !is.null(ncol(x))) {
#' #' 		if ((ncol(x) != ncol(y)) && !is.null(y)) {
#' #' 		  stop("New data is structurally unsimilar.")
#' #' 		}
#' #' 		for (i in 1:ncol(x)) res <- res + fn.H1(x = x[, i], y = y[, i])
#' #' 	}
#' #' 	else res <- fn.H1(x, y)
#' #' 	return(res)
#' #' }
#' #'
#' #' #' @rdname kernel
#' #' #' @export
#' #' fnH2 <- function(x, y = NULL) {
#' #'   # Centred Canonical kernel function. This is the kernel used, as opposed to
#' #'   # the uncentred one.
#' #'   rownames(x) <- colnames(x) <- rownames(y) <- colnames(y) <- NULL
#' #'   x <- scale(x, scale = FALSE)  # centre the variables
#' #'   if (is.null(y)) {
#' #'     tmp <- tcrossprod(x)
#' #'   } else {
#' #'     if (is.vector(y)) y <- matrix(y, ncol = ncol(x))
#' #'     else y <- as.matrix(y)
#' #'     y <- sweep(y, 2, attr(x ,"scaled:center"), "-")
#' #'     tmp <- tcrossprod(y, x)
#' #'   }
#' #'   class(tmp) <- "Canonical"
#' #'   tmp
#' #' }
#' #'
#' #' #' @rdname kernel
#' #' #' @export
#' #' fnH3 <- function(x, y = NULL, gamma = NULL) {
#' #'   if (is.null(gamma)) gamma <- 0.5
#' #'   if (is.vector(x)) x <- matrix(x, ncol = 1)
#' #'   x <- as.matrix(x)
#' #'   n <- nrow(x)
#' #'
#' #'   # fnNorm <- function(x) {
#' #'   #   sum(abs(x) ^ normtype) ^ (normtype * gamma / normtype)
#' #'   # }
#' #'
#' #'   A <- matrix(0, n, n)
#' #'   index.mat <- upper.tri(A)
#' #'   index <- which(index.mat, arr.ind = TRUE)
#' #'   xcrossprod <- tcrossprod(x)
#' #'   tmp1 <- diag(xcrossprod)[index[, 1]]
#' #'   tmp2 <- diag(xcrossprod)[index[, 2]]
#' #'   tmp3 <- xcrossprod[index]
#' #'   A[index.mat] <- tmp1 + tmp2 - 2 * tmp3
#' #'   A <- A + t(A)
#' #'   A <- abs(A) ^ gamma
#' #'   rvec <- apply(A, 1, sum)
#' #'   s <- sum(rvec)
#' #'
#' #'   if (is.null(y)) {
#' #'     rvec1 <- tcrossprod(rvec, rep(1, n))
#' #'     tmp <- (A - rvec1 / n - t(rvec1) / n + s / (n ^ 2)) / (-2)
#' #'   }
#' #'   else{
#' #'     if (is.vector(y)) y <- matrix(y, ncol = 1)
#' #'     else y <- as.matrix(y)
#' #'     m <- nrow(y)
#' #'
#' #'     rvec1 <- tcrossprod(rep(1, m), rvec)
#' #'     B <- matrix(0, m, n)
#' #'     indexy <- expand.grid(1:m, 1:n)
#' #'     ynorm <- apply(y, 1, function(x) sum(x ^ 2))
#' #'     xycrossprod <- tcrossprod(y, x)
#' #'     tmp1 <- ynorm[indexy[, 1]]
#' #'     tmp2 <- diag(xcrossprod)[indexy[, 2]]
#' #'     tmp3 <- as.numeric(xycrossprod)
#' #'     B[, ] <- tmp1 + tmp2 - 2 * tmp3
#' #'     neg.B <- B[B < 0]
#' #'     if (length(neg.B) > 0) {
#' #'       warning(c("These numbers are negative:", neg.B, ". Positive values taken for root."))
#' #'     }
#' #'     B <- abs(B) ^ gamma
#' #'     qvec <- apply(B, 1, sum)
#' #'     qvec1 <- tcrossprod(qvec, rep(1, n))
#' #'     tmp <- (B - qvec1 / n - rvec1 / n + s / (n ^ 2)) / (-2)
#' #'   }
#' #'   class(tmp) <- paste("FBM", gamma, sep = ",")
#' #'   tmp
#' #' }
