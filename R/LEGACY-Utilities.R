#' ################################################################################
#' #
#' #   iprior: Linear Regression using I-priors
#' #   Copyright (C) 2016  Haziq Jamil
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
#' # KERNEL LOADER HELPER FUNCTIONS -----------------------------------------------
#'
#' indxFn <- index_fn_B
#'
#' addZeroesIntr3Plus <- tab_intr_3plus
#'
#' whereInt <- where_int
#'
#' whichIntr3Plus <- which_intr_3plus
#'
#' splitKernel <- function(kernel) {
#'   # Helper function to split the FBMs from the Hurst coefficients, if any
#'   paste(lapply(strsplit(kernel, ","), function(x) x[1]))
#' }
#'
#' splitHurst <- function(kernel) {
#'   # Helper function to split the FBMs from the Hurst coefficients, if any
#'   suppressWarnings(
#'     tmp <- as.numeric(paste(lapply(strsplit(kernel, ","), function(x) x[2])))
#'   )
#'   # tmp[is.na(tmp)] <- 0.5
#'   tmp
#' }
#'
#' # UTILITY FUNCTIONS ------------------------------------------------------------
#'
#' # Although this one should be replaced by get_Hl() in the future.
#' #' @export
#' .hMatList <- function(x, kernel, intr, no.int, gamma, intr.3plus, rootkern,
#'                       xstar = vector("list", p)) {
#'   # Helper function for creation of list of H matrices. Used in Kernel_loader.r
#'   # and predict.R
#'   p <- length(x)
#'
#'   # Check how many Hurst coefficients are provided -----------------------------
#'   # if (length(gamma) < sum(isFBM(kernel))) {
#'   #   warning("Number of Hurst coefficients supplied is less than the number of FBM kernels used.", call. = FALSE)
#'   # }
#'   # if (length(gamma) > p) {
#'   #   stop("Number of Hurst coefficients supplied is more than the number of FBM kernels used.", call. = FALSE)
#'   # }
#'
#'   suppressWarnings(
#'     H <- mapply(canPeaFBM, x = x, kernel = as.list(kernel), gamma = gamma,
#'                 y = xstar, rootkern = rootkern, SIMPLIFY = FALSE)
#'   )
#'   if (!is.null(intr)) {
#'     # Add in interactions, if any.
#'     for (j in 1:no.int) {
#'       H[[p + j]] <- H[[intr[1, j]]] * H[[intr[2, j]]]
#'       class(H[[p + j]]) <- paste(class(H[[intr[1,j]]]), class(H[[intr[2,j]]]),
#'                                  sep = " x ")
#'     }
#'   }
#'   if (!is.null(intr.3plus) & length(intr.3plus) > 0) {
#'     no.int.3plus <- ncol(intr.3plus)
#'     for (j in 1:no.int.3plus) {
#'       H[[p + j + no.int]] <- Reduce("*", H[intr.3plus[, j]])
#'       intr.3plus.class <- sapply(H[intr.3plus[, j]], class)
#'       class(H[[p + j + no.int]]) <- paste(intr.3plus.class, collapse = " x ")
#'     }
#'   }
#'   H
#' }
#'
#' #' @export
#' .lambdaExpand <- function(x = lambda, env = ipriorEM.env) {
#'   # Expands lambda from length l to correct size q = p + no.int, first by
#'   # expanding the higher order terms (if any), and then by adding the
#'   # interaction lambdas after that.
#'   lambda.tmp <- rep(NA, q)
#'   for (i in 1:q) {
#'     if (isTRUE(probit)) {
#'       lambda.tmp[i] <- x[as.numeric(order[i])]
#'     } else {
#'       if (isHOrd(order[i])) {
#'         j.and.pow <- splitHOrd(order[i])
#'         j <- j.and.pow[1]
#'         pow <- j.and.pow[2]
#'         lambda.tmp[i] <- x[as.numeric(j)] ^ as.numeric(pow)
#'       }
#'       else lambda.tmp[i] <- x[as.numeric(order[i])]
#'     }
#'   }
#'   if (parsm && no.int > 0) {
#'     for (j in 1:no.int) {
#'       add.lam <- lambda.tmp[intr[1, j]] * lambda.tmp[intr[2, j]]
#'       lambda.tmp[p + j] <- add.lam
#'     }
#'   }
#'   if (no.int.3plus > 0) {
#'     for (j in 1:no.int.3plus) {
#'       lambda.tmp[p + j + no.int] <- Reduce("*", lambda.tmp[intr.3plus[, j]])
#'     }
#'   }
#'   assign("lambda", lambda.tmp, envir = env)
#' }
#'
#' #' @export
#' .lambdaContract <- function(x = lambda, env = ipriorEM.env) {
#'   # The opposite of lambdaExpand(). Looks for model$order vector and extracts
#'   # only the l lambdas.
#'   assign("lambda", x[whereOrd(order)], envir = env)
#' }
#'
#'
#' is.ipriorKernel_Nystrom <- function(x) inherits(x, "ipriorKernel_Nystrom")
#'
#' isNystrom <- function(x) {
#'   if (!is.list(x$Nystrom)) return(x$Nystrom)
#'   else return(TRUE)
#' }
#'
#' isHOrd <- function(x) {
#'   # Tests whether x contains ^ indicating higher order term.
#'   grepl("\\^", x)
#' }
#'
#' whereOrd <- function(x) {
#'   # Index of non-higher order terms.
#'   grep("\\^", x, invert = TRUE)
#' }
#'
#' lenHOrd <- function(x) {
#'   # How many higher order terms have been specified?
#'   length(grep("\\^", x))
#' }
#'
#' splitHOrd <- function(x) {
#'   # Gets the level 1 index and the power it is raised to
#'   strsplit(x, "\\^")[[1]]
#' }
#'
#' isCan <- function(x) x == "Canonical"
#'
#' isPea <- function(x) x == "Pearson"
#'
#' isFBM <- function(x) grepl("FBM", x)
#'
#' is.ipriorX <- function(x) inherits(x, "ipriorX")
#'
#' testXForm <- function(x) {
#'   # Tests whether object x is a data frame fitted using formula interface.
#'   xform <- FALSE
#'   if (length(x) == 1) {
#'     if (is.data.frame(x[[1]])) {
#'       xform <- !is.null(attr(x[[1]], "terms"))
#'     }
#'   }
#'   xform
#' }
#'
#'
#' canPeaFBM <- function(x, kernel, gamma, y, rootkern = FALSE) {
#'   if (isCan(kernel)) res <- fnH2(x, y)
#'   if (isPea(kernel)) res <- fnH1(x, y)
#'   if (isFBM(kernel)) res <- fnH3(x, y, gamma)
#'   if (rootkern) {
#'     classres <- paste0("r", class(res))
#'     res <- fastSquareRoot2(res)
#'     class(res) <- classres
#'     res
#'   } else {
#'     res
#'   }
#' }
#'
#' #' @export
#' .reorder_ipriorKernel <- function(object, Nys.samp = NULL) {
#'   if (is.null(Nys.samp)) Nys.samp <- sample(seq_len(object$n))
#'   # y and X
#'   object$Y <- object$Y[Nys.samp]
#'   tmp <- lapply(object$x, rwa_1, smp = Nys.samp)
#'   mostattributes(tmp) <- attributes(object$x)
#'   object$x <- tmp
#'
#'   # Hl
#'   if (!is.ipriorKernel_Nystrom(object))
#'     object$Hl <- lapply(object$Hl, rwa_2, smp = Nys.samp)
#'
#'   # In BlockBstuff
#'   if (!is.null(object$BlockBstuff$H2l))
#'     object$BlockBstuff$H2l <- lapply(object$BlockBstuff$H2l, rwa_2, smp = Nys.samp)
#'   if (!is.null(object$BlockBstuff$Hsql))
#'     object$BlockBstuff$Hsql <- lapply(object$BlockBstuff$Hsql, rwa_2, smp = Nys.samp)
#'   if (!is.null(object$BlockBstuff$Pl))
#'     object$BlockBstuff$Pl <- lapply(object$BlockBstuff$Pl, rwa_2, smp = Nys.samp)
#'   if (!is.null(object$BlockBstuff$Psql))
#'     object$BlockBstuff$Psql <- lapply(object$BlockBstuff$Psql, rwa_2, smp = Nys.samp)
#'   if (!is.null(object$BlockBstuff$Sl))
#'     object$BlockBstuff$Sl <- lapply(object$BlockBstuff$Sl, rwa_2, smp = Nys.samp)
#'
#'   object
#' }
#'
#' rwa_1 <- function(z, smp) {
#'   # Reorder with attributes
#'   res <- z[smp, ]
#'   mostattributes(res) <- attributes(z)
#'   res
#' }
#'
#' rwa_2 <- function(z, smp) {
#'   # Reorder with attributes
#'   res <- z[smp, smp]
#'   mostattributes(res) <- attributes(z)
#'   res
#' }
#'
#' triangIndex <- function(k){
#'   # Function to list row and column index of upper triangular matrix including
#'   # diagonals.
#'   w <- 1:k
#'   cbind(
#'     row = rep(w, times = length(w):1 ) ,
#'     col = unlist(lapply(1:length(w), function(x) c(NA,w)[-(0:x)]))
#'   )
#' }
