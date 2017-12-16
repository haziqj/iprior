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
#' #'DEPRECATED Load the kernel matrices of an I-prior model
#' #'
#' #'Prepare the kernel matrices according to a user available model options list.
#' #'This is then passed to the \code{\link{iprior}} function for model fitting.
#' #'Both formula and non-formula input are supported.
#' #'
#' #'When using non-formula to load the model, the explanatory variables can either
#' #'be vectors, matrices or data frames. These need to be entered one by one in
#' #'the function call, separated by commas. This is because each entry will have
#' #'one scale parameter attached to it. Like the \code{\link{iprior}} function,
#' #'grouping the scale parameters can only be done using non-formula input (see
#' #'examples).
#' #'
#' #'Sometimes, the model to be fitted can be quite complex and heavy for the EM
#' #'algorithm. Loading the data into an \code{ipriorKernel} object does the heavy
#' #'matrix matrix operations upfront, and passed on to the EM routine when
#' #'\code{\link{iprior}} is called.
#' #'
#' #'One advantage of having a saved \code{ipriorKernel} object is that we are able
#' #'to use any R optimiser and maximise the log-likelihood of the I-prior model in
#' #'conjunction with \code{\link{logLik}} or \code{\link{deviance}} functions.
#' #'
#' #'@param y Vector of response variables.
#' #'@param ... Only used for when fitting using non-formula, enter the variables
#' #'  (vectors or matrices) separated by commas. No other options applicable here.
#' #'@param model List of model options. Not used for \code{ipriorKernel} or
#' #'  \code{ipriorModel} objects. Available options are:
#' #'  \describe{\item{\code{kernel}}{Vector of character strings of either
#' #'  \code{"Canonical"}, \code{"FBM"}, or \code{"Pearson"}. Defaults to
#' #'  \code{"Canonical"} for continuous variables, and \code{"Pearson"} for factor
#' #'  type objects. To specify a Hurst coefficient, use \code{"FBM,<value>"};
#' #'  otherwise the default of 0.5 is used. Alternatively, see \code{Hurst} option
#' #'  below.} \item{\code{Hurst}}{Set the value of the Hurst coefficient for all
#' #'  \code{FBM} kernels used, rather than one by one. This is a value between 0
#' #'  and 1, and defaults to 0.5.}\item{\code{order}}{Character vector of length
#' #'  equal to the number of explanatory variables used, indicating specification
#' #'  of higher order scale parameters. The syntax is \code{"a^b"}, for parameter
#' #'  \code{a} raised to the power \code{b}. For regular order terms, then just
#' #'  input "a".} \item{\code{parsm}}{Logical, defaults to \code{TRUE}. Set to
#' #'  \code{FALSE} to assign one scale parameter for all kernel matrices.}
#' #'  \item{\code{one.lam}}{Logical, defaults to \code{FALSE}. Only relevant when
#' #'  using the formula call. Should all the variable share the same scale
#' #'  parameter?}\item{\code{rootkern}}{Logical, defaults to \code{FALSE}. Setting to \code{TRUE} is equivalent to Gaussian process regression.}}
#' #'
#' #'  These options are also available, but are only relevant when calling using
#' #'  non-formula: \describe{\item{\code{yname}}{Character vector to set the name
#' #'  of the response variable. It is set to the object name which contains the
#' #'  response variables by default.} \item{\code{xname}}{Character vector to set
#' #'  the name of the explanatory variables. This is also set to the object name
#' #'  by default.} \item{\code{interactions}}{Character vector to specify the
#' #'  interaction terms. When using formulas, this is specified automatically.
#' #'  Syntax is \code{"a:b"} to indicate variable \code{a} interacts with variable
#' #'  \code{b}.}}
#' #'@param formula The formula to fit when using formula interface.
#' #'@param data Data frame containing variables when using formula interface.
#' #'
#' #'@return A list of 11 items. Some of the more important ones are described
#' #'  below. \describe{ \item{\code{Y}}{The response variable.}
#' #'  \item{\code{x}}{The explanatory variables in list form. Each element of the
#' #'  list corresponds to each variable. If \code{one.lam = TRUE} was called, then
#' #'  you should see a single element in this list.} \item{\code{Hl}}{A list of
#' #'  the kernel matrices calculated from the explanatory variables according to
#' #'  the model options.} \item{\code{n, p, l, r, no.int, q}}{These are,
#' #'  respectively, the sample size, the number of explanatory variables, the
#' #'  number of unique scale parameters, the number of higher order terms, the
#' #'  number of interacting variables, and the number of kernel matrices.}}
#' #'
#' #'  The rest of the list are unimportant to the end-user, but they are passed to
#' #'  the EM routine via a call to \code{\link{iprior}}.
#' #'
#' #'@name kernL_old
#' #'@export
#' .kernL <- function(y, ..., model = list()) UseMethod(".kernL")
#'
#' #' @export
#' .kernL.default <- function(y, ..., model = list()) {
#'   x <- list(...)  # don't list if updating ipriorKernel
#'   if (length(x) == 1 && is.ipriorX(x[[1]])) x <- unlist(x, recursive = FALSE)
#'   if (testXForm(x)) x <- unlist(x, recursive = FALSE)
#'   x <- lapply(x, as.matrix)
#'   class(x) <- "ipriorX"
#'   n <- length(y)
#'   p <- length(x)
#'
#'   # Model options and checks ---------------------------------------------------
#'   mod <- list(kernel = "Canonical", Hurst = NULL, interactions = NULL,
#'               parsm = TRUE, one.lam = FALSE, yname = NULL, xname = NULL,
#'               order = as.character(1:p), intr.3plus = NULL, rootkern = FALSE,
#'               probit = FALSE, Nys.kern = FALSE, Nys.samp = NULL)
#'   mod_names <- names(mod)
#'   mod[(model_names <- names(model))] <- model
#'   if (length(noNms <- model_names[!model_names %in% mod_names])) {
#'     warning("Unknown names in model options: ", paste(noNms, collapse = ", "),
#'             call. = FALSE)
#'   }
#'
#'   # This part is for Nystrom (when called from ipriorNystrom function) ---------
#'   if (as.numeric(mod$Nys.kern) > 0) {
#'     if (is.null(mod$Nys.samp)) mod$Nys.samp <- sample(seq_len(n))
#'     y <- y[mod$Nys.samp]
#'     tmp <- lapply(x, rwa_1, smp = mod$Nys.samp)  # defined in .reorder_ipriorKernel()
#'     mostattributes(tmp) <- attributes(x)
#'     x <- tmp
#'     tmp <- lapply(x, rwa_1, smp = seq_len(mod$Nys.kern))
#'     mostattributes(tmp) <- attributes(x)
#'     x.Nys <- tmp
#'   }
#'
#'   # This part is for categorical response models -------------------------------
#'   y.levels <- NULL
#'   if (is.factor(y)) {
#'     mod$probit <- TRUE
#'     tmp <- .checkLevels(y)  # Utilities.R
#'     y <- tmp$y
#'     y.levels <- tmp$levels
#'   }
#'
#'   # What types of kernels? -----------------------------------------------------
#'   if (length(mod$kernel) < p && length(mod$kernel) > 1) {
#'     warning(paste0("Incomplete kernel specification (not of length ", p, ")"),
#'             call. = FALSE)
#'   }
#'   if (length(mod$kernel) > p && length(mod$kernel) > 1) {
#'     warning(paste0("Too many kernel options specification (not of length ", p, ")"),
#'             call. = FALSE)
#'   }
#'   Hurst <- kernel <- rep(NA, p)
#'   suppressWarnings(kernel[] <- splitKernel(mod$kernel))
#'   # The next two lines ensure that the Pearson kernel is used for factors
#'   whichPearson <- unlist(lapply(x, function(x) {is.factor(x) | is.character(x)}))
#'   kernel[whichPearson] <- "Pearson"
#'   check.kern <- match(kernel, c("FBM", "Canonical", "Pearson"))
#'   if (any(is.na(check.kern))) {
#'     stop("kernel should be one of \"Canonical\", \"Pearson\", or \"FBM\".",
#'          call. = FALSE)
#'   }
#'   suppressWarnings(Hurst[] <- splitHurst(mod$kernel))
#'   if (!is.null(mod$Hurst)) {
#'     # User has set a single Hurst coefficient for all FBM kernels
#'     if (any(!is.na(Hurst))) {
#'       warning("Overriding Hurst setting.", call. = FALSE)
#'     }
#'     suppressWarnings(Hurst[] <- mod$Hurst)
#'   }
#'   Hurst[is.na(Hurst)] <- 0.5
#'   mod$Hurst <- Hurst
#'   mod$kernel <- kernel
#'
#'   # Check for higher order terms -----------------------------------------------
#'   mod$order <- as.character(mod$order)
#'   suppressWarnings(hord.check <- all(sort(as.numeric(mod$order)) == 1:p))
#'   if (!hord.check) {
#'     hord.check1 <- length(mod$order) != p
#'     hord.check2 <- any(grepl("\\^", mod$order))
#'     if (hord.check1 | !hord.check2) {
#'       stop("Incorrect prescription of higher order terms.", call. = FALSE)
#'     }
#'   }
#'   r <- lenHOrd(mod$order)
#'
#'   # Set up interactions, p and q -----------------------------------------------
#'   names(mod)[3] <- "intr"  #rename to something simpler
#'   if (!is.null(mod$intr)) {
#'     # Interactions present
#'     if (!is.matrix(mod$intr)) {
#'       # Not fitted using formula
#'       intr.check1 <- is.character(mod$intr)
#'       intr.check2 <- all(grepl(":", mod$intr))
#'       if (!intr.check1 | !intr.check2) {
#'         stop("Incorrect prescription of interactions.", call. = FALSE)
#'       }
#'       ind.intr.3plus <- whichIntr3Plus(mod$intr)
#'       mod$intr.3plus <- mod$intr[ind.intr.3plus]
#'       mod$intr <- mod$intr[!ind.intr.3plus]
#'       mod$intr <- sapply(strsplit(mod$intr, ":"), as.numeric)
#'     }
#'     if (length(mod$intr) == 0) {
#'       mod[match("intr", names(mod))] <- list(NULL)
#'       no.int <- 0
#'     } else {
#'       no.int <- ncol(mod$intr)
#'     }
#'   } else {
#'     # No interactions
#'     no.int <- 0L
#'   }
#'   if (any(mod$intr > p | mod$intr < 1)) {
#'     stop("Prescribed interactions out of bounds.")
#'   }
#'   q <- p + no.int
#'   if (!mod$parsm) {
#'     l <- q
#'     mod$order <- as.character(1:l)
#'   } else {
#'     l <- p - r
#'   }
#'   # For clarity, the definitions of p, q, r, and l are
#'   # p = Number of x variables used
#'   # l = Number of unique lambdas (= q when parsm = FALSE)
#'   # r = Number of higher order terms
#'   # q = Length of expanded lambda = p + no.int
#'   # h = length(H.mat)
#'
#'   # More interactions ----------------------------------------------------------
#'   no.int.3plus <- 0
#'   if (!is.null(mod$intr.3plus) & length(mod$intr.3plus) > 0) {
#'     if (!is.matrix(mod$intr.3plus)) {
#'       mod$intr.3plus <- addZeroesIntr3Plus(mod$intr.3plus)
#'     }
#'     no.int.3plus <- ncol(mod$intr.3plus)
#'     q <- q + no.int.3plus
#'   }
#'
#'   # Set up names for x variables -----------------------------------------------
#'   if (is.null(mod$xname)) mod$xname <- names(x)
#'   else names(x) <- mod$xname[1:p]
#'   suppressWarnings(cond1 <- is.null(mod$xname))
#'   suppressWarnings(cond2 <- any(names(x) == ""))
#'   suppressWarnings(cond3 <- any(is.na(names(x))))
#'   cl <- match.call(); cl[[1L]] <- as.name("kernL")
#'   if (cond1 | cond2 | cond3) {
#'     m <- match(c("y", "model", "control"), names(cl), 0L)
#'     xnamefromcall <- as.character(cl[-m])[-1]
#'     mod$xname <- xnamefromcall
#'   }
#'   suppressWarnings(here <- which((names(x) != "") & !is.na(names(x))))
#'   mod$xname[here] <- names(x)[here]
#'   names(x) <- mod$xname[1:p]
#'
#'   # Set up name for y variable -------------------------------------------------
#'   ynamefromcall <- as.character(cl[2])
#'   check.yname <- is.null(mod$yname)
#'   if (check.yname) mod$yname <- ynamefromcall
#'
#'   # The following chunk checks whether the prescriped level 1 terms are --------
#'   # in order -------------------------------------------------------------------
#'   ord.ind <- whereOrd(mod$order)
#'   order.noh <- mod$order[ord.ind]
#'   hord.check3 <- any(order.noh != as.character(1:l))
#'   hord.check4 <- any(sort(as.numeric(order.noh)) != 1:l)
#'   if (hord.check3 | hord.check4) {
#'     warning("Incorrect prescription of level 1 terms - automatically fixed.", call. = FALSE)
#'   }
#'   mod$order[ord.ind] <- as.character(1:l)
#'   # Next check if higher order terms' kernels similar to the level 1 terms. Test
#'   # when parsm = TRUE.
#'   if (r > 0 && mod$parsm) {
#'     order.h <- mod$order[-ord.ind]
#'     index.h <- which(isHOrd(mod$order))
#'     for (i in 1:r) {
#'       j <- as.numeric(splitHOrd(order.h[i]))[1]
#'       if (mod$kernel[j] != mod$kernel[index.h[i]]) {
#'         warning(paste("Kernel for variable", mod$xname[index.h[i]], "not the same as ", mod$xname[j]), call. = FALSE)
#'       }
#'     }
#'   }
#'
#'   # Set up list of H matrices --------------------------------------------------
#'   # note: hMatList() is in Utitilities.R
#'   if (as.numeric(mod$Nys.kern) > 0) {
#'     Hl <- .hMatList(x = x, kernel = mod$kernel, intr = mod$intr, no.int = no.int,
#'                     gamma = mod$Hurst, intr.3plus = mod$intr.3plus,
#'                     rootkern = mod$rootkern, xstar = x.Nys)
#'   } else {
#'     Hl <- .hMatList(x = x, kernel = mod$kernel, intr = mod$intr, no.int = no.int,
#'                     gamma = mod$Hurst, intr.3plus = mod$intr.3plus,
#'                     rootkern = mod$rootkern)
#'   }
#'   h <- length(Hl)
#'   names(Hl) <- mod$xname[1:h]
#'   if (length(mod$xname) < h && !mod$one.lam && !is.null(mod$intr)) {
#'     for (i in 1:no.int) {
#'       mod$xname <- c(mod$xname, paste(mod$xname[mod$intr[1, i]],
#'                                       mod$xname[mod$intr[2, i]], sep = ":"))
#'     }
#'   }
#'   if (length(mod$xname) < h && !mod$one.lam && !is.null(mod$intr.3plus)) {
#'     for (i in 1:no.int.3plus) {
#'       mod$xname <- c(mod$xname, paste(mod$xname[mod$intr.3plus[, i]],
#'                                       collapse = ":"))
#'     }
#'   }
#'   names(Hl) <- mod$xname
#'
#'   # Set up names for lambda parameters -----------------------------------------
#'   mod$lamnamesx <- mod$xname[whereOrd(mod$order)]
#'
#'   # Block B update function ----------------------------------------------------
#'   intr <- mod$intr
#'   environment(indxFn) <- environment()
#'   H2l <- Hsql <- Pl <- Psql <- Sl <- ind <- ind1 <- ind2 <- NULL
#'   BlockB <- function(k, x = lambda) NULL
#'   if (r == 0L & no.int.3plus == 0L & as.numeric(mod$Nys.kern) == 0L) {
#'     # No need to do all the below Block B stuff if higher order terms involved.
#'     # Also no need if Nys.kern option called.
#'     if (q == 1L) {
#'       Pl <- Hl
#'       Psql <- list(fastSquare(Pl[[1]]))
#'       Sl <- list(matrix(0, nrow = n, ncol = n))
#'     } else {
#'       # Next, prepare the indices required for indxFn().
#'       z <- 1:h
#'       ind1 <- rep(z, times = (length(z) - 1):0)
#'       ind2 <- unlist(lapply(2:length(z), function(x) c(NA, z)[-(0:x)]))
#'       # Prepare the cross-product terms of squared kernel matrices. This is a
#'       # list of q_choose_2.
#'       for (j in 1:length(ind1)) {
#'         H2l.tmp <- Hl[[ind1[j]]] %*% Hl[[ind2[j]]]
#'         H2l[[j]] <- H2l.tmp + t(H2l.tmp)
#'       }
#'
#'       if (!is.null(intr) && mod$parsm) {
#'         # CASE: Parsimonious interactions only ---------------------------------
#'         for (k in z) {
#'           Hsql[[k]] <- fastSquare(Hl[[k]])
#'           if (k <= p) ind[[k]] <- indxFn(k)  # only create indices for non-intr
#'         }
#'         BlockB <- function(k, x = lambda) {
#'           # Calculate Psql instead of directly P %*% P because this way
#'           # is < O(n^3).
#'           indB <- ind[[k]]
#'           lambda.P <- c(1, x[indB$k.int.lam])
#'           Pl[[k]] <<- Reduce("+", mapply("*", Hl[c(k, indB$k.int)], lambda.P,
#'                                          SIMPLIFY = FALSE))
#'           Psql[[k]] <<- Reduce("+", mapply("*", Hsql[indB$Psq],
#'                                            c(1, x[indB$Psq.lam] ^ 2),
#'                                            SIMPLIFY = FALSE))
#'           if (!is.null(indB$P2.lam1)) {
#'             lambda.P2 <- c(rep(1, sum(indB$P2.lam1 == 0)), x[indB$P2.lam1])
#'             lambda.P2 <- lambda.P2 * x[indB$P2.lam2]
#'             Psql[[k]] <<- Psql[[k]] +
#'               Reduce("+", mapply("*", H2l[indB$P2], lambda.P2, SIMPLIFY = FALSE))
#'           }
#'           lambda.PRU <- c(rep(1, sum(indB$PRU.lam1 == 0)), x[indB$PRU.lam1])
#'           lambda.PRU <- lambda.PRU * x[indB$PRU.lam2]
#'           Sl[[k]] <<- Reduce("+", mapply("*", H2l[indB$PRU], lambda.PRU,
#'                                        SIMPLIFY = FALSE))
#'         }
#'       } else {
#'         # CASE: Multiple lambda with no interactions, or with non-parsimonious -
#'         # interactions ---------------------------------------------------------
#'         for (k in 1:q) {
#'           Pl[[k]] <- Hl[[k]]
#'           Psql[[k]] <- fastSquare(Pl[[k]])
#'         }
#'         BlockB <- function(k, x = lambda) {
#'           ind <- which(ind1 == k | ind2 == k)
#'           Sl[[k]] <<- Reduce("+", mapply("*", H2l[ind], x[-k], SIMPLIFY = FALSE))
#'         }
#'       }
#'     }
#'   }
#'
#'   BlockBstuff <- list(H2l = H2l, Hsql = Hsql, Pl = Pl, Psql = Psql, Sl = Sl,
#'                       ind1 = ind1, ind2 = ind2, ind = ind, BlockB = BlockB)
#'   kernelLoaded <- list(Y = y, x = x, Hl = Hl, n = n, p = p, l = l, r = r,
#'                        no.int = no.int, q = q, Nystrom = FALSE, y.levels = y.levels,
#'                        BlockBstuff = BlockBstuff, model = mod, call = cl,
#'                        no.int.3plus = no.int.3plus)
#'   class(kernelLoaded) <- "ipriorKernel_old"
#'   if (as.numeric(mod$Nys.kern) > 0L)
#'     class(kernelLoaded) <- c("ipriorKernel_Nystrom")
#'   kernelLoaded
#' }
#'
#' #' @rdname kernL_old
#' #' @export
#' .kernL.formula <- function(formula, data, model = list(), ...) {
#'   mf <- model.frame(formula = formula, data = data)
#'   tt <- terms(mf)
#'   Terms <- delete.response(tt)
#'   x <- model.frame(Terms, mf)
#'   y <- model.response(mf)
#'   yname <- names(attr(tt, "dataClasses"))[1]
#'   xname <- names(x)
#'   xnl <- length(xname)
#'
#'   # For interactions -----------------------------------------------------------
#'   interactions <- NULL
#'   tmpo <- attr(tt, "order")
#'   # if (any(tmpo > 2)) {
#'   #   stop("iprior does not currently work with higher order interactions.")
#'   # }
#'   tmpf <- attr(tt, "factors")
#'   tmpf2 <- as.matrix(tmpf[-1, tmpo == 2])  # this obtains 2nd order interactions
#'   int2 <- apply(tmpf2, 2, function(x) which(x == 1))
#'   if (any(tmpo == 2)) interactions <- int2
#'
#'   # > 2-way interactions -------------------------------------------------------
#'   intr.3plus <- NULL
#'   tmpf3 <- as.matrix(tmpf[-1, tmpo > 2])
#'   int3 <- apply(tmpf3, 2, whereInt)
#'   if (any(tmpo > 2)) intr.3plus <- int3
#'
#'   # Deal with one.lam option ---------------------------------------------------
#'   one.lam <- FALSE
#'   if (any(names(model) == "one.lam")) one.lam <- model$one.lam
#'   if (one.lam) {
#'     if (!is.null(interactions)) {
#'       stop("Cannot use option one.lam = TRUE with interactions.", call. = FALSE)
#'     }
#'     if (ncol(x) == 1) {
#'       message("Option one.lam = TRUE used with a single covariate anyway.")
#'     }
#'     attributes(x)$terms <- attributes(x)$names <- NULL
#'     if (xnl <= 3) {
#'       xname <- paste(xname, collapse = " + ")
#'     } else {
#'       xname <- paste(xname[1], "+ ... +", xname[xnl])
#'     }
#'     x <- as.data.frame(x)
#'   }
#'
#'   kernelLoaded <- .kernL(y = y, x, model = c(model,
#'                                             list(interactions = interactions,
#'                                                  intr.3plus = intr.3plus,
#'                                                  yname = yname, xname = xname)))
#'
#'   # Changing the call to simply kernL ------------------------------------------
#'   cl <- match.call()
#'   cl[[1L]] <- as.name("kernL")
#'   m <- match(c("formula", "data"), names(cl), 0L)
#'   cl <- cl[c(1L, m)]
#'   kernelLoaded$call <- cl
#'   names(kernelLoaded$call)[2] <- "formula"
#'   kernelLoaded$terms <- tt
#'   kernelLoaded$formula <- formula
#'   kernelLoaded
#' }
#'
#' #' @export
#' print.ipriorKernel_old <- function(x, ...) {
#'   cat("\n")
#'   # if (x$model$kernel == 'Canonical') CanOrFBM <- 'Canonical' else CanOrFBM <-
#'   # paste0('Fractional Brownian Motion with Hurst coef. ', x$gamfbm) kerneltypes <-
#'   # c(CanOrFBM, 'Pearson', paste(CanOrFBM, '& Pearson')) if (all(x$whichPearson))
#'   # cat(kerneltypes[2], 'RKHS loaded') else { if (!all(x$whichPearson) &&
#'   # !any(x$whichPearson)) cat(kerneltypes[1], 'RKHS loaded') else
#'   # cat(kerneltypes[3], 'RKHS loaded') } if (x$q == 1 | x$model$one.lam) cat(',
#'   # with a single scale parameter.\n') else cat(', with', x$q, 'scale
#'   # parameters.\n')
#'   if (isTRUE(x$model$probit)) cat("Categorical response variables\n")
#'   cat("Sample size = ", x$n, "\n")
#'   cat("Number of x variables, p = ", x$p, "\n")
#'   cat("Number of scale parameters, l = ", x$l, "\n")
#'   cat("Number of interactions = ", x$no.int + x$no.int.3plus, "\n")
#'   if (x$model$rootkern) cat("\nInfo on root H matrix:\n\n")
#'   else cat("\nInfo on H matrix:\n\n")
#'   str(x$Hl, ...)
#'   cat("\n")
#' }
