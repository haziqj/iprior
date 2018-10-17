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

#' Load the kernel matrices for I-prior models
#'
#' @param y Vector of response variables
#' @param ... Only used when fitting using non-formula, enter the variables
#'   (vectors or matrices) separated by commas.
#' @param formula The formula to fit when using formula interface.
#' @param data Data frame containing variables when using formula interface.
#' @param kernel Character vector indicating the type of kernel for the
#'   variables. Available choices are: \itemize{ \item{\code{"linear"} -
#'   (default) for the linear kernel} \item{\code{"canonical"} - alternative
#'   name for \code{"linear"}} \item{\code{"fbm"}, \code{"fbm,0.5"} - for the
#'   fBm kernel with Hurst coefficient 0.5 (default)} \item{\code{"se"},
#'   \code{"se,1"} - for the SE kernel with lengthscale 1 (default)}
#'   \item{\code{"poly"}, \code{"poly2"}, \code{"poly2,0"} - for the polynomial
#'   kernel of degree 2 with offset 0 (default)} \item{\code{"pearson" - for the
#'   Pearson kernel}}} The \code{kernel} argument can also be a vector of length
#'   equal to the number of variables, therefore it is possible to specify
#'   different kernels for each variables. Note that factor type variables are
#'   assigned the Pearson kernel by default, and that non-factor types can be
#'   forced to use the Pearson kernel (not recommended).
#' @param interactions Character vector to specify the interaction terms. When
#'   using formulas, this is specified automatically, so is not required. Syntax
#'   is \code{"a:b"} to indicate variable \code{a} interacts with variable
#'   \code{b}.
#' @param est.lambda Logical. Estimate the scale parameters? Defaults to
#'   \code{TRUE}.
#' @param est.hurst Logical. Estimate the Hurst coefficients for fBm kernels?
#'   Defaults to \code{FALSE}.
#' @param est.lengthscale Logical. Estimate the lengthscales for SE kernels?
#'   Defaults to \code{FALSE}.
#' @param est.offset Logical. Estimate the offsets for polynomial kernels?
#'   Defaults to \code{FALSE}.
#' @param est.psi Logical. Estimate the error precision? Defaults to
#'   \code{TRUE}.
#' @param fixed.hyp Logical. If \code{TRUE}, then no hyperparameters are
#'   estimated, i.e. all of the above \code{est.x} are set to \code{FALSE}, and
#'   vice versa. If \code{NULL} (default) then all of the \code{est.x} defaults
#'   are respected.
#' @param lambda Initial/Default scale parameters. Relevant especially if
#'   \code{est.lambda = FALSE}.
#' @param psi Initial/Default value for error precision. Relevant especially if
#'   \code{est.psi = FALSE}.
#' @param nystrom Either logical or an integer indicating the number of Nystrom
#'   samples to take. Defaults to \code{FALSE}. If \code{TRUE}, then
#'   approximately 10\% of the sample size is used for the Nystrom
#'   approximation.
#' @param nys.seed The random seed for the Nystrom sampling. Defaults to
#'   \code{NULL}, which means the random seed is not fixed.
#' @param model DEPRECATED.
#' @param one.lam Logical. When using formula input, this is a convenient way of
#'   letting the function know to treat all variables as a single variable (i.e.
#'   shared scale parameter). Defaults to \code{FALSE}.
#' @param train.samp (Optional) A vector indicating which of the data points
#'   should be used for training, and the remaining used for testing.
#' @param test.samp (Optional) Similar to \code{train.samp}, but on test samples
#'   instead.
#'
#' @return An \code{ipriorKernel} object which contains the relevant material to
#'   be passed to the \code{iprior} function for model fitting.
#'
#' @seealso \link[=iprior]{iprior}
#'
#' @examples
#'
#' str(ToothGrowth)
#' (mod <- kernL(y = ToothGrowth$len,
#'                supp = ToothGrowth$supp,
#'                dose = ToothGrowth$dose,
#'                interactions = "1:2"))
#' kernL(len ~ supp * dose, data = ToothGrowth)  # equivalent formula call
#'
#' # Choosing different kernels
#' str(stackloss)
#' kernL(stack.loss ~ ., stackloss, kernel = "fbm")  # all fBm kernels
#' kernL(stack.loss ~ ., stackloss, kernel = "FBm")  # cApS dOn't MatTeR
#' kernL(stack.loss ~ ., stackloss,
#'        kernel = c("linear", "se", "poly3"))  # different kernels
#'
#' # Sometimes the print output is too long, can use str() options here
#' print(mod, strict.width = "cut", width = 30)
#'
#' @export
kernL <- function(
  y, ..., kernel = "linear", interactions = NULL, est.lambda = TRUE,
  est.hurst = FALSE, est.lengthscale = FALSE, est.offset = FALSE,
  est.psi = TRUE, fixed.hyp = NULL, lambda = 1, psi = 1, nystrom = FALSE,
  nys.seed = NULL, model = list(), train.samp, test.samp
) UseMethod("kernL")

#' @export
kernL.default <- function(y, ..., kernel = "linear", interactions = NULL,
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                           psi = 1, nystrom = FALSE, nys.seed = NULL,
                           model = list(), train.samp, test.samp) {
  # Checks ---------------------------------------------------------------------
  if (is.list(model) & length(model) > 0) {
    stop("\'model\' option is deprecated. Use the arguments directly instead. See \'?kernL\' for details.", call. = FALSE)
  }
  kernel <- tolower(kernel)
  Xl <- list(...)
  # It is common to make the mistake and type kernels instead of kernel. This
  # corrects it.
  Xl.kernel.mistake <- match("kernels", names(Xl))
  if ("kernels" %in% names(Xl)) {
    kernel <- Xl[[Xl.kernel.mistake]]
    Xl[[Xl.kernel.mistake]] <- NULL
  }
  # Check formula
  Xl.formula <- match("Xl.formula", names(Xl))
  formula.method <- FALSE
  if ("Xl.formula" %in% names(Xl)) {
    Xl <- Xl[[Xl.formula]]
    formula.method <- TRUE
  }
  # Get names
  xname <- names(Xl)
  yname <- attr(y, "yname")
  # Check for probit model
  if (is.factor(y)) {
    tmp <- get_y_and_levels(y)
    y <- tmp$y
    y.levels <- tmp$levels
  } else {
    y.levels <- NULL
  }
  one.lam <- attr(Xl, "one.lam")

  # Were training samples specified? -------------------------------------------
  train.check <- FALSE
  if (!missing(train.samp) | !missing(test.samp)) {
    if (!missing(train.samp) & !missing(test.samp)) {
      stop("Use either train.samp or test.samp only.", call. = FALSE)
    }
    if (missing(train.samp)) train.samp <- seq_along(y)[-test.samp]
    if (all(train.samp %in% seq_along(y))) {
      train.check <- TRUE
      test.samp <- seq_along(y)[-train.samp]
      y.test <- y[test.samp]
      Xl.test <- lapply(Xl, function(x) {
        if (is.matrix(x) | is.data.frame(x)) return(x[test.samp, , drop = FALSE])
        else return(x[test.samp])
      })
      y <- y[train.samp]
      Xl <- lapply(Xl, function(x) {
        if (is.matrix(x) | is.data.frame(x)) return(x[train.samp, , drop = FALSE])
        else return(x[train.samp])
      })
    } else {
      warning("Training samples incorrectly specified.", call. = FALSE)
    }
    if (!is.null(one.lam)) attr(Xl, "one.lam") <- attr(Xl.test, "one.lam") <-
        one.lam
  }

  # Get intercept --------------------------------------------------------------
  y <- scale(y, scale = FALSE)  # centre variables
  intercept <- attr(y, "scaled:center")

  # Meta -----------------------------------------------------------------------
  n <- length(y)
  p <- length(Xl)
  if (is.null(xname)) {
    xname <- paste0("X", seq_len(p))
  } else {
    where.blanks <- grepl("^$", xname)
    xname[where.blanks] <- paste0("X", which(where.blanks))
  }
  if (is.null(yname)) yname <- "y"

  # For Nystrom method: Reorder data and create Xl.Nys -------------------------
  if (as.numeric(nystrom) > 0 & as.numeric(nystrom) != n) {
    if (as.numeric(nystrom) == 1) nystrom <- floor(0.1 * n)
    if (!is.null(nys.seed)) set.seed(nys.seed)
    nys.samp <- sample(seq_along(y))
    y.tmp <- y[nys.samp]
    mostattributes(y.tmp) <- attributes(y)
    y <- y.tmp
    tmp <- lapply(Xl, reorder_x, smp = nys.samp)
    mostattributes(tmp) <- attributes(Xl)
    Xl <- tmp
    Xl.nys <- lapply(Xl, reorder_x, smp = seq_len(nystrom))
    mostattributes(Xl.nys) <- attributes(Xl)
    nys.check <- TRUE
  } else {
    nys.check <- FALSE
  }

  # What types of kernels? -----------------------------------------------------
  if (length(kernel) < p && length(kernel) > 1) {
    warning(paste0("Incomplete kernel specification (not of length ", p, ")"),
            call. = FALSE)
  }
  if (length(kernel) > p && length(kernel) > 1) {
    warning(paste0("Too many kernel options specification (not of length ", p,
                   ")"),
            call. = FALSE)
  }
  kernels <- rep(NA, p)
  suppressWarnings(kernels[1:p] <- kernel)
  # The next two lines ensure that the Pearson kernel is used for factors
  which.pearson <- unlist(lapply(Xl, function(x) {is.factor(x) |
      is.character(x)}))
  kernels <- correct_pearson_kernel(kernels, which.pearson)

  # Interactions ---------------------------------------------------------------
  intr <- intr.3plus <- NULL
  no.int.3plus <- no.int <- 0
  if (!is.null(interactions)) {
    if (isTRUE(formula.method)) {
      intr <- interactions$intr
      intr.3plus <- interactions$intr.3plus
    } else {
      ind.intr.3plus <- which_intr_3plus(interactions)
      intr.3plus <- interactions[ind.intr.3plus]
      intr <- interactions[!ind.intr.3plus]
      intr <- sapply(strsplit(intr, ":"), as.numeric)
    }
    if (length(intr > 0)) no.int <- ncol(intr)
  }
  if (length(intr.3plus) == 0) intr.3plus <- NULL
  if (!is.null(intr.3plus)) {
    if (!isTRUE(formula.method)) intr.3plus <- tab_intr_3plus(intr.3plus)
    no.int.3plus <- ncol(intr.3plus)
  }

  if (isTRUE(nys.check)) {
    Hl <- get_Hl(Xl, Xl.nys, kernels, lambda)
  } else {
    Hl <- get_Hl(Xl, list(NULL), kernels, lambda)
  }
  kernels <- get_kernels_from_Hl(Hl)

  if (!is.null(fixed.hyp)) {
    if (isTRUE(fixed.hyp)) {
      est.lambda <- est.hurst <- est.lengthscale <- est.offset <- est.psi <- FALSE
    }
    if (!isTRUE(fixed.hyp)) {
      est.lambda <- est.hurst <- est.lengthscale <- est.offset <- est.psi <- TRUE
    }
  }
  estl <- list(est.lambda = est.lambda, est.hurst = est.hurst,
               est.lengthscale = est.lengthscale, est.offset = est.offset,
               est.psi = est.psi)

  names(lambda) <- names(psi) <- NULL  # need to clean the names otherwise weird
                                       # things happen
  param <- kernel_to_param(kernels, lambda)
  poly.deg <- param$deg
  thetal <- param_to_theta(param, estl, log(psi))
  thetal$n.theta <- length(thetal$theta)

  nystroml <- NULL
  if (isTRUE(as.logical(nystrom))) {
    nystroml <- list(nys.samp = nys.samp, nys.seed = nys.seed,
                     nys.size = as.numeric(nystrom))
  }

  BlockBStuff <- NULL
  # Only calculate BlockBStuff when using the closed form EM algorithm.
  # |                 | EM.closed == TRUE |
  # |-----------------|------------------:|
  # | poly.deg        |              NULL |
  # | est.lambda      |              TRUE |
  # | est.hurst       |             FALSE |
  # | est.lengthscale |             FALSE |
  # | est.offset      |             FALSE |
  # | est.psi         |              TRUE |
  # | intr.3plus      |              NULL |
  BlockB.cond <- (
    all(is.na(poly.deg)) & !isTRUE(est.hurst) & !isTRUE(est.lengthscale) &
      !isTRUE(est.offset) & (isTRUE(est.lambda) | isTRUE(est.psi)) &
      !isTRUE(nys.check) & is.null(intr.3plus)
  )
  if (isTRUE(BlockB.cond)) {
    BlockBStuff <- BlockB_fn(Hl, intr, n, p)
  }

  res <- list(
    # Data
    y = y, Xl = Xl, Hl = Hl, intercept = intercept,
    # Model
    kernels = kernels, which.pearson = which.pearson,
    poly.deg = poly.deg, thetal = thetal, estl = estl,
    intr = intr, intr.3plus = intr.3plus, nystroml = nystroml,
    BlockBStuff = BlockBStuff,
    # Meta
    n = n, p = p, no.int = no.int, no.int.3plus = no.int.3plus,
    xname = xname, yname = yname, formula = NULL, terms = NULL,
    y.levels = y.levels
  )
  if (isTRUE(train.check)) {
    res$y.test <- y.test
    res$Xl.test <- Xl.test
  }
  res$call <- fix_call_default(match.call(), "kernL")  # fix function call

  class(res) <- "ipriorKernel"
  res
}

#' @rdname kernL
#' @export
kernL.formula <- function(formula, data, kernel = "linear", one.lam = FALSE,
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, fixed.hyp = NULL, lambda = 1,
                           psi = 1, nystrom = FALSE, nys.seed = NULL,
                           model = list(), train.samp, test.samp, ...) {
  list2env(formula_to_xy(formula = formula, data = data, one.lam = one.lam),
           envir = environment())
  res <- kernL.default(y = y, Xl.formula = Xl, interactions = interactions,
                       kernel = kernel, est.lambda = est.lambda,
                       est.hurst = est.hurst,
                       est.lengthscale = est.lengthscale,
                       est.offset = est.offset, est.psi = est.psi,
                       fixed.hyp = fixed.hyp, lambda = lambda, psi = psi,
                       nystrom = nystrom, nys.seed = nys.seed, model = model,
                       train.samp = train.samp, test.samp = test.samp, ...)
  res$yname <- yname
  res$formula <- formula
  res$terms <- tt
  res$call <- fix_call_formula(match.call(), "kernL")  # fix function call

  res
}

#' @export
print.ipriorKernel <- function(x, units = "auto", standard = "SI", ...) {
  tmp <- expand_Hl_and_lambda(x$Hl, seq_along(x$Hl), x$intr, x$intr.3plus)

  # if (isTRUE(x$probit)) {
  #   cat("Categorical response variables\n")
  # } else if (is.ipriorKernel_nys(x)) {
  #   cat("Nystrom kernel approximation ()\n")
  # }

  n <- x$n
  if (is.ipriorKernel_cv(x)) {
    n.test <- length(x$y.test)
    n <- paste0(n + n.test, " (", n, " train + ", n.test, " test)")
  }
  cat("Sample size:", n, "\n")
  cat("No. of covariates:", length(x$Xl), "\n")
  # cat("No. of interactions:", x$no.int + x$no.int.3plus, "\n")
  cat("Object size: ")
  print(object.size(x), units = units, standard = standard)

  cat("\n")
  cat("Kernel matrices:\n")
  for (i in seq_along(tmp$Hl)) {
    cat("", i, print_kern(tmp$Hl[[i]], ...), "\n")
  }
  cat("\n")
  cat("Hyperparameters to estimate:\n")
  if (x$thetal$n.theta > 0)
    cat(paste(names(x$thetal$theta), collapse = ", "))
  else
    cat("none")
  cat("\n")

  cat("\n")
  methods <- c("direct", "em", "canonical", "mixed", "fixed")
  poss.method <- NULL
  for (i in seq_along(methods)) {
    suppressWarnings(tmp <-  iprior_method_checker(x, methods[i]))
    poss.method <- c(poss.method, names(which(tmp)))
  }
  poss.method <- gsub("em.closed", "em", poss.method)
  poss.method <- gsub("em.reg", "em", poss.method)
  poss.method <- gsub("nystrom", "direct", poss.method)
  if (is.nystrom(x)) {
    poss.method <- paste(poss.method, "(Nystrom)")
  }
  if (is.categorical(x)) {
    poss.method <- c(poss.method, "iprobit (recommended)")
  }
  poss.method <- paste0(unique(poss.method), collapse = ", ")
  cat("Estimation methods available:\n")
  cat(poss.method)
  cat("\n")
}

print_kern <- function(x, ...) {
  # Helper function to prettify the print output of kernel matrices.
  #
  # Args: x is a kernel matrix obtained from one of the kern_x() functions.
  # Additional ... are passed to str().
  #
  # Returns: Prettified kernel matrix str() print ouput.
  kern.type <- attr(x, "kernel")
  res <- capture.output(str(x, ...))[1]
  res <- gsub(" num", kern.type, res)
  res
}
