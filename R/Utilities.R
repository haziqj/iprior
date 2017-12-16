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

#' Cut a numeric vector to a certain number of decimal places
#'
#' @param x A numeric vector.
#' @param k The number of decimal places.
#'
#' @return A character vector with the correct number of decimal places.
#'
#' @examples
#' decimal_place(pi, 3)
#' decimal_place(c(exp(1), pi, sqrt(2)), 4)
#'
#' @export
decimal_place <- function(x, k = 2) format(round(x, k), nsmall = k)

#' @rdname decimal_place
#' @export
dec_plac <- decimal_place

#' Convert \code{difftime} class into \code{time} class
#'
#' @param x A \code{difftime} object.
#'
#' @return A \code{time} object which contains the time difference and units.
#'
#' @export
as.time <- function(x) {
  time <- as.numeric(x)
  unit <- attr(x, "units")
  structure(list(time = time, unit = unit), class = "time")
}

#' @export
print.time <- function(x, ...) {
  cat(x$time, x$unit)
}

#' Check the structure of the hyperparameters of an I-prior model
#'
#' @param object An \code{ipriorMod} object or an \code{ipriorKernel} object.
#'
#' @return A printout of the structure of the hyperparameters.
#'
#' @export
check_theta <- function(object) {
  if (is.ipriorMod(object)) object <- object$ipriorKernel
  if (is.ipriorKernel(object)) {
    res <- names(object$thetal$theta)

    ind.lam <- grep("lambda", res)
    if (length(ind.lam) == 1) res[ind.lam] <- "log(lambda)"

    ind.hur <- grep("hurst", res)
    res[ind.hur] <- paste0("qnorm(", res[ind.hur], ")")

    ind.len <- grep("lengthscale", res)
    res[ind.len] <- paste0("log(", res[ind.len], ")")

    ind.off <- grep("offset", res)
    res[ind.off] <- paste0("log(", res[ind.off], ")")

    ind.psi <- grep("psi", res)
    res[ind.psi] <- "log(psi)"

    if (length(res) == 0) {
      cat("none")
    }
    else {
      cat(paste0("theta consists of ", length(res), ":\n"))
      cat(paste(res, collapse = ", "))
    }
  }
}

check_and_get_ipriorKernel <- function(object, assign.to.env = FALSE) {
  # Helper function to check whether object is of ipriorMod or ipriorKernel
  # class, and if so, replaces the object in environment with ipriorKernel.
  #
  # Args: An ipriorMod or ipriorKernel object; logical assign.to.env.
  #
  # Returns: Replacement of object with ipriorKernel object if necessary, or
  # assignment of ipriorKernel object to environment.
  if (is.ipriorMod(object) | is.ipriorKernel(object$ipriorKernel)) {
    if (isTRUE(assign.to.env)) {
      list2env(object$ipriorKernel, parent.frame())
    } else {
      assign(deparse(substitute(object)), object$ipriorKernel,
             envir = parent.frame())
    }
  } else if (is.ipriorKernel(object)) {
    if (isTRUE(assign.to.env)) {
      list2env(object, parent.frame())
    } else {
      assign(deparse(substitute(object)), object, envir = parent.frame())
    }
  } else {
    stop("Input an I-prior object.", call. = FALSE)
  }
}

check_and_get_ipriorMod <- function(object, assign.to.env = FALSE) {
  # Helper function to check whether object is of ipriorMod class.
  #
  # Args: An ipriorMod or ipriorKernel object; logical assign.to.env.
  #
  # Returns: Nothing - just checks. Unless assign.to.env is TRUE.
  if (is.ipriorMod(object)) {
    if (isTRUE(assign.to.env)) list2env(object$ipriorKernel, parent.frame())
  } else {
    stop("Input an ipriorMod object.", call. = FALSE)
  }
}

#' Test \code{iprior} objects
#'
#' Test whether an object is an \code{ipriorMod}, \code{ipriorKernel}, or either
#' object with Nystrom method enabled.
#'
#' @param x An \code{ipriorMod} or \code{ipriorKernel} object.
#'
#' @return Logical.
#'
#' @name is.iprior_x
NULL

#' @rdname is.iprior_x
#' @export
is.ipriorMod <- function(x) inherits(x, "ipriorMod")

#' @rdname is.iprior_x
#' @export
is.ipriorKernel <- function(x) inherits(x, "ipriorKernel")

is.ipriorKernel_old <- function(x) inherits(x, "ipriorKernel_old")

is.ipriorKernel_nys <- function(x) {
  if (is.ipriorMod(x)) x <- x$ipriorKernel
  if (is.ipriorKernel(x)) {
    !is.null(x$nystroml)
  } else {
    return(FALSE)
  }
}

is.ipriorKernel_cv <- function(x) {
  if (is.ipriorMod(x)) x <- x$ipriorKernel
  if (is.ipriorKernel(x)) {
    return(!is.null(x$y.test) & !is.null(x$Xl.test))
  } else {
    return(FALSE)
  }
}

#' @export
.is.ipriorKernel_cv <- is.ipriorKernel_cv

#' @rdname is.iprior_x
#' @export
is.nystrom <- is.ipriorKernel_nys

is.categorical <- function(x) {
  # Checks whether iprobit fitting is possible. This just checks whether the
  # response variables were factor type.
  #
  # Args: An ipriorMod, ipriorKernel or even an iprobitMod_x object (see iprobit
  # package for details.)
  #
  # Returns: Logical.
  check_and_get_ipriorKernel(x)
  !is.null(x$y.levels)
}

#' @export
.is.categorical <- is.categorical

#' Test kernel attributes
#'
#' Test whether an object uses a specific type of kernel.
#'
#' @param x An \code{ipriorMod} object, \code{ipriorKernel} object, a kernel
#'   matrix generated from one of the \code{kern_x()} functions, or even simply
#'   just a character vector.
#'
#' @return Logical.
#'
#' @name is.kern_x
NULL

is.kern_type <- function(x, type) {
  if (is.ipriorMod(x)) kernel_type <- x$ipriorKernel$kernels
  else if (is.ipriorKernel(x)) kernel_type <- x$kernels
  else if (!is.null(attributes(x)$kernel)) kernel_type <- attributes(x)$kernel
  else kernel_type <- x
  if (!is.null(kernel_type)) grepl(type, kernel_type)
  else return(FALSE)
}

#' @rdname is.kern_x
#' @export
is.kern_linear <- function(x) is.kern_type(x, type = "linear")

#' @rdname is.kern_x
#' @export
is.kern_canonical <- is.kern_linear

#' @rdname is.kern_x
#' @export
is.kern_fbm <- function(x) is.kern_type(x, type = "fbm")

#' @rdname is.kern_x
#' @export
is.kern_pearson <- function(x) is.kern_type(x, type = "pearson")

#' @rdname is.kern_x
#' @export
is.kern_se <- function(x) is.kern_type(x, type = "se")

#' @rdname is.kern_x
#' @export
is.kern_poly <- function(x) is.kern_type(x, type = "poly")

is.theta_lambda <- function(x) {
  # Helper function to determine whether or not a given set of hyperparameter
  # consists only of lambdas.
  #
  # Args: an ipriorMod or ipriorKernel object.
  #
  # Returns: Logical.
  if (is.ipriorMod(x)) x <- x$ipriorKernel
  if (is.ipriorKernel(x)) {
    theta <- names(x$thetal$theta)
    any.hurst       <- any(grepl("hurst"      , theta))
    any.lengthscale <- any(grepl("lengthscale", theta))
    any.offset      <- any(grepl("offset"     , theta))
    any.lambda      <- any(grepl("lambda"     , theta))
    return(
      all(!any.hurst, !any.lengthscale, !any.offset, any.lambda)
    )
  } else {
    stop("Not an ipriorX object.", call. = FALSE)
  }
}

#' Emulate \code{ggplot2} default colour palette
#'
#' Emulate \code{ggplot2} default colour palette. \code{ipriorColPal} and
#' \code{ggColPal} are DEPRECATED.
#'
#' This is the default colour scale for categorical variables in \code{ggplot2}.
#' It maps each level to an evenly spaced hue on the colour wheel. It does not
#' generate colour-blind safe palettes.
#'
#' \code{ipriorColPal()} used to provide the colour palette for the
#' \code{iprior} package, but this has been changed \code{ggplot2}'s colour
#' palette instead.
#'
#' @param x The number of colours required.
#' @param h Range of hues to use, in [0, 360].
#' @param c Chroma (intensity of colour), maximum value varies depending on
#'   combination of hue and luminance.
#' @param l Luminance (lightness), in [0, 100].
#'
#' @export
gg_colour_hue <- function(x, h = c(0, 360) + 15, c = 100, l = 65) {
  hues <- seq(h[1], h[2], length = x + 1)
  grDevices::hcl(h = hues, c = c, l = l)[1:x]
}

#' @rdname gg_colour_hue
#' @export
gg_color_hue <- gg_colour_hue

#' @rdname gg_colour_hue
#' @export
gg_col_hue <- gg_colour_hue

#' @rdname gg_colour_hue
#' @export
ipriorColPal <- function(x) {
  warning("Deprecated. Use gg_colour_hue() instead.")
  gg_colour_hue(x)
}
# ipriorColPal <- function(x) {
#   colx <- c(RColorBrewer::brewer.pal(9, "Set1")[-9],
#             RColorBrewer::brewer.pal(8, "Dark2"))
#   colx[6] <- RColorBrewer::brewer.pal(8, "Set2")[6]
#   colx[x]
# }

#' @rdname gg_colour_hue
#' @export
ggColPal <- function(x) {
  warning("Deprecated. Use gg_colour_hue() instead.")
  gg_colour_hue(x)
}

get_y_and_levels <- function(y) {
  # Function used for categorical response model to obtains the levels in the
  # ys.
  #
  # Args: Categorical response variables y.
  #
  # Returns: A list of the numerical values (1, 2, 3, ...) and the levels.
  list(y = as.numeric(y), levels = levels(y))
}

#' @export
.checkLevels <- get_y_and_levels

#' @export
.get_y_and_levels <- get_y_and_levels

fix_call_default <- function(cl = match.call(), new.name = "iprior") {
  # Replace the default call name with a new name. When using the default call,
  # it is possible that some of the X names are blank. This fixes that too.
  #
  # Args: The call and the new.name.
  #
  # Returns: The fixed call.
  cl[[1L]] <- as.name(new.name)
  where.blanks <- grepl("^$", names(cl))[-(1:2)]
  names(cl)[-(1:2)][where.blanks] <- paste0("X", which(where.blanks))
  # names(cl)[2] <- ""  # get rid of "y ="
  cl
}

#' @export
.fix_call_default <- fix_call_default

fix_call_formula <- function(cl = match.call(), new.name = "iprior") {
  # Replace the formula call name with a new name.
  #
  # Args: The call and the new.name.
  #
  # Returns: The fixed call.
  cl[[1L]] <- as.name(new.name)
  cl
}

#' @export
.fix_call_formula <- fix_call_formula

formula_to_xy <- function(formula, data, one.lam) {
  # Convert formula entry to y, X entry.
  #
  # Args: The formula, data frame and a logical option for one.lam.
  #
  # Returns: A list containing the data y and X, interactions instructions, x
  # and y names, and model terms (tt).
  mf <- model.frame(formula = formula, data = data)
  tt <- terms(mf)
  Terms <- delete.response(tt)
  x <- model.frame(Terms, mf)
  y <- model.response(mf)
  xname <- names(x)
  yname <- names(attr(tt, "dataClasses"))[1]
  x <- as.list(x)
  attr(x, "terms") <- NULL
  # attr(y, "yname") <- yname

  # Interactions ---------------------------------------------------------------
  interactions <- NULL
  tmpo <- attr(tt, "order")
  tmpf <- attr(tt, "factors")
  tmpf2 <- as.matrix(tmpf[-1, tmpo == 2])  # this obtains 2nd order interactions
  int2 <- apply(tmpf2, 2, function(x) which(x == 1))
  if (any(tmpo == 2)) interactions <- int2
  intr.3plus <- NULL
  tmpf3 <- as.matrix(tmpf[-1, tmpo > 2])
  int3 <- apply(tmpf3, 2, where_int)
  if (any(tmpo > 2)) intr.3plus <- int3
  interactions <- list(intr = interactions, intr.3plus = intr.3plus)

  # Deal with one.lam option ---------------------------------------------------
  if (isTRUE(one.lam)) {
    # Writes x and xname to env.
    list2env(deal_with_one.lam(x, interactions), envir = environment())
  }

  list(y = y, Xl = x, interactions = interactions, xname = xname, yname = yname,
       tt = tt)
}

#' @export
.formula_to_xy <- formula_to_xy

deal_with_one.lam <- function(x, interactions) {
  # Helper function to convert list of X according to one.lam option.
  #
  # Args: List of covariates and interactions information from ipriorKernel.
  #
  # Returns: Update Xl.
  xname <- names(x)
  xnl <- length(xname)
  if (!all(sapply(interactions, is.null))) {
    stop("Cannot use option one.lam = TRUE with interactions.", call. = FALSE)
  }
  if (length(x) == 1) {
    message("Option one.lam = TRUE used with a single covariate anyway.")
  }
  attributes(x)$terms <- attributes(x)$names <- NULL
  if (xnl <= 3) {
    xname <- paste(xname, collapse = " + ")
  } else {
    xname <- paste(xname[1], "+ ... +", xname[xnl])
  }
  x <- list(matrix(unlist(x), ncol = length(x)))
  names(x) <- xname
  attr(x, "one.lam") <- TRUE

  list(x = x, xname = xname)
}

#' @export
.deal_with_one.lam <- deal_with_one.lam

terms_to_xy <- function(object, newdata) {
  # Args: An ipriorKernel object.
  tt <- object$terms
  Terms <- delete.response(tt)
  x <- model.frame(Terms, newdata)
  y <- NULL
  if (any(colnames(newdata) == object$yname))
    y <- model.extract(model.frame(tt, newdata), "response")

  # Deal with one.lam option ---------------------------------------------------
  one.lam <- attr(object$Xl, "one.lam")
  if (isTRUE(one.lam)) {
    # Writes x and xname to env.
    list2env(deal_with_one.lam(x, object$interactions), envir = environment())
  }

  list(Xl = x, y = y)
}

#' @export
.terms_to_xy <- terms_to_xy

fastSquareRoot2 <- function(x) {
  # Function to quickly find a square root of a matrix from its
  # eigendecomposition.
  #
  # Args: x a square matrix.
  #
  # Returns: x ^ {1/2}.
  tmp <- eigenCpp(x)
  tmp$vec %*% tcrossprod(diag(sqrt(abs(tmp$val))), tmp$vec)
}

.onUnload <- function(libpath) {
  # Whenever you use C++ code in your package, you need to clean up after
  # yourself when your package is unloaded.
  library.dynam.unload("iprior", libpath)
}
