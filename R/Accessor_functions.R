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

#' Accessor functions for \code{ipriorMod} objects.
#'
#' @param object An \code{ipriorMod} object.
#' @param units Units for object size.
#' @param standard Standard for object size.
#' @param theta (Optional) Value of hyperparameters to evaluate the kernel
#'   matrix.
#' @param newdata (Optional) If not supplied, then a square, symmetric kernel
#'   matrix is returned using the data as input points. Otherwise, the kernel
#'   matrix is evaluated with respect to this set of data as well. It must be a
#'   list of vectors/matrices with similar dimensions to the original data.
#' @param error.type (Optional) Report the mean squared error of prediction
#'   (\code{"MSE"}), or the root mean squared error of prediction
#'   (\code{"RMSE"})
#'
#' @name Accessors
NULL

#' @describeIn Accessors Obtain the intercept.
#' @export
get_intercept <- function(object) {
  check_and_get_ipriorKernel(object)
  object$intercept
}

#' @describeIn Accessors Obtain the response variables.
#' @export
get_y <- function(object) {
  check_and_get_ipriorKernel(object)
  if (is.categorical(object)) {
    warning("Categorical variables were numerised.", call. = FALSE)
  }
  res <- as.numeric(object$y + get_intercept(object))
  names(res) <- rownames(object$y)
  res
}

#' @describeIn Accessors Obtain the object size of the I-prior model.
#' @export
get_size <- function(object, units = "kB", standard = "SI") {
  check_and_get_ipriorKernel(object)
  print(object.size(object), units = units, standard = standard)
}

#' @describeIn Accessors Obtain the hyerparameters of the model (both estimated and fixed ones).
#' @export
get_hyp <- function(object) {
  check_and_get_ipriorMod(object)
  res <- object$param.full
  if (length(res) > 0) return(res)
  else cat("NA")
}

#' @describeIn Accessors Obtain the scale parameters used.
#' @export
get_lambda <- function(object) {
  tmp <- get_hyp(object)
  res <- tmp[grep("lambda", names(tmp))]
  if (length(res) > 0) return(res)
  else cat("NA")
}

#' @describeIn Accessors Obtain the error precision.
#' @export
get_psi <- function(object) {
  tmp <- get_hyp(object)
  res <- tmp[grep("psi", names(tmp))]
  if (length(res) > 0) return(res)
  else cat("NA")
}

#' @describeIn Accessors Obtain the lengthscale for the SE kernels used.
#' @export
get_lengthscale <- function(object) {
  tmp <- get_hyp(object)
  res <- tmp[grep("lengthscale", names(tmp))]
  if (length(res) > 0) return(res)
  else cat("NA")
}

#' @describeIn Accessors Obtain the Hurst coefficient of the fBm kernels used.
#' @export
get_hurst <- function(object) {
  tmp <- get_hyp(object)
  res <- tmp[grep("hurst", names(tmp))]
  if (length(res) > 0) return(res)
  else cat("NA")
}

#' @describeIn Accessors Obtain the offset parameters for the polynomial kernels used.
#' @export
get_offset <- function(object) {
  tmp <- get_hyp(object)
  res <- tmp[grep("offset", names(tmp))]
  if (length(res) > 0) return(res)
  else cat("NA")
}

#' @describeIn Accessors Obtain the degree of the polynomial kernels used.
#' @export
get_degree <- function(object) {
  if (is.kern_poly(object)) {
    tmp <- get_kernels(object)
    return(get_polydegree(tmp))
  } else {
    cat("NA")
  }
}

#' @describeIn Accessors Obtain the standard errors of the estimated hyperparameters.
#' @export
get_se <- function(object) {
  check_and_get_ipriorMod(object)
  expand_theta(object$se, object$ipriorKernel$thetal$theta.drop, NA)
}

#' @describeIn Accessors Obtain the kernels used.
#' @export
get_kernels <- function(object) {
  if (is.ipriorMod(object)) theta <- object$theta
  if (is.ipriorKernel(object)) theta <- object$thetal$theta
  check_and_get_ipriorKernel(object)
  param.tab <- theta_to_param(theta, object)
  res <- param.tab$kernel
  names(res) <- object$xname
  res
}

#' @describeIn Accessors Obtain the kernel matrix of the I-prior model.
#' @export
get_kern_matrix <- function(object, theta = NULL, newdata) {
  if (is.ipriorMod(object)) {
    # estl <- object$ipriorKernel$estl
    # til.cond <- (
    #   !isTRUE(estl$est.hurst) & !isTRUE(estl$est.lengt) & !isTRUE(estl$est.offs)
    # )
    if (missing(newdata)) {
      res <- get_Hlam(object$ipriorKernel, object$theta, FALSE)
    } else {
      list2env(before_predict(object, newdata), envir = environment())  # writes xstar to env
      res <- get_Htildelam(object$ipriorKernel, object$theta, xstar)
    }
    return(res)
  } else if (is.ipriorKernel(object)) {
    # estl <- object$estl
    # til.cond <- (
    #   !isTRUE(estl$est.hurst) & !isTRUE(estl$est.lengt) & !isTRUE(estl$est.offs)
    # )
    if (missing(newdata)) {
      res <- get_Hlam(object, object$thetal$theta, FALSE)
    } else {
      list2env(before_predict(list(ipriorKernel = object), newdata),
               envir = environment())  # writes xstar to env
      res <- get_Htildelam(object, object$thetal$theta, xstar)
    }
    return(res)
  }
}

#' @describeIn Accessors Obtain the training mean squared error.
#' @export
get_prederror <- function(object, error.type = c("RMSE", "MSE")) {
  check_and_get_ipriorMod(object)
  error.type <- match.arg(toupper(error.type), c("RMSE", "MSE"))
  if (error.type == "MSE") fun <- function(x) x
  if (error.type == "RMSE") fun <- function(x) sqrt(x)
  train.error <- c("Training" = object$train.error)
  if (is.ipriorKernel_cv(object)) {
    res <- fun(c(train.error, "Test" = object$test$test.error))
  } else {
    res <- fun(train.error)
  }
  names(res) <- paste(names(res), error.type)
  res
}

get_mse <- function(object) get_prederror(object, "MSE")

get_rmse <- function(object) get_prederror(object, "RMSE")

#' @describeIn Accessors Obtain information on which hyperparameters were
#'   estimated and which were fixed.
#' @export
get_estl <- function(object) {
  check_and_get_ipriorKernel(object)
  unlist(object$estl)
}

#' @describeIn Accessors Obtain the estimation method used.
#' @export
get_method <- function(object) {
  check_and_get_ipriorMod(object)
  cat(object$est.method)
}

#' @describeIn Accessors Obtain the convergence information.
#' @export
get_convergence <- function(object) {
  check_and_get_ipriorMod(object)
  cat(object$est.conv)
}

#' @describeIn Accessors Obtain the number of iterations performed.
#' @export
get_niter <- function(object) {
  check_and_get_ipriorMod(object)
  niter <- object$niter
  maxit <- object$control$maxit
  cat("Iterations:", paste0(niter, "/", maxit, "."))
}

#' @describeIn Accessors Obtain the time taken to complete the estimation
#'   procedure.
#' @export
get_time <- function(object) {
  check_and_get_ipriorMod(object)
  object$time
}

#' @describeIn Accessors Extract the theta value at convergence. Note that this
#'   is on an unrestricted scale (see the vignette for details).
#' @export
get_theta <- function(object) {
  check_and_get_ipriorMod(object)
  object$theta
}
