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

iprior_method_checker <- function(object, method) {
  # Depending on the I-prior model, only certain types of estimation methods are
  # available. E.g. if estimating Hurst parameter for fBm kernel, then the EM
  # closed-form is not possible.
  #
  # Args: An ipriorKernel object, and the method input by the user in iprior().
  #
  # Returns: A logical vector of length seven, with only one entry being TRUE,
  # directing the iprior() function to the correct iprior_x() estimation
  # routine.
  res <- rep(FALSE, 7)
  names(res) <- c("fixed", "canonical", "em.closed", "em.reg", "direct",
                  "nystrom", "mixed")

  if (object$thetal$n.theta == 0 | method == "fixed") {
    res["fixed"] <- TRUE
    # if (method != "fixed") warning("No hyperparameters estimated. Using fixed estimation method.", call. = FALSE)
  } else if (is.ipriorKernel_nys(object)) {
    res["nystrom"] <- TRUE
  } else if (method == "canonical") {
    if (all(is.kern_linear(object$kernels))) {
      res["canonical"] <- TRUE
    } else {
      res["direct"] <- TRUE
      warning("Non-linear kernels used. Using direct estimation method.", call. = FALSE)
    }
    if (object$no.int > 0) {
      res["direct"] <- TRUE
      warning("Canonical method not possible with interactions. Using direct estimation method.", call. = FALSE)
    }
  } else if (method == "em") {
    if (is.null(object$BlockBStuff)) {
      res["em.reg"] <- TRUE
    } else {
      if (!is.null(object$y.levels)) res["em.closed"] <- TRUE
      else if (!isTRUE(object$estl$est.lambda) | !isTRUE(object$estl$est.psi))
        res["direct"] <- TRUE
      else
        res["em.closed"] <- TRUE
    }
  } else if (method == "mixed") {
    res["mixed"] <- TRUE
  } else {
    res["direct"] <- TRUE
  }

  res
}

#' @export
.iprior_method_checker <- iprior_method_checker

update_control <- function(arg.list, default.list) {
  # In most of the iprior_x functions, a list of control options is required.
  # This helper function checks whether the list supplied contains the allowed
  # set of control options, and if so, overwrites them.
  #
  # Args: arg.list is the control list supplied from the function, and
  # default.list is the allowed list of controls.
  #
  # Returns: An appropriate list.
  default_names <- names(default.list)
  default.list[(arg_names <- names(arg.list))] <- arg.list
  if (length(noNms <- arg_names[!arg_names %in% default_names])) {
    warning("Unknown names in control options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  default.list
}

#' @export
.update_control <- update_control

convert_se <- function(se, theta, object) {
  # Converts the standard errors for theta into the standard errors for param
  # using the delta method. These are then used in the summary() method for
  # ipriorMod objects.
  #
  # Args: se for theta, the optimised theta value, and the ipriorKernel object.
  #
  # Returns: A vector similar to se, but with delta-transformed values.
  theta.names <- names(object$thetal$theta)
  res <- se

  types <- c("lambda", "hurst", "offset", "lengthscale", "psi")
  for (i in seq_along(types)) {
    type <- types[i]
    ind <- grep(type, theta.names)
    if (type == "lambda" & length(ind) == 1)
      res[ind] <- se[ind] * exp(theta[ind])
    else if (type == "hurst")
      res[ind] <- se[ind] * dnorm(theta[ind])
    else
      res[ind] <- se[ind] * exp(theta[ind])
  }

  res
}

#' @export
.convert_se <- convert_se

get_Hlam <- function(object, theta, theta.is.lambda = FALSE) {
  # Obtain the scaled kernel matrix Hlam.
  #
  # Args: An ipriorKernel object, theta values to calculate Hlam at, and a
  # logical argument theta.is.lambda. If theta.is.lambda = TRUE, then it is just
  # a matter of taking the sumproduct of Hl and theta. Otherwise, Hl will be
  # re-calculated everytime this function is called.
  #
  # Returns: The scaled kernel matrix. This has a "kernel" attribute indicating
  # which kernels were used to generate it.
  if (isTRUE(theta.is.lambda)) {
    kernels <- object$kernels
    lambda <- theta
    lambda.only <- TRUE
  } else {
    tmp <- theta_to_param(theta, object)
    kernels <- tmp$kernels
    lambda <- tmp$lambda
    lambda.only <- is.theta_lambda(object)
  }

  if (isTRUE(lambda.only)) {
    Hl <- object$Hl
  } else {
    if (is.nystrom(object)) {
      Hl <- get_Hl(object$Xl, get_Xl.nys(object), kernels = kernels,
                   lambda = lambda)
    } else {
      Hl <- get_Hl(object$Xl, list(NULL), kernels = kernels, lambda = lambda)
    }
  }
  lambda[is.kern_poly(kernels)] <- 1
  calc_Hlam(Hl, lambda[seq_along(Hl)], object)
}

get_Htildelam <- function(object, theta, xstar, theta.is.lambda = FALSE) {
  # A helper function to calculate Hlam given a new set of data. This is similar
  # to get_Hlam().
  #
  # Args: An ipriorKernel object, theta values to calculate Hlam at, and a list
  # of new data.
  #
  # Returns: Assuming m new data points, then an m x n matrix is returned.
  if (isTRUE(theta.is.lambda)) {
    kernels <- object$kernels
    lambda <- theta
  } else {
    tmp <- theta_to_param(theta, object)
    kernels <- tmp$kernels
    lambda <- tmp$lambda
  }

  # if (is.nystrom(object)) {
  #   m <- object$nystroml$nys.size
  #   Hlam <- get_Hlam(object, theta)  # [A B]
  #   A <- Hlam[1:m, 1:m]
  #   Hl <- get_Hl(get_Xl.nys(object), xstar, kernels, lambda)
  #   Hlam.new <- calc_Hlam(Hl, lambda, object)
  #   res <- Hlam.new %*% solve(A, Hlam)
  # } else {
    Hl <- get_Hl(object$Xl, xstar, kernels = kernels, lambda = lambda)
    lambda[is.kern_poly(kernels)] <- 1
    res <- calc_Hlam(Hl, lambda, object)
  # }

  res
}

calc_Hlam <- function(Hl, lambda, object) {
  # Helper function to get the sumproduct of the list of kernel matrices and the
  # scale parameters.
  #
  # Args: Hl list of kernel matrices, scale parameters lambda, and ipriorKernel
  # object.
  #
  # Returns: Hlam matrix.
  expand_Hl_and_lambda(Hl, lambda, object$intr, object$intr.3plus,
                       environment())
  res <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
  kernels.to.add <- sapply(Hl, attr, "kernel")
  attr(res ,"kernel") <- paste(kernels.to.add, collapse = " + ")
  res
}

#' @export
.get_Hlam <- get_Hlam

#' @export
.get_Htildelam <- get_Htildelam

#' @export
.calc_Hlam <- calc_Hlam

eigen_Hlam <- function(Hlam, env = NULL) {
  # The routine for the eigendecomposition of the scaled kernel matrix.
  #
  # Args: The scaled kernel matrix Hlam, and the optional environment in which
  # to assign the results.
  #
  # Returns: A list of the eigen values (u) and vectors (V) of Hlam if env =
  # NULL, otherwise u and V are assigned to env.
  tmp <- eigenCpp(Hlam)
  res <- list(u = tmp$val, V = tmp$vec)
  if (is.null(env)) return(res)
  else list2env(res, env)
}

#' @export
.eigen_Hlam <- eigen_Hlam

eigen_Hlam_nys <- function(Hlam, env = NULL) {
  # The routine for the eigendecomposition of the scaled kernel matrix using
  # Nystrom approximation.
  #
  # Args: The (truncated) scaled kernel matrix Hlam, which is m x n in size,
  # where m = object$nys.size
  #
  # Returns: A list of the eigen values (u) and vectors (V) of Hlam if env =
  # NULL, otherwise u and V are assigned to env.
  m <- nrow(Hlam)
  A <- Hlam[, seq_len(m)]
  B <- Hlam[, -seq_len(m)]
  tmp1 <- eigenCpp(A)
  # print(c(lambda, psi))
  # if (any(tmp1$val + 1e-7 < 0)) print(tmp1$val + 1e-7)
  U <- tmp1$vectors
  C.tmp <- U * rep(1 / sqrt(tmp1$val + 1e-7), each = nrow(U))
  C <- C.tmp %*% crossprod(U, B)
  Q <- A + tcrossprod(C)
  tmp2 <- eigenCpp(Q)
  # if (any(tmp2$val + 1e-7 < 0)) print(tmp2$val + 1e-7)
  u <- tmp2$values + 1e-7
  R <- tmp2$vectors
  V <- rbind(A, t(B)) %*% tcrossprod(C.tmp, U) %*%
    (R * rep(1 / sqrt(u), each = nrow(R)))

  res <- list(u = u, V = V)
  if (is.null(env)) return(res)
  else list2env(res, env)
}

#' @export
.eigen_Hlam_nys <- eigen_Hlam_nys

A_times_a <- function(u, V, a) {
  # Calculate A %*% a from the eigendecomposition of A. Mostly used to calculate
  # Vy^{-1} %*% a without having to invert Vy, because the eigendecompostion of
  # Vy would have already been obtained.
  #
  # Args: Eigenvalues u, eigenvectors V (of some matrix A), and a vector a.
  #
  # Returns: A vector A %*% a.
  (V * rep(u, each = nrow(V))) %*% crossprod(V, a)
}

#' @export
.A_times_a <- A_times_a

em_loop_logical <- function() {
  # Helper function to determine when to stop the while loop for the EM
  # algorithm.
  #
  # Args: none, but reads loglik, niter, maxit, stop.crit from the parent
  # environment.
  #
  # Returns: If niter == 0 return TRUE because must complete 1 iteration. If
  # niter == 1 then stop if maxit == 1, otherwise continue (nothing to compare).
  # If niter > 1 then just check whether maxit reached or stop.crit reached.
  ll.diff <- loglik[niter] - loglik[niter - 1]
  crit1 <- (niter != maxit)
  crit2 <- (abs(ll.diff) > stop.crit)
  if (niter == 0) {
    return(TRUE)
  } else if (niter == 1) {
    return(crit1)
  } else {
    if (ll.diff < 0) {
      warning(paste0("Log-likelihood decreased at iteration ", niter),
              call. = FALSE)
    }
    return(crit1 & crit2)
  }
}

#' @export
.em_loop_logical <- em_loop_logical

get_w <- function(u, V, Vy.inv.y, psi) {
  # Helper function to obtain posterior mean of I-prior random effects. It is
  # calculated as psi * Hlam %*% Vy.inv %*% y.
  #
  # Args: u and V are the eigendecomposition of Hlam, Vy.inv.y and psi are
  # self-explanatory.
  #
  # Returns: Numeric w.
  as.numeric(psi * A_times_a(u, V, Vy.inv.y))
}

get_y.hat <- function(u, V, w) {
  # Obtain fitted values after estimation procedures. This is calculated as
  # y.hat = Hlam %*% w. Note that intercepts have NOT been added.
  #
  # Args: u and V are the eigendecomposition of Hlam, and the rest are
  # self-explanatory.
  #
  # Returns: Numeric y.hat.
  as.numeric(A_times_a(u, V, w))
}
