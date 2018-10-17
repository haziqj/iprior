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

# The main difficulty with estimating I-prior models is the number of
# hyperparameters to estimate depend on which kernels are used. In kernL2(),
# kernels are specified by the user, which are then applied to the data x. From
# this, we create a "param" table which specifies the kernels used for each
# variable. The param table must be converted to the "theta" vector which is the
# vector of parameters to estimate. This file contains helper functions to
# convert between kernels, param, and theta. In particular, we have: 1.
# param_to_kernel(); 2. theta_to_kernel(); 3. param_to_theta(); 4.
# kernel_to_param(); 5. theta_to_param(); 6. theta_to_collapsed_param(); 7.
# theta_to_psi().
#
# We then also have several other helper functions, in particular
# expand_Hl_and_lambda() and get_Hl() which are also used in predict().

param_to_kernel <- function(param) {
  # Convert param table to vector of kernels.
  #
  # Args: param table.
  #
  # Returns: Vector of kernels.
  as.character(param$kernels)
}

#' @export
.param_to_kernel <- param_to_kernel

theta_to_kernel <- function(theta, object) {
  # Convert theta to vector of kernels.
  #
  # Args: theta and an ipriorKernel object.
  #
  # Returns: Vector of kernels.
  param <- theta_to_param(mod$theta, mod$kernL)
  param_to_kernel(param)
}

#' @export
.theta_to_kernel <- theta_to_kernel

param_to_theta <- function(param, est.list, logpsi = 0) {
  # Convert param table to the theta vector. Note that theta is designed so that
  # the values are unbounded, i.e. hurst is Phi^{-1}(hurst), lengthscale is
  # log(lengthscale), etc. theta_to_param() reverses this for final
  # presentation.
  #
  # Args: A param table, the list of parameters to be estimated and optional
  # logpsi value.
  #
  # Returns: theta is a vector of parameters to be passed to optim or EM for
  # optimisation, including the logpsi value.
  param <- param[, seq_len(4)]  # lambda, hurst, lengthscale, offset
  param$hurst <- qnorm(param$hurst)
  param$lengthscale <- log(param$lengthscale)
  param$offset <- log(param$offset)
  if (nrow(param) == 1) param$lambda <- log(param$lambda)

  tmp <- collapse_param(param)
  theta.full <- c(tmp$param, psi = logpsi)
  param.na <- tmp$na
  tmp <- reduce_theta(theta.full, est.list)
  theta.reduced <- tmp$theta.reduced
  theta.drop <- tmp$theta.drop
  theta.omitted <- tmp$theta.omitted

  list(theta = theta.reduced, param.na = param.na, theta.drop = theta.drop,
       theta.omitted = theta.omitted)
}

#' @export
.param_to_theta <- param_to_theta

kernel_to_param <- function(kernels, lambda) {
  # Convert vector of kernels to a param table.
  #
  # Args: kernels is a p-vector of kernels to apply on the data. lambda are the
  # scale parameters.
  #
  # Returns: The param table.
  param <- as.data.frame(matrix(NA, ncol = 5, nrow = length(kernels)))
  names(param) <- c("lambda", "hurst", "lengthscale", "offset", "degree")
  param$lambda <- lambda
  kernel.match <- c("fbm", "se", "poly", "poly")
  for (i in seq_along(kernel.match)) {
    res <- rep(NA, length(kernels))
    where_kernel <- grepl(kernel.match[i], kernels)
    for (j in seq_len(sum(where_kernel))) {
      res[where_kernel][j] <- get_hyperparam(kernels[where_kernel][j])
      if (i == 4) res[where_kernel][j] <- get_polydegree(kernels[where_kernel][j])
    }
    param[, i + 1] <- res
  }

  res <- cbind(param, kernels)
  res$kernels <- as.character(res$kernels)
  res
}

#' @export
.kernel_to_param <- kernel_to_param

reduce_theta <- function(theta.full, est.list) {
  # The user may specify for some of the hyperparameters to not be estimated,
  # therefore theta vector should be smaller than the full set of
  # hyperparameters. This function trims theta correctly.
  #
  # Args: theta.full is the original and full theta vector. est.list is
  # generated in kernL2() which is a list telling us which parameters are
  # estimated and which are not.
  #
  # Returns: The trimmed theta vector.
  theta.full.orig <- theta.full

  # Estimate Hurst coefficient? ------------------------------------------------
  est.lambda <- est.list$est.lambda
  ind.lambda <- grepl("lambda", names(theta.full))
  theta.full[ind.lambda][!est.lambda] <- NA

  # Estimate Hurst coefficient? ------------------------------------------------
  est.hurst <- est.list$est.hurst
  ind.hurst <- grepl("hurst", names(theta.full))
  theta.full[ind.hurst][!est.hurst] <- NA

  # Estimate lengthscale in SE kernel? -----------------------------------------
  est.l <- est.list$est.lengthscale
  ind.l <- grepl("lengthscale", names(theta.full))
  theta.full[ind.l][!est.l] <- NA

  # Estimate offset in polynomial kernel? --------------------------------------
  est.c <- est.list$est.offset
  ind.c <- grepl("offset", names(theta.full))
  theta.full[ind.c][!est.c] <- NA

  # Estimate error precision psi? ----------------------------------------------
  est.psi <- est.list$est.psi
  ind.psi <- grepl("psi", names(theta.full))
  theta.full[ind.psi][!est.psi] <- NA

  theta.drop <- is.na(theta.full)
  theta.reduced <- theta.full[!theta.drop]
  theta.omitted <- theta.full.orig[theta.drop]

  list(theta.reduced = theta.reduced, theta.omitted = theta.omitted,
       theta.drop = theta.drop)
}

#' @export
.reduce_theta <- reduce_theta

expand_theta <- function(theta.reduced, theta.drop, theta.omitted) {
  # This convertes the reduced theta vector back to the full theta vector.
  # Useful in other functions such as theta_to_param().
  #
  # Args: theta.reduced is obvious. theta.drop is information on which of the
  # full theta was not estimated. theta.omitted contains the values for which
  # theta component that was not estimated. These are all generated by
  # param_to_theta().
  #
  # Returns: The full theta.
  theta.full <- theta.drop
  theta.full[!theta.drop] <- theta.reduced
  theta.full[theta.drop] <- theta.omitted
  theta.full
}

#' @export
.expand_theta <- expand_theta

collapse_param <- function(param) {
  # Args: A param table.
  #
  # Returns: A vectorised form of param with the na values removed. The param.na
  # values are output here too.
  #
  # Notes: Used as a helper function in param_to_theta(), and also useful for
  # final presentation of the parameters as this function names the parameters
  # too.
  res <- na.omit(unlist(param[, 1:4]))
  param.names <- names(res)
  na <- as.numeric(na.action(res))
  res <- as.numeric(res)

  param.digits <- gsub("[^[:digit:]]", "", param.names)
  param.names <- gsub("[[:digit:]]", "", param.names)
  param.names <- paste0(param.names, "[", param.digits, "]")
  param.names <- gsub("[[]]", "", param.names)
  names(res) <- param.names

  list(param = res, na = na)
}

#' @export
.collapse_param <- collapse_param

theta_to_param <- function(theta, object) {
  # Args: A vector of parameters to be optimised, including logpsi. object must
  # be either a ipriorKernel type object, or a list containing param.na,
  # which.pearson and poly.deg.
  #
  # Returns: A param table.
  #
  # Notes: The logpsi value is removed. To obtain this use theta_to_psi(). If
  # object is specified, then the param.na, which.pearson and poly.deg are
  # obtained from object.
  param.na <- object$thetal$param.na
  which.pearson <- object$which.pearson
  poly.deg <- object$poly.deg
  theta <- expand_theta(theta, object$thetal$theta.drop,
                        object$thetal$theta.omitted)
  theta <- theta[-length(theta)]

  full.length <- length(c(theta, param.na))
  param <- matrix(NA, ncol = 4, nrow = full.length / 4)
  tmp <- c(param)
  tmp[-param.na] <- theta
  param[] <- tmp
  param <- cbind(param, degree = poly.deg)

  param <- as.data.frame(param)
  names(param) <- c("lambda", "hurst", "lengthscale", "offset", "degree")

  param$hurst <- pnorm(param$hurst)
  param$lengthscale <- exp(param$lengthscale)
  param$offset <- exp(param$offset)
  if (nrow(param) == 1) param$lambda <- exp(param$lambda)
  param$kernels <- correct_pearson_kernel(
    apply(param, 1, param_translator), which.pearson
  )

  param
}

#' @export
.theta_to_param <- theta_to_param

param_translator <- function(x) {
  # Helper function in theta_to_param().
  #
  # Args: Row vector from param table.
  #
  # Returns: The kernel used.
  hyperparam <- x[-1]
  if (!is.na(hyperparam[1]))
    return(paste0("fbm,", hyperparam[1]))
  if (!is.na(hyperparam[2]))
    return(paste0("se,", hyperparam[2]))
  if (!is.na(hyperparam[3]))
    return(paste0("poly", hyperparam[4], ",", hyperparam[3]))
  "linear"
}

#' @export
.param_translator <- param_translator

correct_pearson_kernel <- function(x, which.pearson) {
  # When using theta_to_param(), unable to identify which data x uses the
  # Pearson kernel. This helper function corrects it by reading from the logical
  # which.pearson vector.
  #
  # Args: The kernel vector and which.pearson (logical), indicating which of the
  # x position uses the Pearson kernel.
  #
  # Returns: The corrected kernel vector.
  x[which.pearson] <- "pearson"
  x
}

#' @export
.correct_pearson_kernel <- correct_pearson_kernel

kernel_translator <- function(x, y = NULL, kernel, lam.poly = 1) {
  # Used as a helper function in get_Hl() to output list of kernel
  # matrices in kernL2() and predict(). For future expansion, add new kernels
  # here.
  #
  # Args: x, y (optional) data and kernel a character vector indicating which
  # kernel to apply x and y on. lam.poly is the scale for polynomial kernels.
  #
  # Returns: A kernel matrix.
  if (grepl("linear", kernel)) return(kern_linear(x, y))
  if (grepl("canonical", kernel)) return(kern_linear(x, y))
  if (grepl("fbm", kernel)) {
    if (grepl(",", kernel)) {
      hurst <- get_hyperparam(kernel)
      return(kern_fbm(x, y, gamma = hurst))
    } else {
      return(kern_fbm(x, y))
    }
  }
  if (grepl("se", kernel)) {
    if (grepl(",", kernel)) {
      lengthscale <- get_hyperparam(kernel)
      return(kern_se(x, y, l = lengthscale))
    } else {
      return(kern_se(x, y))
    }
  }
  if (grepl("poly", kernel)) {
    if (grepl(",", kernel)) {
      offset <- get_hyperparam(kernel)
    } else {
      offset <- 0
    }
    degree <- get_polydegree(kernel)
    return(kern_poly(x, y, c = offset, d = degree, lam.poly = lam.poly))
  }
  if (grepl("pearson", kernel)) return(kern_pearson(x, y))

  stop("Incorrect kernel specification or unsupported kernel.",
       call. = FALSE)
}

#' @export
.kernel_translator <- kernel_translator

theta_to_collapsed_param <- function(theta, object) {
  # This is a wrapper function for theta_to_param(). It is useful to get
  # the vector of parameters directly from theta.
  #
  # Args: theta (usually theta.full) and the ipriorKernel object.
  #
  # Returns: A vector of parameters.
  param <- theta_to_param(theta, object)
  c(collapse_param(param)$param, theta_to_psi(theta, object))
}

#' @export
.theta_to_collapsed_param <- theta_to_collapsed_param

theta_to_psi <- function(theta, object) {
  # Obtains psi from the theta vector, or if not estimated, from the
  # ipriorKernel object.
  #
  # Args: A vector of parameters to be optimised, including logpsi, and an
  # ipriorKernel object.
  #
  # Returns: psi, the error precision.
  theta <- expand_theta(theta, object$thetal$theta.drop,
                        object$thetal$theta.omitted)
  logpsi <- theta[length(theta)]
  exp(logpsi)
}

get_hyperparam <- function(x) {
  # Obtain the hyperparameters of kernels. Note: the exported version is
  # get_hyp().
  #
  # Args: x is either a character vector of the kernel specification (e.g. x <-
  # "fbm,0.5"), or it is the output from one of the kern_x() functions.
  #
  # Returns: The hyperparameter associated with the kernel (value that comes
  # after the comma).
  if (!is.null(attributes(x)$kernel)) x <- attributes(x)$kernel
  as.numeric(unlist(strsplit(x, ","))[2])
}

#' @export
.get_hyperparam <- get_hyperparam

get_polydegree <- function(x) {
  # Obtain the degree of the polynomial kernel. Note: the exported version is
  # get_degree().
  #
  # Args: x is either a character vector of the kernel specification (e.g. x <-
  # "fbm,0.5"), or it is the output from one of the kern_x() functions.
  #
  # Returns: Integer.
  if (!is.null(attributes(x)$kernel)) x <- attributes(x)$kernel
  degree <- unlist(strsplit(x, ","))[1]
  degree <- unlist(strsplit(degree, "poly"))
  if (length(degree) == 1) degree <- 2
  else degree <- as.numeric(degree[2])
  degree
}

#' @export
.get_polydegree <- get_polydegree

get_kernels_from_Hl <- function(x) {
  # Obtain the kernels used to generate the list of kernel matrices.
  #
  # Args: x is a list of kernel matrices.
  #
  # Returns: A character vector of kernels.
  sapply(x, function(x) attributes(x)$kernel)
}

#' @export
.get_kernels_from_Hl <- get_kernels_from_Hl

get_Xl.nys <- function(object) {
  # When using the Nystrom method, we calculate the kernel matrix using
  # get_Hl(Xl, Xl.nys) which produces a smaller m x n matrix. This is the helper
  # function to obtain Xl.nys.
  #
  # Args: An ipriorKernel object with Nystrom option called.
  #
  # Returns: Xl.nys, a list similar to Xl but only the first m rows are returned
  # (after sampling that is).
  nys.check <- is.ipriorKernel_nys(object)
  if (isTRUE(nys.check)) {
    Xl.nys <- lapply(object$Xl, reorder_x,
                     smp = seq_len(object$nystroml$nys.size))
    mostattributes(Xl.nys) <- attributes(object$Xl)
    return(Xl.nys)
  } else {
    stop("Nystrom option not called.", call. = FALSE)
  }
}

#' @export
.get_Xl.nys <- get_Xl.nys

get_Hl <- function(Xl, yl = list(NULL), kernels, lambda) {
  # Obtain the list of kernel matrices. Except for polynomial kernels, these are
  # not scaled, i.e. not Hlam matrices.
  #
  # Args: List of data Xl and yl (optional), vector of same length of kernel
  # characters to instruct kernel_translator() which kernels to apply each of
  # the x and y. lambda is needed for the polynomial kernels.
  #
  # Returns: List of kernel matrices.
  mapply(kernel_translator, Xl, yl, kernels, lambda, SIMPLIFY = FALSE)
}

#' @export
.get_Hl <- get_Hl

expand_Hl_and_lambda <- function(Hl, lambda, intr, intr.3plus, env = NULL) {
  # Helper function to expand Hl (list of kernel matrices) and lambda (scale
  # parameters) according to any interactions specification.
  #
  # Args: List of kernel matrices (Hl) and vector of scale parameters (lambda).
  # intr and intr.3plus are obtained from kernL2 which specifies the
  # interactions. env is the environment in which to assign the results.
  #
  # Returns: A list of expanded Hl and lambda. If env = NULL, then this list is
  # returned. If an environment is specified, it is directly assigned there.
  p <- length(Hl)
  no.int <- max(0, ncol(intr))
  no.int.3plus <- max(0, ncol(intr.3plus))
  if (!is.null(intr) & length(intr) > 0) {
    for (j in seq_len(no.int)) {
      ind1 <- intr[1, j]; ind2 <- intr[2, j]
      lambda[p + j] <- lambda[ind1] * lambda[ind2]
      Hl[[p + j]] <- Hl[[ind1]] * Hl[[ind2]]
      attr(Hl[[p + j]], "kernel") <- paste(attr(Hl[[ind1]], "kernel"),
                                           attr(Hl[[ind2]], "kernel"),
                                           sep = " x ")
    }
  }
  if (!is.null(no.int.3plus) & length(intr.3plus) > 0) {
    for (j in seq_len(no.int.3plus)) {
      lambda[p + no.int + j] <- Reduce("*", lambda[intr.3plus[,j]])
      Hl[[p + no.int + j]] <- Reduce("*", Hl[intr.3plus[, j]])
      intr.3plus.kernels <- sapply(Hl[intr.3plus[, j]], attr, "kernel")
      attr(Hl[[p + no.int + j]], "kernel") <- paste(intr.3plus.kernels,
                                                    collapse = " x ")
    }
  }

  if (is.null(env)) return(list(Hl = Hl, lambda = lambda))
  else {
    assign("Hl", Hl, envir = env)
    assign("lambda", lambda, envir = env)
  }
}

#' @export
.expand_Hl_and_lambda <- expand_Hl_and_lambda

BlockB_fn <- function(Hl, intr, n, p) {
  # This is the function which returns the BlockBStuff required for closed-form
  # EM algorithm.
  #
  # Args: The list of kernel matrices, (two-way) interaction specifications
  # intr, sample size n and number of variables p.
  #
  # Returns: BlockBStuff - a bunch of precalculated kernel matrices and a
  # function BlockB which tells iprior_em_closed() how to manipulate these
  # matrices.

  # Initialise -----------------------------------------------------------------
  Hl <- expand_Hl_and_lambda(Hl, seq_along(Hl), intr, NULL)$Hl
  environment(index_fn_B) <- environment()
  H2l <- Hsql <- Pl <- Psql <- Sl <- ind <- ind1 <- ind2 <- NULL
  BlockB <- function(k, x = lambda) NULL

  if (length(Hl) == 1L) {
    # CASE: Single lambda ------------------------------------------------------
    Pl <- Hl
    Psql <- list(fastSquare(Pl[[1]]))
    Sl <- list(matrix(0, nrow = n, ncol = n))
    BB.msg <- "Single lambda"
  } else {
    # Next, prepare the indices required for index_fn_B().
    z <- seq_along(Hl)
    ind1 <- rep(z, times = (length(z) - 1):0)
    ind2 <- unlist(lapply(2:length(z), function(x) c(NA, z)[-(0:x)]))
    # Prepare the cross-product terms of squared kernel matrices
    for (j in seq_along(ind1)) {
      H2l.tmp <- Hl[[ind1[j]]] %*% Hl[[ind2[j]]]
      H2l[[j]] <- H2l.tmp + t(H2l.tmp)
    }

    if (!is.null(intr)) {
      # CASE: Parsimonious interactions only ---------------------------------
      for (k in z) {
        Hsql[[k]] <- fastSquare(Hl[[k]])
        if (k <= p) ind[[k]] <- index_fn_B(k)  # only create indices for non-intr
      }
      BlockB <- function(k, x = lambda) {
        # Calculate Psql instead of directly P %*% P because this way
        # is < O(n^3).
        indB <- ind[[k]]
        lambda.P <- c(1, x[indB$k.int.lam])
        Pl[[k]] <<- Reduce("+", mapply("*", Hl[c(k, indB$k.int)], lambda.P,
                                       SIMPLIFY = FALSE))
        Psql[[k]] <<- Reduce("+", mapply("*", Hsql[indB$Psq],
                                         c(1, x[indB$Psq.lam] ^ 2),
                                         SIMPLIFY = FALSE))
        if (!is.null(indB$P2.lam1)) {
          lambda.P2 <- c(rep(1, sum(indB$P2.lam1 == 0)), x[indB$P2.lam1])
          lambda.P2 <- lambda.P2 * x[indB$P2.lam2]
          Psql[[k]] <<- Psql[[k]] +
            Reduce("+", mapply("*", H2l[indB$P2], lambda.P2, SIMPLIFY = FALSE))
        }
        lambda.PRU <- c(rep(1, sum(indB$PRU.lam1 == 0)), x[indB$PRU.lam1])
        lambda.PRU <- lambda.PRU * x[indB$PRU.lam2]
        Sl[[k]] <<- Reduce("+", mapply("*", H2l[indB$PRU], lambda.PRU,
                                       SIMPLIFY = FALSE))
      }
      BB.msg <- "Multiple lambda with parsimonious interactions"
    } else {
      # CASE: Multiple lambda with no interactions, or with non-parsimonious -
      # interactions ---------------------------------------------------------
      for (k in seq_along(Hl)) {
        Pl[[k]] <- Hl[[k]]
        Psql[[k]] <- fastSquare(Pl[[k]])
      }
      BlockB <- function(k, x = lambda) {
        ind <- which(ind1 == k | ind2 == k)
        Sl[[k]] <<- Reduce("+", mapply("*", H2l[ind], x[-k], SIMPLIFY = FALSE))
      }
      BB.msg <- "Multiple lambda with no interactions"
    }
  }

  list(H2l = H2l, Hsql = Hsql, Pl = Pl, Psql = Psql, Sl = Sl, ind1 = ind1,
       ind2 = ind2, ind = ind, BlockB = BlockB, BB.msg = BB.msg)
}

index_fn_B <- function(k) {
  # Indexer helper function used to create indices for H2l. This is a helper
  # function for BlockB_fn()
  #
  # Args: k is the index (out of p scale parameters) which is currently looped.
  # Note that intr, ind1 and ind2 are created in kernL().
  #
  # Returns: A list of various information which helps manipulates the kernel
  # matrices in BlockB_fn().
  ind.int1 <- intr[1, ] == k; ind.int2 <- intr[2, ] == k	# locating var/kernel matrix
  ind.int <- which(ind.int1 | ind.int2)  # of interactions (out of 1:no.int)
  k.int <- ind.int + p	# which kernel matrix has interactions involves k
  k.int.lam <- c(intr[1, ][ind.int2], intr[2, ][ind.int1])	# which has interaction with k?
  nok <- (1:p)[-k]	# all variables excluding k
  k.noint <- which(!(ind.int1 | ind.int2)) + p	# the opposite of k.int

  # P.mat %*% R.mat + R.mat %*% P.mat indices ----------------------------------
  grid.PR1 <- expand.grid(k, nok)
  za <- apply(grid.PR1, 1, findH2, ind1 = ind1, ind2 = ind2)
  grid.PR2 <- expand.grid(k.int, nok)
  zb <- apply(grid.PR2, 1, findH2, ind1 = ind1, ind2 = ind2)
  grid.PR.lam <- expand.grid(k.int.lam, nok)

  # P.mat %*% U.mat + U.mat %*% P.mat indices ----------------------------------
  grid.PU1 <- expand.grid(k, k.noint)
  zc <- apply(grid.PU1, 1, findH2, ind1 = ind1, ind2 = ind2)
  grid.PU2 <- expand.grid(k.int, k.noint)
  zd <- apply(grid.PU2, 1, findH2, ind1 = ind1, ind2 = ind2)
  grid.PU.lam <- expand.grid(k.int.lam, k.noint)

  # P.mat %*% P.mat indices ----------------------------------------------------
  grid.Psq <- t(combn(c(k, k.int), 2))
  ze <- apply(grid.Psq, 1, findH2, ind1 = ind1, ind2 = ind2)
  grid.Psq.lam <- NULL
  if (length(k.int.lam) > 0) grid.Psq.lam <- t(combn(c(0, k.int.lam), 2))

  list(
    k.int     = k.int,
    k.int.lam = k.int.lam,
    PRU       = c(za, zc, zb, zd),
    PRU.lam1  = c(rep(0, length(nok) + length(k.noint)),
                  grid.PR.lam[,1],
                  grid.PU.lam[,1]),
    PRU.lam2  = c(nok, k.noint, grid.PR.lam[,2], grid.PU.lam[,2]),
    Psq       = c(k, k.int),
    Psq.lam   = k.int.lam,
    P2        = ze,
    P2.lam1   = grid.Psq.lam[,1],
    P2.lam2   = grid.Psq.lam[,2]
  )
}

findH2 <- function(z, ind1, ind2){
  # This function finds position of H2 (cross-product terms of H). Used in
  # index_fn_B().
  #
  # Args: z is a dataframe created from expand.grid(), while ind1 and ind2 are
  # the indices of the cross-product matrices. These are created in BlockB_fn().
  #
  # Returns: A logical vector for use in index_fn_B().
  x <- z[1]; y <- z[2]
  which((ind1 == x & ind2 == y) | (ind2 == x & ind1 == y))
}

tab_intr_3plus <- function(x) {
  # Helper function which tabulates >2-way interactions.
  #
  # Args: x is a character vector of the form "1:2:3" etc. which basically says
  # that variable 1 interacts with variable 2 and 3.
  #
  # Returns: A matrix of interactions information.
  p <- max(sapply(strsplit(x, ":"), length))
  sapply(strsplit(x, ":"), function(x) {
    p_ <- length(x)
    if (p_ < p) as.numeric(c(x, rep(0, p - p_)))
    else as.numeric(x)
  })
}

where_int <- function(x) {
  # Function to determine where interactions are. Used in formula_to_xy().
  #
  # Args: x is a table obtained from the terms of a formula fit.
  #
  # Returns: The same table x but with replaced values.
  tmp <- x > 0
  x[tmp] <- which(x > 0)
  x[!tmp] <- 0
  x
}

which_intr_3plus <- function(x) {
  # Function to determine the positions of the >2-way interactions.
  #
  # Args: x is a character vector of the form "1:2:3" etc. which basically says
  # that variable 1 interacts with variable 2 and 3.
  #
  # Returns: A logical vector, with TRUE indicating that that position has a
  # >2-way interaction.
  sapply(strsplit(x, ""), function(x) length(x) > 3)
}

reorder_x <- function(x, smp) {
  # For the Nystrom method, the data (vector or matrix) needs to be reordered
  # according to the random sample.
  #
  # Args: x (vector or matrix), and smp is the sequence of reordering.
  #
  # Returns: A reordered vector or matrix x.
  if (is.null(nrow(x))) res <- x[smp]
  else res <- x[smp, ]
  mostattributes(res) <- attributes(x)
  res
}
