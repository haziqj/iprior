kernL2 <- function(...) UseMethod("kernL2")

kernL2.default <- function(y, ..., kernel = "linear", interactions = NULL,
                           fixed.hyp = FALSE,
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, lambda = 1, psi = 1) {
  Xl <- list(...)
  # It is common to make the mistake and type kernels instead of kernel. This
  # correct it.
  Xl.kernel.mistake <- match("kernels", names(Xl))
  if ("kernels" %in% names(Xl)) {
    kernel <- Xl[[Xl.kernel.mistake]]
    Xl[[Xl.kernel.mistake]] <- NULL
  }
  if (is.factor(y)) {
    probit <- TRUE
  } else {
    probit <- FALSE
    y <- scale(y, scale = FALSE)  # centre variables
  }

  # Meta -----------------------------------------------------------------------
  n <- length(y)
  p <- length(Xl)

  # What types of kernels? -----------------------------------------------------
  if (length(kernel) < p && length(kernel) > 1) {
    warning(paste0("Incomplete kernel specification (not of length ", p, ")"),
            call. = FALSE)
  }
  if (length(kernel) > p && length(kernel) > 1) {
    warning(paste0("Too many kernel options specification (not of length ", p, ")"),
            call. = FALSE)
  }
  kernels <- rep(NA, p)
  suppressWarnings(kernels[1:p] <- kernel)
  # The next two lines ensure that the Pearson kernel is used for factors
  which.pearson <- unlist(lapply(Xl, function(x) {is.factor(x) | is.character(x)}))
  kernels <- correct_pearson_kernel(kernels, which.pearson)

  Hl <- get_Hl(Xl, list(NULL), kernels, lambda)
  kernels <- get_kernels_from_Hl(Hl)

  if (isTRUE(fixed.hyp)) {
    est.lambda <- est.hurst <- est.lengthscale <- est.offset <- est.psi <- FALSE
  }
  est.list <- list(est.lambda = est.lambda, est.hurst = est.hurst,
                   est.lengthscale = est.lengthscale, est.offset = est.offset,
                   est.psi = est.psi)

  param <- kernel_to_param(kernels, lambda)
  tmp <- param_to_theta(param, est.list, log(psi))
  theta <- tmp$theta
  nt <- length(theta)
  param.na <- tmp$param.na
  theta.drop <- tmp$theta.drop
  theta.omitted <- tmp$theta.omitted
  poly.degree <- param$degree

  res <- list(
    y = y, Xl = Xl, n = n, p = p, nt = nt, kernels = kernels, Hl = Hl,
    which.pearson = which.pearson, param.na = param.na, probit = probit,
    poly.degree = poly.degree, theta = theta, theta.drop = theta.drop,
    theta.omitted = theta.omitted, est.list = est.list
  )
  class(res) <- "ipriorKernel2"
  res
}

print.ipriorKernel2 <- function(x) {
  cat("Sample size:", x$n, "\n")
  cat("No. of covariates:", length(x$Xl), "\n")
  cat("\n")
  cat("Kernel matrices:\n")
  for (i in seq_along(x$Hl)) {
    cat("", i, print_kern(x$Hl[[i]]), "\n")
  }
  cat("\n")
  cat("Hyperparameters to estimate:\n")
  if (x$nt > 0)
    cat(paste(names(x$theta), collapse = ", "))
  else
    cat("none")
}

print_kern <- function(x) {
  kern.type <- attr(x, "kernel")
  res <- capture.output(str(x))[1]
  res <- gsub(" num", kern.type, res)
  res
}

get_Hl <- function(Xl, yl = list(NULL), kernels, lambda) {
  # Args: List of data Xl and yl (optional), vector of same length of kernel
  # characters to instruct kernel_translator() which kernels to apply each of
  # the x and y. lambda is needed for the polynomial kernels.
  #
  # Output: List of kernel matrices.
  #
  # Notes: Except for polynomial kernels, these are not Hlam matrices.
  mapply(kernel_translator, Xl, yl, kernels, lambda, SIMPLIFY = FALSE)
}

kernel_to_param <- function(kernels, lambda) {
  # Args: kernels is a p-vector of kernels to apply on the data. lambda are the
  # scale parameters.
  #
  # Output: The param table.
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

param_to_theta <- function(param, est.list, logpsi = 0) {
  # Args: A param table, the list of parameters to be estimated and optional
  # logpsi value.
  #
  # Output: theta is a vector of parameters to be passed to optim or EM for
  # optimisation, including the logpsi value.
  #
  # Notes: theta is designed so that the values are unbounded, i.e. hurst is
  # Phi^{-1}(hurst), lengthscale is log(lengthscale), etc. theta_to_param()
  # reverses this for final presentation.
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

reduce_theta <- function(theta.full, est.list) {
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

expand_theta <- function(theta.reduced, theta.drop, theta.omitted) {
  theta.full <- theta.drop
  theta.full[!theta.drop] <- theta.reduced
  theta.full[theta.drop] <- theta.omitted
  theta.full
}

collapse_param <- function(param) {
  # Args: A param table.
  #
  # Output: A vectorised form of param with the na values removed. The param.na
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

theta_to_param <- function(theta, object) {
  # Args: A vector of parameters to be optimised, including logpsi. object must
  # be either a ipriorKernel2 type object, or a list containing param.na,
  # which.pearson and poly.degree.
  #
  # Output: A param table.
  #
  # Notes: The logpsi value is removed. To obtain this use theta_to_psi(). If
  # object is specified, then the param.na, which.pearson and poly.degree are
  # obtained from object.
  param.na <- object$param.na
  which.pearson <- object$which.pearson
  poly.degree <- object$poly.degree
  theta <- expand_theta(theta, object$theta.drop, object$theta.omitted)
  theta <- theta[-length(theta)]

  full.length <- length(c(theta, param.na))
  param <- matrix(NA, ncol = 4, nrow = full.length / 4)
  tmp <- c(param)
  tmp[-param.na] <- theta
  param[] <- tmp
  param <- cbind(param, degree = poly.degree)

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

theta_to_psi <- function(theta, object) {
  # Args: A vector of parameters to be optimised, including logpsi.
  #
  # Output: psi, the error precision.
  theta <- expand_theta(theta, object$theta.drop, object$theta.omitted)
  logpsi <- theta[length(theta)]
  exp(logpsi)
}

param_translator <- function(x) {
  # Args: Row vector from param table.
  #
  # Output: The kernel used.
  #
  # Notes: Used as a helper function in theta_to_param().
  hyperparam <- x[-1]
  if (!is.na(hyperparam[1]))
    return(paste0("fbm,", hyperparam[1]))
  if (!is.na(hyperparam[2]))
    return(paste0("se,", hyperparam[2]))
  if (!is.na(hyperparam[3]))
    return(paste0("poly", hyperparam[4], ",", hyperparam[3]))
  "linear"
}

correct_pearson_kernel <- function(x, which.pearson) {
  # Args: The kernel vector and which.pearson (logical), indicating which of the
  # x position uses the Pearson kernel.
  #
  # Output: The corrected kernel vector.
  #
  # Notes: When using theta_to_param(), unable to identify which data x uses the
  # Pearson kernel. This helper function corrects it by reading from the logical
  # which.pearson vector.
  x[which.pearson] <- "pearson"
  x
}

kernel_translator <- function(x, y = NULL, kernel, lam.poly = 1) {
  # Args: x, y (optional) data and kernel a character vector indicating which
  # kernel to apply x and y on. lam.poly is the scale for polynomial kernels.
  #
  # Output: A kernel matrix.
  #
  # Notes: Used as a helper function in get_Hl() to output list of kernel
  # matrices in kernL2() and predict(). For future expansion, add new kernels
  # here.
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

theta_to_collapsed_param <- function(theta, object) {
  # Args: theta (usually theta.full) and the ipriorKernel object.
  #
  # Output: A vector of parameters.
  #
  # Notes: This is a wrapper function for theta_to_param(). It is useful to get
  # the vector of parameters directly from theta.
  param <- theta_to_param(theta, object)
  c(collapse_param(param)$param, theta_to_psi(theta, object))
}
