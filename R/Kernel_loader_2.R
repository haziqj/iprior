kernL2 <- function(...) UseMethod("kernL2")

kernL2.default <- function(y, ..., kernel = "linear") {
  X <- list(...)

  # Meta -----------------------------------------------------------------------
  n <- length(y)
  p <- length(X)

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
  kernels[1:p] <- kernel
  # The next two lines ensure that the Pearson kernel is used for factors
  which.pearson <- unlist(lapply(X, function(x) {is.factor(x) | is.character(x)}))
  kernels <- correct_pearson_kernel(kernels, which.pearson)

  lambda <- rep(1, p)
  Hl <- mapply(kernel_translator, X, list(NULL), kernels, lambda, SIMPLIFY = FALSE)
  kernels <- get_kernels_from_Hl(Hl)


}

kernel_to_param <- function(kernels, lambda) {
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

  cbind(param, kernels)
}

param_to_theta <- function(x) {
  x <- x[, -ncol(x)]
  x$hurst <- qnorm(x$hurst)
  x$lengthscale <- log(x$lengthscale)
  x$offset <- log(x$offset)
  if (nrow(x) == 1) x$lambda <- log(x$lambda)
  res <- na.omit(c(as.matrix(x))) # don't forget psi
  list(theta = as.numeric(res), na = na.action(res))
}

theta_to_param <- function(theta, na.info, which.pearson) {
  # don't forget psi
  full.length <- length(c(theta, na.info))
  param <- matrix(NA, ncol = 5, nrow = full.length / 5)
  tmp <- c(param)
  tmp[-na.info] <- theta
  param[] <- tmp

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

param_translator <- function(x) {
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
  x[which.pearson] <- "pearson"
  x
}

kernel_translator <- function(x, y = NULL, kernel, lam.poly = 1) {
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
