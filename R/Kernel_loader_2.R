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
  kernels[which.pearson] <- "pearson"

  Hl <- mapply(kernel_translator, X, list(NULL), kernels, SIMPLIFY = FALSE)

  list(Hl)
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
