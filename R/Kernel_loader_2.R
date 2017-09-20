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
  Xl.formula <- match("Xl.formula", names(Xl))
  formula.method <- FALSE
  if ("Xl.formula" %in% names(Xl)) {
    Xl <- Xl[[Xl.formula]]
    formula.method <- TRUE
  }
  xname <- names(Xl)
  yname <- attr(y, "yname")
  if (is.factor(y)) {
    probit <- TRUE
  } else {
    probit <- FALSE
    y <- scale(y, scale = FALSE)  # centre variables
  }

  # Meta -----------------------------------------------------------------------
  n <- length(y)
  p <- length(Xl)
  if (is.null(xname)) xname <- paste0("X", seq_len(p))
  if (is.null(yname)) yname <- "y"

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
  if (!is.null(intr.3plus)) {
    if (!isTRUE(formula.method)) intr.3plus <- add_zeroes_intr_3plus(intr.3plus)
    no.int.3plus <- ncol(intr.3plus)
  }

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
    theta.omitted = theta.omitted, est.list = est.list, intr = intr,
    intr.3plus = intr.3plus, no.int = no.int, no.int.3plus = no.int.3plus,
    xname = xname, yname = yname, formula = NULL
  )

  class(res) <- "ipriorKernel2"
  res
}

kernL2.formula <- function(formula, data, kernel = "linear", ...) {
  mf <- model.frame(formula = formula, data = data)
  tt <- terms(mf)
  Terms <- delete.response(tt)
  x <- model.frame(Terms, mf)
  y <- model.response(mf)
  yname <- names(attr(tt, "dataClasses"))[1]
  xname <- names(x)
  x <- as.list(x)
  attr(x, "terms") <- NULL
  attr(y, "yname") <- yname

  # Interactions ---------------------------------------------------------------
  interactions <- NULL
  tmpo <- attr(tt, "order")
  tmpf <- attr(tt, "factors")
  tmpf2 <- as.matrix(tmpf[-1, tmpo == 2])  # this obtains 2nd order interactions
  int2 <- apply(tmpf2, 2, function(x) which(x == 1))
  if (any(tmpo == 2)) interactions <- int2
  intr.3plus <- NULL
  tmpf3 <- as.matrix(tmpf[-1, tmpo > 2])
  int3 <- apply(tmpf3, 2, whereInt)
  if (any(tmpo > 2)) intr.3plus <- int3
  interactions <- list(intr = interactions, intr.3plus = intr.3plus)

  res <- kernL2.default(y = y, Xl.formula = x, kernel = kernel,
                        interactions = interactions, ...)
  res$formula <- formula
  res$terms <- tt
  res
}

print.ipriorKernel2 <- function(x) {
  tmp <- expand_Hl_and_lambda(x$Hl, seq_along(x$Hl), x$intr, x$intr.3plus)

  cat("Sample size:", x$n, "\n")
  cat("No. of covariates:", length(x$Xl), "\n")
  cat("No. of interactions:", x$no.int + x$no.int.3plus, "\n")

  cat("\n")
  cat("Kernel matrices:\n")
  for (i in seq_along(tmp$Hl)) {
    cat("", i, print_kern(tmp$Hl[[i]]), "\n")
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
