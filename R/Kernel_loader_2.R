kernL2 <- function(...) UseMethod("kernL2")

kernL2.default <- function(y, ..., kernel = "linear", interactions = NULL,
                           fixed.hyp = FALSE, nystrom = FALSE, nys.seed = NULL,
                           est.lambda = TRUE, est.hurst = FALSE,
                           est.lengthscale = FALSE, est.offset = FALSE,
                           est.psi = TRUE, lambda = 1, psi = 1) {
  Xl <- list(...)
  # It is common to make the mistake and type kernels instead of kernel. This
  # corrects it.
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

  # For Nystrom method: Reorder data and create Xl.Nys -------------------------
  if (as.numeric(nystrom) > 0 & as.numeric(nystrom) != n) {
    if (as.numeric(nystrom) == 1) nystrom <- floor(0.1 * n)
    if (!is.null(nys.seed)) set.seed(nys.seed)
    nys.samp <- sample(seq_along(y))
    y <- y[nys.samp]
    tmp <- lapply(Xl, reorder_x, smp = nys.samp)  # defined in .reorder_ipriorKernel()
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
    if (!isTRUE(formula.method)) intr.3plus <- add_zeroes_intr_3plus(intr.3plus)
    no.int.3plus <- ncol(intr.3plus)
  }

  if (isTRUE(nys.check)) {
    Hl <- get_Hl(Xl, Xl.nys, kernels, lambda)
  } else {
    Hl <- get_Hl(Xl, list(NULL), kernels, lambda)
  }
  kernels <- get_kernels_from_Hl(Hl)

  if (isTRUE(fixed.hyp)) {
    est.lambda <- est.hurst <- est.lengthscale <- est.offset <- est.psi <- FALSE
  }
  estl <- list(est.lambda = est.lambda, est.hurst = est.hurst,
               est.lengthscale = est.lengthscale, est.offset = est.offset,
               est.psi = est.psi)

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
  # |                 | EM.closed == TRUE |
  # |-----------------|------------------:|
  # | poly.deg        |                NA |
  # | est.lambda      |              TRUE |
  # | est.hurst       |             FALSE |
  # | est.lengthscale |             FALSE |
  # | est.offset      |             FALSE |
  # | est.psi         |              TRUE |
  # | nys.check        |             FALSE |
  BlockB.cond <- (
    all(is.na(poly.deg)) & !isTRUE(est.hurst) & !isTRUE(est.lengthscale) &
      !isTRUE(est.offset) & (isTRUE(est.lambda) | isTRUE(est.psi)) &
      !isTRUE(nys.check)
  )
  if (isTRUE(BlockB.cond)) {
    BlockBStuff <- BlockB_fn(Hl, intr, n, p)
  }

  res <- list(
    # Data
    y = y, Xl = Xl, Hl = Hl,
    # Model
    kernels = kernels, which.pearson = which.pearson, probit = probit,
    poly.deg = poly.deg, thetal = thetal, estl = estl,
    intr = intr, intr.3plus = intr.3plus, nystroml = nystroml,
    BlockBStuff = BlockBStuff,
    # Meta
    n = n, p = p, no.int = no.int, no.int.3plus = no.int.3plus,
    xname = xname, yname = yname, formula = NULL, terms = NULL
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

  # if (isTRUE(x$probit)) {
  #   cat("Categorical response variables\n")
  # } else if (is.ipriorKernel_nys(x)) {
  #   cat("Nystrom kernel approximation ()\n")
  # }

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
  if (x$thetal$n.theta > 0)
    cat(paste(names(x$thetal$theta), collapse = ", "))
  else
    cat("none")
}

print_kern <- function(x) {
  kern.type <- attr(x, "kernel")
  res <- capture.output(str(x))[1]
  res <- gsub(" num", kern.type, res)
  res
}

BlockB_fn <- function(Hl, intr, n, p) {
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
    # Next, prepare the indices required for indxFn().
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

# |                 | EM.closed == TRUE |
# |-----------------|------------------:|
# | poly.deg        |              NULL |
# | est.lambda      |              TRUE |
# | est.hurst       |             FALSE |
# | est.lengthscale |             FALSE |
# | est.offset      |             FALSE |
# | est.psi         |              TRUE |

get_Xl.nys <- function(object) {
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
