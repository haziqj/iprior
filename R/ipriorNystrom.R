#' @export
ipriorNystrom <- function(...) {
  UseMethod("ipriorNystrom")
}

#' Nystrom method to fit I-prior models
#'
#' @param y Response variable
#'
#' @param ... Covariates
#' @param size Size of Nystrom sample
#' @param seed Choose a random seed, defaults to NULL
#' @param lambda.init NULL
#' @param psi.init NULL
#' @param model Model type
#' @param formula formula
#' @param data data
#'
#' @name ipriorNystrom
#' @export
ipriorNystrom.default <- function(y, ..., size, seed = NULL, lambda.init = NULL,
                                  psi.init = NULL, model = list()) {
  x.names <- as.character(as.list(match.call(expand.dots = FALSE))$...)

  # Take subsamples and prepare kernel -----------------------------------------
  if (is.ipriorKernel_Nystrom(y)) {
    ipriorKernel <- y
  } else {
    if (!is.null(seed)) set.seed(seed)
    Nys.samp <- sample(seq_along(y))
    model$Nys.kern <- size
    model$Nys.samp <- Nys.samp
    if (is.null(model$xname)) model$xname <- x.names
    model$kernel <- "FBM"
    ipriorKernel <- kernL(y, ..., model = model)
  }

  # Initialise parameters ------------------------------------------------------
  lambda <- lambda.init
  psi <- psi.init
  if (is.null(lambda)) lambda <- abs(rnorm(ipriorKernel$l, sd = 0.1))
  else {
    if (length(lambda) != ipriorKernel$l) {
      stop(paste("Incorrect dimension of lambda initial values. vector of
                 length", ipriorKernel$l, "required."), call. = FALSE)
    }
  }
  if (is.null(psi)) psi <- abs(rnorm(1))
  theta <- c(lambda, psi)

  # Minimise deviance = -2 * log-likelihood ------------------------------------
  res.optim <- optim(theta, deviance, object = ipriorKernel, method = "L-BFGS-B",
                     control = list(trace = 1), hessian = TRUE)

  # Change the call to "ipriorNystrom" -----------------------------------------------
  cl <- match.call()
  ynamefromcall <- as.character(cl[2])
  check.yname <- is.null(ipriorKernel$model$yname)
  if (check.yname) ipriorKernel$model$yname <- ynamefromcall
  cl[[1L]] <- as.name("ipriorNystrom")
  m <- match(c("control"), names(cl), 0L)
  if (any(m > 0)) cl <- cl[-m]

  # Results --------------------------------------------------------------------
  theta <- exp(res.optim$par)
  names(theta) <- c(paste0("lambda[", seq_len(ipriorKernel$l), "]"), "psi")
  lambda <- theta[-length(theta)]
  psi <- theta[length(theta)]
  alpha <- mean(ipriorKernel$Y); names(alpha) <- "alpha"
  coefficients <- c(alpha, theta)
  # tmp <- Nystrom_eigen(ipriorKernel, lambda, psi)
  # Vy.inv <- tmp$V %*% diag(1 / tmp$z) %*% t(tmp$V)
  # se.alpha <- 1 / sum(Vy.inv)
  se <- sqrt(c(0, diag(solve(res.optim$hessian + diag(1e-8, length(theta))))))
  loglik.adj <- logLik(ipriorKernel, res.optim$par, Nys.adj = TRUE)
  res <- list(alpha = alpha, lambda = lambda, psi = psi, coefficients = coefficients,
              se = se, ipriorKernel = ipriorKernel, loglik.adj = loglik.adj,
              loglik = -2 * res.optim$value, optim = res.optim,
              call = cl)
  class(res) <- "ipriorMod_Nystrom"
  res
}

#' @rdname ipriorNystrom
#' @export
ipriorNystrom.formula <- function(formula, data, size, seed = NULL, lambda.init = NULL,
                                  psi.init = NULL, model = list()) {
  # Take subsamples and prepare kernel -----------------------------------------
  if (!is.null(seed)) set.seed(seed)
  Nys.samp <- sample(seq_len(nrow(data)))
  model$Nys.kern <- size
  model$Nys.samp <- Nys.samp
  model$kernel <- "FBM"
  ipriorKernel <- kernL(formula, data, model = model)

  # Pass to default function ---------------------------------------------------
  res <- ipriorNystrom.default(y = ipriorKernel, size = size, seed = seed,
                               lambda.init = lambda.init, psi.init = psi.init)

  # Changing the call to simply "ipriorNystrom" --------------------------------
  cl <- match.call()
  cl[[1L]] <- as.name("ipriorNystrom")
  m <- match(c("formula", "data"), names(cl), 0L)
  cl <- cl[c(1L, m)]
  res$call <- cl
  names(res$call)[2] <- "formula"
  res$formula <- formula

  res
}
