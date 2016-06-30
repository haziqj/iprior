#' Load the kernel matrices for \code{iprior} fit
#'
#' @param formula the model formula to fit
#' @param data data frame containing variables
#' @param model list of model options
#' @param control list of control options for EM algorithm and output
#'
#' @return an object of class iprior
#'
#' @examples
#' str(ToothGrowth)
#' mod <- kernL(y = ToothGrowth$len, supp = ToothGrowth$supp,
#'              dose = ToothGrowth$dose, model = list(interactions="1:2"))
#' summary(mod)
#' \dontrun{mod <- kernL(len ~ supp * dose, data = ToothGrowth)  # also accepts
#' formula}
#'
#' @export
kernL <- function(y, ..., model = list()) UseMethod("kernL")

#' @rdname kernL
#' @export
kernL.default <- function(y, ..., model = list()) {
  x <- list(...)
  if (any(sapply(x, is.list))) x <- unlist(x, recursive = FALSE)
  N <- length(y)
  p <- length(x)
  whichPearson <- unlist(lapply(x, is.factor))

  # Model options and checks ---------------------------------------------------
  mod <- list(kernel = "Canonical", Hurst = 0.5, interactions = NULL,
              parsm = TRUE, one.lam = FALSE, yname = "y", xname = NULL,
              silent = TRUE)
  mod_names <- names(mod)
  mod[(model_names <- names(model))] <- model
  if (length(noNms <- model_names[!model_names %in% mod_names])) {
    warning("Unknown names in model options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  mod$kernel <- match.arg(mod$kernel, c("Canonical", "FBM"))

  ## Set up interactions, p and q ----------------------------------------------
  names(mod)[3] <- "intr"  #rename to something simpler
  if (!is.null(mod$intr)) {
    if (!is.matrix(mod$intr)) {
      # not fitted using formula
      if (!is.character(mod$intr))
        stop("Incorrect prescription of interactions.")
      mod$intr <- sapply(strsplit(mod$intr, ":"), as.numeric)
    }
    no.int <- ncol(mod$intr)
    if (mod$parsm) q <- p
    else q <- p + no.int
  } else {
    # no interactions
    no.int <- 0L
    q <- p
  }
  if (mod$one.lam) {
    # only relevant when fitted using formula
    if (q == 1)
      message("Option one.lam = TRUE used with a single covariate anyway.")
    p <- q <- 1
  }
  if (any(mod$intr > p | mod$intr < 1)) {
    stop("Prescribed interactions out of bounds.")
  }

  # Set up names for x variables -----------------------------------------------
  if (is.null(mod$xname)) mod$xname <- names(x)
  else names(x) <- mod$xname[1:p]
  cond1 <- is.null(mod$xname)
  cond2 <- any(names(x) == "")
  cond3 <- any(is.na(names(x)))
  if (suppressWarnings(cond1 | cond2 | cond3)) {
    cl <- match.call()
    m <- match(c("y", "model", "control"), names(cl), 0L)
    xnamefromcall <- as.character(cl[-m])[-1]
    mod$xname <- xnamefromcall
  }
  suppressWarnings(here <- which((names(x) != "") & !is.na(names(x))))
  mod$xname[here] <- names(x)[here]
  names(x) <- mod$xname[1:p]

  # Set up list of H matrices --------------------------------------------------
  H.mat <- hMatList(x, mod$kernel, whichPearson, mod$intr, no.int, mod$Hurst)
  names(H.mat) <- mod$xname[1:length(H.mat)]
  if (length(mod$xname) < length(H.mat) && !mod$one.lam) {
    for (i in 1:ncol(mod$intr)) {
      mod$xname <- c(mod$xname, paste(mod$xname[mod$intr[1, i]],
                                      mod$xname[mod$intr[2, i]], sep = ":"))
    }
    names(H.mat) <- mod$xname
  }

  ### Set up progress bar ------------------------------------------------------
  if (!mod$silent)
    pb <- txtProgressBar(min = 0, max = 1, style = 3, width = 47)
  pb.count <- 0

  ### Block B update function --------------------------------------------------
  intr <- mod$intr
  environment(indx.fn) <- environment()
  H.mat2 <- H.matsq <- P.mat <- P.matsq <- S.mat <- ind <- ind1 <- ind2 <- NULL
  if (mod$one.lam | q == 1) {
    P.mat <- Reduce("+", mapply("*", H.mat, 1, SIMPLIFY = FALSE))
    P.matsq <- list(fastSquare(P.mat))
    if (mod$one.lam)
      mod$xname <- paste0("(", paste(names(H.mat), collapse = " + "), ")")
    H.mat <- P.mat <- list(P.mat)
    x <- list(as.data.frame(x))
    names(H.mat) <- names(x) <- mod$xname
    if (mod$one.lam)
      whichPearson <- FALSE  # never use Pearson kernel with one.lam
    BlockB <- function(k) NULL
    S.mat <- list(matrix(0, nrow = N, ncol = N))
    if (!mod$silent)
      setTxtProgressBar(pb, 1)
  } else {
    # Next, prepare the indices (also required for indx.fn).
    z <- 1:(p + no.int)
    ind1 <- rep(z, times = (length(z) - 1):0)
    ind2 <- unlist(lapply(2:length(z), function(x) c(NA, z)[-(0:x)]))
    if (!mod$silent) {
      pb <- txtProgressBar(min = 0, max = length(c(ind1, z)), style = 3,
                           width = 47)
    }
    # Prepare the cross-product terms of squared kernel matrices. This is a list
    # of (p+no.int)_choose_2.
    for (j in 1:length(ind1)) {
      H.mat2[[j]] <- H.mat[[ind1[j]]] %*% H.mat[[ind2[j]]] +
                     H.mat[[ind2[j]]] %*% H.mat[[ind1[j]]]
      pb.count <- pb.count + 1
      if (!mod$silent) setTxtProgressBar(pb, pb.count)
    }

    if (!is.null(intr) && mod$parsm) {
      # CASE: Parsimonious interactions only -----------------------------------
      for (k in z) {
        H.matsq[[k]] <- fastSquare(H.mat[[k]])
        if (k <= p)
          ind[[k]] <- indx.fn(k)
        pb.count <- pb.count + 1
        if (!mod$silent)
          setTxtProgressBar(pb, pb.count)
      }
      BlockB <- function(k) {
        indB <- ind[[k]]
        lambda.P <- c(1, lambda[indB$k.int.lam])
        P.mat[[k]] <<- Reduce("+", mapply("*", H.mat[c(k, indB$k.int)],
                                          lambda.P, SIMPLIFY = FALSE))
        P.matsq[[k]] <<- Reduce("+", mapply("*", H.matsq[indB$Psq],
                                            c(1, lambda[indB$Psq.lam] ^ 2),
                                            SIMPLIFY = FALSE))
        if (!is.null(indB$P2.lam1)) {
          lambda.P2 <- c(rep(1, sum(indB$P2.lam1 == 0)), lambda[indB$P2.lam1])
          lambda.P2 <- lambda.P2 * lambda[indB$P2.lam2]
          P.matsq[[k]] <<- P.matsq[[k]] +
                           Reduce("+", mapply("*", H.mat2[indB$P2], lambda.P2,
                                              SIMPLIFY = FALSE))
        }
        lambda.PRU <- c(rep(1, sum(indB$PRU.lam1 == 0)), lambda[indB$PRU.lam1])
        lambda.PRU <- lambda.PRU * lambda[indB$PRU.lam2]
        S.mat[[k]] <<- Reduce("+", mapply("*", H.mat2[indB$PRU], lambda.PRU,
                                          SIMPLIFY = FALSE))
      }
    } else {
      # CASE: Multiple lambda with no interactions, or with non-parsimonious ---
      # interactions -----------------------------------------------------------
      for (k in 1:q) {
        P.mat[[k]] <- H.mat[[k]]
        P.matsq[[k]] <- fastSquare(P.mat[[k]])
        pb.count <- pb.count + 1
        if (!mod$silent) setTxtProgressBar(pb, pb.count)
      }
      BlockB <- function(k) {
        ind <- which(ind1 == k | ind2 == k)
        S.mat[[k]] <<- Reduce("+", mapply("*", H.mat2[ind], lambda[-k],
                                          SIMPLIFY = FALSE))
      }
    }
  }
  if (!mod$silent) close(pb)
  mod <- mod[-8]  #remove silent control

  BlockBstuff <- list(H.mat2 = H.mat2, H.matsq = H.matsq, P.mat = P.mat,
                      P.matsq = P.matsq, S.mat = S.mat, ind1 = ind1,
                      ind2 = ind2, ind = ind, BlockB = BlockB)
  kernelLoaded <- list(Y = y, x = x, H.mat = H.mat, N = N, p = p, q = q,
                       no.int = no.int, whichPearson = whichPearson,
                       BlockBstuff = BlockBstuff, model = mod)
  class(kernelLoaded) <- "ipriorKernel"
  kernelLoaded
}

#' @export
kernL.formula <- function(formula, data, model = list(), ...) {
  mf <- model.frame(formula = formula, data = data)
  tt <- terms(mf)
  Terms <- delete.response(tt)
  x <- model.frame(Terms, mf)
  y <- model.response(mf)

  # For interactions -----------------------------------------------------------
  interactions <- NULL
  tmpo <- attr(tt, "order")
  if (any(tmpo > 2)) {
    stop("iprior does not currently work with higher order interactions.")
  }
  tmpf <- attr(tt, "factors")
  tmpf2 <- as.matrix(tmpf[-1, tmpo == 2])  #this obtains 2nd order interactions
  int2 <- apply(tmpf2, 2, function(x) which(x == 1))
  if (any(tmpo == 2)) interactions <- int2

  yname <- names(attr(tt, "dataClasses"))[1]
  xname <- attr(tt, "term.labels")
  kernelLoaded <- kernL(y = y, x, model = c(model,
                                            list(interactions = interactions,
                                                 yname = yname, xname = xname)))
  kernelLoaded
}

#' @export
print.ipriorKernel <- function(x, ...) {
  cat("\n")
  # if (x$model$kernel == 'Canonical') CanOrFBM <- 'Canonical' else CanOrFBM <-
  # paste0('Fractional Brownian Motion with Hurst coef. ', x$gamfbm) kerneltypes <-
  # c(CanOrFBM, 'Pearson', paste(CanOrFBM, '& Pearson')) if (all(x$whichPearson))
  # cat(kerneltypes[2], 'RKHS loaded') else { if (!all(x$whichPearson) &&
  # !any(x$whichPearson)) cat(kerneltypes[1], 'RKHS loaded') else
  # cat(kerneltypes[3], 'RKHS loaded') } if (x$q == 1 | x$model$one.lam) cat(',
  # with a single scale parameter.\n') else cat(', with', x$q, 'scale
  # parameters.\n')
  cat("Sample size = ", x$N, "\n")
  cat("Number of scale parameters, p = ", x$q, "\n")
  cat("Number of interactions = ", x$no.int, "\n")
  cat("\nInfo on H matrix:\n\n")
  str(x$H.mat)
  cat("\n")
}
