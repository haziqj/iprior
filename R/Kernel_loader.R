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
  if (testXForm(x)) x <- unlist(x, recursive = FALSE)
  x <- lapply(x, as.matrix)
  N <- length(y)
  p <- length(x)
  whichPearson <- unlist(lapply(x, function(x) is.factor(x) | is.character(x) ))

  # Model options and checks ---------------------------------------------------
  mod <- list(kernel = "Canonical", Hurst = 0.5, interactions = NULL,
              parsm = TRUE, one.lam = FALSE, yname = "y", xname = NULL,
              silent = TRUE, order = as.character(1:p))
  mod_names <- names(mod)
  mod[(model_names <- names(model))] <- model
  if (length(noNms <- model_names[!model_names %in% mod_names])) {
    warning("Unknown names in model options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  mod$kernel <- match.arg(mod$kernel, c("Canonical", "FBM"))

  # Check for higher order terms -----------------------------------------------
  mod$order <- as.character(mod$order)
  hord.check <- all(mod$order == as.character(1:p))
  if (!hord.check) {
    hord.check1 <- length(mod$order) != p
    hord.check2 <- any(grepl("\\^", mod$order))
    if (hord.check1 | !hord.check2) {
      stop("Incorrect prescription of higher order terms.", call. = FALSE)
    }
  }
  mod$order <- sepOrd(mod$order)
  r <- length(mod$order$hord)

  # Set up interactions, p and q -----------------------------------------------
  names(mod)[3] <- "intr"  #rename to something simpler
  if (!is.null(mod$intr)) {
    # Interactions present
    if (!is.matrix(mod$intr)) {
      # Not fitted using formula
      intr.check1 <- is.character(mod$intr)
      intr.check2 <- all(grepl(":", mod$intr))
      if (!intr.check1 | !intr.check2) {
        stop("Incorrect prescription of interactions.", call. = FALSE)
      }
      mod$intr <- sapply(strsplit(mod$intr, ":"), as.numeric)
    }
    no.int <- ncol(mod$intr)
  } else {
    # No interactions
    no.int <- 0L
  }
  if (any(mod$intr > p | mod$intr < 1)) {
    stop("Prescribed interactions out of bounds.")
  }
  q <- p + no.int
  if (!mod$parsm) l <- q else l <- p - r
  # For clarity, the definitions of p, q, r, and l are
  # p = Number of H matrices in H.mat minus interactions = l + r
  # l = Number of unique lambdas (= q when parsm = TRUE)
  # r = Number of higher order terms
  # q = Length of expanded lambda = p + no.int
  # h = length(H.mat)

  # Set up names for x variables -----------------------------------------------
  if (is.null(mod$xname)) mod$xname <- names(x)
  else names(x) <- mod$xname[1:p]
  suppressWarnings(cond1 <- is.null(mod$xname))#; if (cond1) print("cond1")
  suppressWarnings(cond2 <- any(names(x) == ""))
  # if (!is.na(cond2)) {
  #   if (cond2) print("cond2")
  # }
  suppressWarnings(cond3 <- any(is.na(names(x))))#; if (cond3) print("cond3")
  if (cond1 | cond2 | cond3) {
    cl <- match.call()
    m <- match(c("y", "model", "control"), names(cl), 0L)
    xnamefromcall <- as.character(cl[-m])[-1]
    mod$xname <- xnamefromcall
  }
  suppressWarnings(here <- which((names(x) != "") & !is.na(names(x))))
  # print(here)
  mod$xname[here] <- names(x)[here]
  names(x) <- mod$xname[1:p]

  # Set up list of H matrices --------------------------------------------------
  Hl <- hMatList(x, mod$kernel, whichPearson, mod$intr, no.int, mod$Hurst)
  h <- length(Hl)
  names(Hl) <- mod$xname[1:h]
  if (length(mod$xname) < h && !mod$one.lam && !is.null(mod$intr)) {
    for (i in 1:no.int) {
      mod$xname <- c(mod$xname, paste(mod$xname[mod$intr[1, i]],
                                      mod$xname[mod$intr[2, i]], sep = ":"))
    }
    names(Hl) <- mod$xname
  }

  ### Set up progress bar ------------------------------------------------------
  if (!mod$silent) pb <- txtProgressBar(min = 0, max = 1, style = 3, width = 47)
  pb.count <- 0

  ### Block B update function --------------------------------------------------
  intr <- mod$intr
  environment(indxFn) <- environment()
  H2l <- Hsql <- Pl <- Psql <- Sl <- ind <- ind1 <- ind2 <- NULL
  if (q == 1) {
    Pl <- Hl
    Psql <- list(fastSquare(Pl[[1]]))
    BlockB <- function(k) NULL
    Sl <- list(matrix(0, nrow = N, ncol = N))
    if (!mod$silent) setTxtProgressBar(pb, 1)
  } else {
    # Next, prepare the indices required for indxFn().
    z <- 1:h
    ind1 <- rep(z, times = (length(z) - 1):0)
    ind2 <- unlist(lapply(2:length(z), function(x) c(NA, z)[-(0:x)]))
    if (!mod$silent) {
      pb <- txtProgressBar(min = 0, max = length(c(ind1, z)), style = 3,
                           width = 47)
    }
    # Prepare the cross-product terms of squared kernel matrices. This is a list
    # of q_choose_2.
    for (j in 1:length(ind1)) {
      H2l[[j]] <- Hl[[ind1[j]]] %*% Hl[[ind2[j]]] +
                  Hl[[ind2[j]]] %*% Hl[[ind1[j]]]
      pb.count <- pb.count + 1
      if (!mod$silent) setTxtProgressBar(pb, pb.count)
    }

    if (!is.null(intr) && mod$parsm) {
      # CASE: Parsimonious interactions only -----------------------------------
      for (k in z) {
        Hsql[[k]] <- fastSquare(Hl[[k]])
        if (k <= p) ind[[k]] <- indxFn(k)  # only create indices for non-intr
        pb.count <- pb.count + 1
        if (!mod$silent) setTxtProgressBar(pb, pb.count)
      }
      BlockB <- function(k) {
        indB <- ind[[k]]
        lambda.P <- c(1, lambda[indB$k.int.lam])
        Pl[[k]] <<- Reduce("+", mapply("*", Hl[c(k, indB$k.int)], lambda.P,
                                       SIMPLIFY = FALSE))
        Psql[[k]] <<- Reduce("+", mapply("*", Hsql[indB$Psq],
                                         c(1, lambda[indB$Psq.lam] ^ 2),
                                         SIMPLIFY = FALSE))
        if (!is.null(indB$P2.lam1)) {
          lambda.P2 <- c(rep(1, sum(indB$P2.lam1 == 0)), lambda[indB$P2.lam1])
          lambda.P2 <- lambda.P2 * lambda[indB$P2.lam2]
          Psql[[k]] <<- Psql[[k]] + Reduce("+", mapply("*", H2l[indB$P2],
                                                       lambda.P2,
                                                       SIMPLIFY = FALSE))
        }
        lambda.PRU <- c(rep(1, sum(indB$PRU.lam1 == 0)), lambda[indB$PRU.lam1])
        lambda.PRU <- lambda.PRU * lambda[indB$PRU.lam2]
        Sl[[k]] <<- Reduce("+", mapply("*", H2l[indB$PRU], lambda.PRU,
                                       SIMPLIFY = FALSE))
      }
    } else {
      # CASE: Multiple lambda with no interactions, or with non-parsimonious ---
      # interactions -----------------------------------------------------------
      for (k in 1:q) {
        Pl[[k]] <- Hl[[k]]
        Psql[[k]] <- fastSquare(Pl[[k]])
        pb.count <- pb.count + 1
        if (!mod$silent) setTxtProgressBar(pb, pb.count)
      }
      BlockB <- function(k) {
        ind <- which(ind1 == k | ind2 == k)
        Sl[[k]] <<- Reduce("+", mapply("*", H2l[ind], lambda[-k],
                                       SIMPLIFY = FALSE))
      }
    }
  }
  if (!mod$silent) close(pb)
  mod <- mod[-8]  #remove silent control

  BlockBstuff <- list(H2l = H2l, Hsql = Hsql, Pl = Pl, Psql = Psql, Sl = Sl,
                      ind1 = ind1, ind2 = ind2, ind = ind, BlockB = BlockB)
  kernelLoaded <- list(Y = y, x = x, Hl = Hl, N = N, p = p, l = l, r = r,
                       no.int = no.int, q = q, whichPearson = whichPearson,
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
  yname <- names(attr(tt, "dataClasses"))[1]
  # xname <- attr(tt, "term.labels")
  xname <- names(x)
  xnl <- length(xname)

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

  # Deal with one.lam option ---------------------------------------------------
  one.lam <- FALSE
  if (any(names(model) == "one.lam")) one.lam <- model$one.lam
  if (one.lam) {
    if (!is.null(interactions)) {
      stop("Cannot use option one.lam = TRUE with interactions.", call. = FALSE)
    }
    if (ncol(x) == 1) {
      message("Option one.lam = TRUE used with a single covariate anyway.")
    }
    attributes(x)$terms <- attributes(x)$names <- NULL
    if (xnl <= 3) {
      xname <- paste(xname, collapse = " + ")
    } else {
      xname <- paste(xname[1], "+ ... +", xname[xnl])
    }
    x <- as.data.frame(x)
  }

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
  cat("Number of x variables, p = ", x$p, "\n")
  cat("Number of scale parameters, l = ", x$l, "\n")
  cat("Number of interactions = ", x$no.int, "\n")
  cat("\nInfo on H matrix:\n\n")
  str(x$Hl)
  cat("\n")
}
