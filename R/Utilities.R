################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2016  Haziq Jamil
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

triangIndex <- function(k){
  # Function to list row and column index of upper triangular matrix including
  # diagonals.
  w <- 1:k
  cbind(
    row = rep(w, times = length(w):1 ) ,
    col = unlist(lapply(1:length(w), function(x) c(NA,w)[-(0:x)]))
  )
}

is.ipriorMod <- function(x) inherits(x, "ipriorMod")

is.ipriorKernel <- function(x) inherits(x, "ipriorKernel")

is.ipriorX <- function(x) inherits(x, "ipriorX")

testXForm <- function(x) {
  # Tests whether object x is a data frame fitted using formula interface.
  xform <- FALSE
  if (length(x) == 1) {
    if (is.data.frame(x[[1]])) {
      xform <- !is.null(attr(x[[1]], "terms"))
    }
  }
  xform
}

isHOrd <- function(x) {
  # Tests whether x contains ^ indicating higher order term.
  grepl("\\^", x)
}

whereOrd <- function(x) {
  # Index of non-higher order terms.
  grep("\\^", x, invert = TRUE)
}

lenHOrd <- function(x) {
  # How many higher order terms have been specified?
  length(grep("\\^", x))
}

splitHOrd <- function(x) {
  # Gets the level 1 index and the power it is raised to
  strsplit(x, "\\^")[[1]]
}

isCan <- function(x) x == "Canonical"

isPea <- function(x) x == "Pearson"

isFBM <- function(x) grepl("FBM", x)

canPeaFBM <- function(x, kernel, gamma, y) {
  if (isCan(kernel)) res <- fnH2(x, y)
  if (isPea(kernel)) res <- fnH1(x, y)
  if (isFBM(kernel)) res <- fnH3(x, y, gamma)
  res
}

whereInt <- function(x) {
  tmp <- x > 0
  x[tmp] <- which(x > 0)
  x[!tmp] <- 0
  x
}

whichIntr3Plus <- function(x) {
  sapply(strsplit(x, ""), function(x) length(x) > 3)
}

addZeroesIntr3Plus <- function(x) {
  p <- max(sapply(strsplit(x, ":"), length))
  sapply(strsplit(x, ":"), function(x) {
    p_ <- length(x)
    if (p_ < p) as.numeric(c(x, rep(0, p - p_)))
    else as.numeric(x)
  })
}

splitKernel <- function(kernel) {
  # Helper function to split the FBMs from the Hurst coefficients, if any
  paste(lapply(strsplit(kernel, ","), function(x) x[1]))
}

splitHurst <- function(kernel) {
  # Helper function to split the FBMs from the Hurst coefficients, if any
  suppressWarnings(
    tmp <- as.numeric(paste(lapply(strsplit(kernel, ","), function(x) x[2])))
  )
  # tmp[is.na(tmp)] <- 0.5
  tmp
}

hMatList <- function(x, kernel, intr, no.int, gamma, intr.3plus,
                     xstar = vector("list", p)) {
  # Helper function for creation of list of H matrices. Used in Kernel_loader.r
  # and predict.R
  p <- length(x)

  # Check how many Hurst coefficients are provided -----------------------------
  # if (length(gamma) < sum(isFBM(kernel))) {
  #   warning("Number of Hurst coefficients supplied is less than the number of FBM kernels used.", call. = FALSE)
  # }
  # if (length(gamma) > p) {
  #   stop("Number of Hurst coefficients supplied is more than the number of FBM kernels used.", call. = FALSE)
  # }

  suppressWarnings(
    H <- mapply(canPeaFBM, x = x, kernel = as.list(kernel),
                gamma = gamma, y = xstar, SIMPLIFY = FALSE)
  )
  if (!is.null(intr)) {
	  # Add in interactions, if any.
		for (j in 1:no.int) {
			H[[p + j]] <- H[[intr[1, j]]] * H[[intr[2, j]]]
			class(H[[p + j]]) <- paste(class(H[[intr[1,j]]]), class(H[[intr[2,j]]]),
			                           sep = " x ")
		}
  }
  if (!is.null(intr.3plus) & length(intr.3plus) > 0) {
    no.int.3plus <- ncol(intr.3plus)
    for (j in 1:no.int.3plus) {
      H[[p + j + no.int]] <- Reduce("*", H[intr.3plus[, j]])
      intr.3plus.class <- sapply(H[intr.3plus[, j]], class)
      class(H[[p + j + no.int]]) <- paste(intr.3plus.class, collapse = " x ")
    }
  }
	H
}

indxFn <- function(k) {
  # Indexer helper function used to create indices for H2l. Note: intr, ind1 and
  # ind2 are created in kernL().
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
  # indxFn(). z is a dataframe created from expand.grid().
  x <- z[1]; y <- z[2]
  which((ind1 == x & ind2 == y) | (ind2 == x & ind1 == y))
}

# flatten <- function(x) {
	# len <- sum(rapply(x, function(x) 1L))
	# y <- vector("list", len)
	# i <- 0L
	# rapply(x, function(x) { i <<- i+1L; y[[i]] <<- x })
	# y
# }

if (getRversion() < "3.3.0") {
  sigma <- function(object, ...) UseMethod("sigma")
}

#' Obtain the standard deviation of the residuals 'sigma'
#'
#' Extract the standard deviation of the residuals. For I-prior models, this is
#' \code{sigma = 1 / sqrt(psi)}.
#'
#' This basically obtains \code{object$sigma}. For \code{R (>= 3.3.0)} then
#' \code{sigma} is an S3 method with the default method coming from the
#' \code{stats} package.
#'
#' @param object An object of class \code{ipriorMod}.
#' @param ... This is not used here.
#'
#' @rawNamespace if (getRversion() >= "3.3.0") importFrom(stats,sigma)
#' @rawNamespace if (getRversion() < "3.3.0") export(sigma)
#' @name sigma
#' @export
sigma.ipriorMod <- function(object, ...) object$sigma


.onUnload <- function(libpath) {
  # Whenever you use C++ code in your package, you need to clean up after
  # yourself when your package is unloaded.
  library.dynam.unload("iprior", libpath)
}

#' Colour palette for \code{iprior} plots
#'
#' This is the colour palette used by the \code{iprior} package. It is based off
#' \code{RColorBrewer::brewer.pal}'s Set 1, Set 2 and Dark 2 palettes.
#'
#' @param x (optional) A vector of maximum length 16.
#'
#' @return The colour palette indexed by \code{x}.
#'
#' @export
ipriorColPal <- function(x) {
  colx <- c(RColorBrewer::brewer.pal(9, "Set1")[-9],
            RColorBrewer::brewer.pal(8, "Dark2"))
  colx[6] <- RColorBrewer::brewer.pal(8, "Set2")[6]
  colx[x]
}

# Hacky way to pass R CMD CHECK "no visible binding" note ----------------------
globalVariables(c("BlockB", "BlockBstuff", "Hl", "Hlam.mat", "Pl", "Psql", "Sl",
                  "V", "Var.Y.inv", "VarY.inv", "W.hat", "Y", "alpha",
                  "force.nlm", "force.regEM", "hlamFn", "ind1", "ind2", "intr",
                  "intr.3plus", "ipriorEM.env", "l", "lambda", "maxit", "model",
                  "n", "nlm", "no.int", "no.int.3plus", "one.lam", "p", "parsm",
                  "psi", "r", "report", "s", "stop.crit", "theta", "u", "w.hat",
                  "x", "x0", "intercept"))
