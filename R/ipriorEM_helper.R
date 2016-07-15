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

lambdaExpand <- function(x = lambda, env = ipriorEM.env){
  # Expands lambda from length l to correct size q = p + no.int, first by
  # expanding the higher order terms (if any), and then by adding the
  # interaction lambdas after that.
  lambda.tmp <- rep(NA, q)
  for (i in 1:q) {
    if (isHOrd(order[i])) {
      j.and.pow <- splitHOrd(order[i])
      j <- j.and.pow[1]
      pow <- j.and.pow[2]
      lambda.tmp[i] <- x[as.numeric(j)] ^ as.numeric(pow)
    }
    else lambda.tmp[i] <- x[as.numeric(order[i])]
  }
  if (parsm && no.int > 0) {
    for (j in 1:no.int) {
      add.lam <- lambda.tmp[intr[1, j]] * lambda.tmp[intr[2, j]]
      lambda.tmp[p + j] <- add.lam
    }
  }
  if (no.int.3plus > 0) {
    for (j in 1:no.int.3plus) {
      lambda.tmp[p + j + no.int] <- Reduce("*", lambda.tmp[intr.3plus[, j]])
    }
  }
  assign("lambda", lambda.tmp, envir = env)
}

lambdaContract <- function(x = lambda, env = ipriorEM.env) {
  # The opposite of lambdaExpand(). Looks for model$order vector and extracts
  # only the l lambdas.
  assign("lambda", x[whereOrd(order)], envir = env)
}

linSolvInv <- function(b = NULL){
  # Function to solve VarY %*% a = b,
  # where VarY = V %*% diag(1 / (u + s)) %*% V
  if (is.null(b)) a <- fastVDiag(V, 1 / (u + s))  # a C++ alternative
  else a <- V %*% (diag(1 / (u + s)) %*% (t(V) %*% b) )
  a
}

BlockA <- function(){
  # Block A update function
  lambdaExpand()
  hlamFn()
  A <- Hlam.mat
  assign("s", 1/psi, envir = parent.frame())
  tmp <- eigenCpp(A)  # a C++ alternative
  assign("u", psi * tmp$val ^ 2, envir = parent.frame())
  assign("V", tmp$vec, envir = parent.frame())
  assign("is.VarYneg", FALSE, envir = parent.frame())
  assign("is.VarYneg", any(u + s < 0), envir = parent.frame())
}

BlockC <- function(){
  # Block C update function
  assign("VarY.inv", linSolvInv(), envir = parent.frame())
  w.hat <- psi * Hlam.mat %*% (VarY.inv %*% matrix(Y - alpha, ncol = 1))
  W.hat <- VarY.inv + tcrossprod(w.hat)
  assign("w.hat", w.hat, envir = parent.frame())
  assign("W.hat", W.hat, envir = parent.frame())
}

logLikEM <- function(){
  # Function to calculate log-likelihood value within the EM routine
  a <- linSolvInv(Y - alpha)
  logdet <- Re(sum(log((u + s)[u + s > 0])))
  log.lik <- -(n / 2) * log(2 * pi) - logdet / 2 - crossprod(Y - alpha, a) / 2
  as.numeric(log.lik)
}

# alphaUpdate <- function() {
#   # DEPRECATED - MLE for alpha is actually mean(Y). This was used to be called
#   # in the EM routines below.
#   tmp.alpha <- crossprod(matrix(1, ncol = 1, nrow = n), Var.Y.inv)
#   alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))
#   assign("alpha", alpha, envir = parent.frame())
# }

psiUpdate <- function() {
  # The common psi update routine in the EM routine below.
  Hlamsq.mat <- fastSquare(Hlam.mat)  # a C++ alternative
  T3 <- crossprod(Y - alpha) + sum(Hlamsq.mat * W.hat) -
    2 * crossprod(Y - alpha, crossprod(Hlam.mat, w.hat))
  psi <- sqrt(max(0, as.numeric(sum(diag(W.hat)) / T3)))
  assign("psi", psi, envir = parent.frame())
}

ipriorEMClosedForm <- function() {
  # The main EM engine in most regular cases, i.e. when there are no higher
  # order terms present. It is fast because the parameters are obtained in
  # closed form.

  # Update for lambda ----------------------------------------------------------
  lambda.tmp <- rep(NA, l)

  for (k in 1:l) {
    lambdaExpand()
    BlockB(k)
    T1 <- sum(Psql[[k]] * W.hat)
    T2 <- 2 * crossprod(Y - alpha, crossprod(Pl[[k]], w.hat)) -
          sum(Sl[[k]] * W.hat)
    lambda.tmp[k] <- as.vector(T2 / (2 * T1))
  }
  assign("lambda", lambda.tmp, envir = parent.frame())

  # Update for psi -------------------------------------------------------------
  environment(psiUpdate) <- environment()
  psiUpdate()
  assign("psi", psi, envir = parent.frame())
}

# ipriorEMOptim <- function() {
#   BlockC()  # obtains VarY.inv and updates w.hat and W.hat
#   theta <- c(lambda, psi)
#   theta.new <- optim(theta, QEstep, method = "L-BFGS-B",
#                      lower = c(rep(-Inf, l), 1e-9),
#                      Y = Y, alpha = alpha, W.hat = W.hat, w.hat = w.hat,
#                      lambdaExpand = lambdaExpand, hlamFn = hlamFn,
#                      env = ipriorEM.env)$par
#   assign("lambda", theta.new[-length(theta.new)], envir = parent.frame())
#   assign("psi", theta.new[length(theta.new)], envir = parent.frame())
#
# }
#
# QEstep <- function(theta, Y, alpha, W.hat, w.hat, lambdaExpand, hlamFn, env) {
#   n <- length(Y)
#   lambda <- theta[-length(theta)]
#   psi <- theta[length(theta)]
#   environment(lambdaExpand) <- environment(hlamFn) <- env
#   lambdaExpand(lambda, env = environment())
#   hlamFn(lambda, env = environment())
#   Var.Y <- psi * fastSquare(Hlam.mat) + diag(1 / psi, n)
#   Q <- psi * crossprod(Y - alpha) + sum(Var.Y * W.hat)
#   Q <- Q - 2 * psi * crossprod(Y - alpha, (Hlam.mat %*% w.hat))
#   as.numeric(Q)
# }

ipriorEMOptim1 <- function() {
  # This is the EM engine when there are higher orders present, and only one
  # lambda present. The optim() routine uses method "Brent" and upper and lower
  # bounds.

  # Update for lambda ----------------------------------------------------------
  assign("lambda", optim(lambda, QEstepLambda, method = "Brent", lower = -1e9,
                         upper = 1e9, Y = Y, alpha = alpha, psi = psi,
                         W.hat = W.hat, w.hat = w.hat,
                         lambdaExpand = lambdaExpand, hlamFn = hlamFn,
                         env = ipriorEM.env, hessian = FALSE)$par,
         envir = ipriorEM.env)

  # Update for psi -------------------------------------------------------------
  environment(psiUpdate) <- environment()
  psiUpdate()
  assign("psi", psi, envir = parent.frame())
}

ipriorEMOptim2 <- function() {
  # This is the EM engine when there are higher orders present, and only one
  # lambda present. The optim() routine uses method "Nelder-Mead".

  # Update for lambda ----------------------------------------------------------
  assign("lambda", optim(lambda, QEstepLambda, Y = Y, alpha = alpha, psi = psi,
                         W.hat = W.hat, w.hat = w.hat,
                         lambdaExpand = lambdaExpand, hlamFn = hlamFn,
                         env = ipriorEM.env, hessian = FALSE)$par,
         envir = ipriorEM.env)

  # Update for psi -------------------------------------------------------------
  environment(psiUpdate) <- environment()
  psiUpdate()
  assign("psi", psi, envir = parent.frame())
}

QEstepLambda <- function(lambda, Y, alpha, psi, W.hat, w.hat, lambdaExpand,
                         hlamFn, env) {
  # The Q function for the E-step.
  n <- length(Y)
  environment(lambdaExpand) <- environment(hlamFn) <- env
  lambdaExpand(lambda, env = environment())
  hlamFn(lambda, env = environment())
  Var.Y <- psi * fastSquare(Hlam.mat) + diag(1 / psi, n)
  Q <- psi * crossprod(Y - alpha) + sum(Var.Y * W.hat)
  Q <- Q - 2 * psi * crossprod(Y - alpha, (Hlam.mat %*% w.hat))
  as.numeric(Q)
}

ipriorEMnlm <- function() {
  # An alternative minimiser for the E-step using nlm.

  # Update for lambda ----------------------------------------------------------
  assign("lambda", nlm(QEstepLambda, lambda, Y = Y, alpha = alpha, psi = psi,
                       W.hat = W.hat, w.hat = w.hat,
                       lambdaExpand = lambdaExpand, hlamFn = hlamFn,
                       env = ipriorEM.env, hessian = FALSE)$estimate,
         envir = ipriorEM.env)

  # Update for psi -------------------------------------------------------------
  environment(psiUpdate) <- environment()
  psiUpdate()
  assign("psi", psi, envir = parent.frame())
}

decimalPlaces <- function(x){
  # This function calculates how many decimal places there are.
  x <- as.numeric(x)
  if ((x %% 1) != 0) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  }
  else return(0)
}

ipriorEMprettyIter <- function(x, iter){
  # Used in ipriorEM.R to make the report of the EM prettier.

  # First need to figure out the correct significant figures -------------------
  tmp <- prettyNum(x, digits = 7, width = 11)
  dig <- 7 + 10 - nchar(gsub(" ", "", tmp))
  for (i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits = dig[i], width = 11)

  # Then trim down excess digits, because sometimes it overshoots to 11 width --
  # Sometimes it works, sometimes it doesn't - so run twice --------------------
  for (k in 1:2) {
    dig1 <- 7 + 10 - nchar(gsub(" ", "", tmp))
    ind <- dig1 < 7
    dig[ind] <- (dig + dig1 - 7)[ind]
    for (i in 1:length(dig)) {
      tmp[i] <- prettyNum(x[i], digits = dig[i], width = 11)
    }
  }

  # After trimming to correct significant figures, just add the trailing zeroes
  # Makes use of function decimalPlaces() --------------------------------------
  miss <- 10 - nchar(gsub(" ", "", tmp))
  ind <- which(miss > 0)
  for (i in ind) {
    if (!is.na(x[i])) {
      tmp[i] <- prettyNum(x[i], digits = dig[i], width = 11,
                          nsmall = ifelse(decimalPlaces(tmp[i]) == 0,
                                          miss[i] - 1,
                                          decimalPlaces(tmp[i]) + miss[i]))
    }
  }

  # What remains must be scientific numbers ------------------------------------
  miss <- 10 - nchar(gsub(" ", "", tmp))
  ind <- which(miss > 0)
  for(i in ind){
    if (!is.na(x[i])) {
      tmp[i] <- formatC(as.numeric(tmp[i]), format = "e",
                        digits = decimalPlaces(tmp[i]) - 4 + miss[i], width = 11)
    }
  }

  # Finally, print result ------------------------------------------------------
  Iter <- format(paste0("Iteration " , iter, ":"), width = 16, just = "left")
  cat(Iter, tmp, "")
}

ipriorEMprettyLoglik <- function(x){
  # Used in ipriorEM.R to make the log-likelihood reporting prettier.
  tmp <- prettyNum(x, digits = 8, width = 10)
  miss <- 10 - nchar(gsub(" ", "", tmp))
  tmp <- format(x, digits = 8, width = 10, nsmall = decimalPlaces(tmp) + miss)
  # One more time --------------------------------------------------------------
  # dig <- 10 - nchar(gsub(" ", "", tmp))
  # tmp <- format(tmp, digits = 8, width = 10, nsmall = 5)
  tmp
}
