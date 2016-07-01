ipriorEMClosedForm <- function() {
  # Update for lambda ----------------------------------------------------------
  BlockC()  # obtains Var.Y.inv and updates w.hat and W.hat
  for (k in 1:q) {
    BlockB(k)
    T1 <- sum(P.matsq[[k]] * W.hat)
    T2 <- 2*crossprod(Y - alpha, crossprod(P.mat[[k]], w.hat)) -
          sum(S.mat[[k]] * W.hat)
    lambda[k] <<- as.vector(T2/(2 * T1))
  }

  # Update for psi -------------------------------------------------------------
  H.mat.lamsq <<- fastSquare(H.mat.lam)  # a C++ alternative
  T3 <- crossprod(Y - alpha) + sum(H.mat.lamsq * W.hat) -
        2 * crossprod(Y - alpha, crossprod(H.mat.lam, w.hat))
  psi <<- sqrt(max(0, as.numeric(sum(diag(W.hat))/T3)))

  # Estimating alpha -----------------------------------------------------------
  # tmp.alpha <- crossprod(x0, Var.Y.inv)
  # alpha <- as.vector(tcrossprod(Y, tmp.alpha) / tcrossprod(x0, tmp.alpha))
}

ipriorEMOptim <- function() {
  BlockC()  # obtains Var.Y.inv and updates w.hat and W.hat
  theta <- c(lambda, psi)
  theta.new <- optim(theta, QEstep, method = "L-BFGS-B",
                     lower = c(rep(-Inf, length(theta) - 1), 1e-9),
                     Y = Y, alpha = alpha, W.hat = W.hat,
                     w.hat = w.hat, lambda.fn = lambda.fn,
                     H.mat.lam.fn = H.mat.lam.fn, env = ipriorEM.env)
  # print(theta.new$count)
  theta.new <- theta.new$par
  lambda <<- theta.new[-length(theta.new)]
  psi <<- theta.new[length(theta.new)]
}

QEstep <- function(theta, Y, alpha, W.hat, w.hat, lambda.fn,
                   H.mat.lam.fn, env) {
  H.mat.lam <- 0
  N <- length(Y)
  lambda <- theta[-length(theta)]
  psi <- theta[length(theta)]
  environment(lambda.fn) <- environment(H.mat.lam.fn) <- env
  lambda.fn(lambda_ = lambda, env = environment())
  H.mat.lam.fn(lambda_ = lambda, env = environment())
  Var.Y <- psi * fastSquare(H.mat.lam) + diag(1 / psi, N)
  Q <- psi * crossprod(Y - alpha) + sum(Var.Y * W.hat)
  Q <- Q - 2 * psi * crossprod(Y - alpha, H.mat.lam %*% w.hat)
  Q
}
