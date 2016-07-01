lambdaExpand <- function(lambda_ = lambda, env = ipriorEM.env){
  assign("lambda", lambdaContract(lambda_, order), envir = env)
  if (parsm && no.int > 0) {
    for (j in 1:no.int) {
      assign("lambda", c(lambda, lambda_[intr[1, j]] * lambda_[intr[2, j]]),
             envir = env)
    }
  }
}

lambdaContract <- function(x, ord) {
  x[whereOrd(ord)]
}

ipriorEMClosedForm <- function() {
  # Update for lambda ----------------------------------------------------------
  BlockC()  # obtains Var.Y.inv and updates w.hat and W.hat
  for (k in 1:l) {
    BlockB(k)
    T1 <- sum(Psql[[k]] * W.hat)
    T2 <- 2 * crossprod(Y - alpha, crossprod(Pl[[k]], w.hat)) -
          sum(Sl[[k]] * W.hat)
    lambda[k] <<- as.vector(T2 / (2 * T1))
  }

  # Update for psi -------------------------------------------------------------
  Hlamsq.mat <<- fastSquare(Hlam.mat)  # a C++ alternative
  T3 <- crossprod(Y - alpha) + sum(Hlamsq.mat * W.hat) -
        2 * crossprod(Y - alpha, crossprod(Hlam.mat, w.hat))
  psi <<- sqrt(max(0, as.numeric(sum(diag(W.hat)) / T3)))

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
                     w.hat = w.hat, lambdaExpand = lambdaExpand,
                     hlamFn = hlamFn, env = ipriorEM.env)
  # print(theta.new$count)
  theta.new <- theta.new$par
  lambda <<- theta.new[-length(theta.new)]
  psi <<- theta.new[length(theta.new)]
}

QEstep <- function(theta, Y, alpha, W.hat, w.hat, lambdaExpand,
                   hlamFn, env) {
  Hlam.mat <- 0
  N <- length(Y)
  lambda <- theta[-length(theta)]
  psi <- theta[length(theta)]
  environment(lambdaExpand) <- environment(Hlam.mat.fn) <- env
  lambdaExpand(lambda_ = lambda, env = environment())
  Hlam.mat.fn(lambda_ = lambda, env = environment())
  Var.Y <- psi * fastSquare(Hlam.mat) + diag(1 / psi, N)
  Q <- psi * crossprod(Y - alpha) + sum(Var.Y * W.hat)
  Q <- Q - 2 * psi * crossprod(Y - alpha, Hlam.mat %*% w.hat)
  Q
}

###
### An internal function of ipriorEM() which does the pretty formatting in the report
###

decimalplaces <- function(x){	#this function calculates how many decimal places there are
  x <- as.numeric(x)
  if((x %% 1) != 0) nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  else return(0)
}

ipriorEMprettyIter <- function(x, iter){
  #first need to figure out the correct significant figures
  tmp <- prettyNum(x, digits=7, width=11)
  dig <- 7 + 10 - nchar(gsub(" ", "", tmp))
  for(i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11)


  #then trim down excess digits, because sometimes it overshoots to 11 width
  dig1 <- 7 + 10 - nchar(gsub(" ", "", tmp))
  ind <- dig1 < 7
  dig[ind] <- (dig + dig1 - 7)[ind]
  for(i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11)

  #sometimes it works, sometimes it doesn't. so run one more time.
  dig1 <- 7 + 10 - nchar(gsub(" ", "", tmp))
  ind <- dig1 < 7
  dig[ind] <- (dig + dig1 - 7)[ind]
  for(i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11)

  #after trimming to correct signif, just add the trailing zeroes.
  #makes use of function decimalplaces()
  miss <- 10 - nchar(gsub(" ", "", tmp))
  ind <- which(miss > 0)
  for(i in ind){
    if(!is.na(x[i])) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11, nsmall=ifelse(decimalplaces(tmp[i])==0, miss[i]-1, decimalplaces(tmp[i])+miss[i]))
  }

  #what remains must be scientific numbers
  miss <- 10 - nchar(gsub(" ", "", tmp))
  ind <- which(miss > 0)
  for(i in ind){
    if(!is.na(x[i])) tmp[i] <- formatC(as.numeric(tmp[i]), format="e", digits=decimalplaces(tmp[i])-4+miss[i], width=11)
  }

  #then print result
  Iter <- format(paste0("Iteration " , iter, ":"), width=16, just="left")
  cat(Iter, tmp, "")
}

ipriorEMprettyLoglik <- function(x){
  tmp <- prettyNum(x, digits=8, width=10)
  miss <- 10 - nchar(gsub(" ", "", tmp))
  tmp <- format(x, digits=8, width=10, nsmall=decimalplaces(tmp)+miss)
  #one more time
  #dig <- 10 - nchar(gsub(" ", "", tmp))#; print(dig); print(decimalplaces(tmp))
  #tmp <- format(tmp, digits=8, width=10, nsmall=5)
  tmp
}
