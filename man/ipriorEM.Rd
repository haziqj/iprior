\name{ipriorEM}
\alias{ipriorEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ipriorEM(x, y, maxit = 50000, delt = 1e-05, report.int = 1000, silent = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{delt}{
%%     ~~Describe \code{delt} here~~
}
  \item{report.int}{
%%     ~~Describe \code{report.int} here~~
}
  \item{silent}{
%%     ~~Describe \code{silent} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, maxit = 50000, delt = 1e-05, report.int = 1000, 
    silent = F) 
{
    X <- x
    Y <- y
    N <- length(Y)
    p <- ncol(X)
    x0 <- rep(1, N)
    lambda <- abs(rnorm(1))
    alpha <- rnorm(1)
    psi <- abs(rnorm(1))
    H.mat <- 0
    for (j in 1:p) {
        H.mat <- H.mat + fn.H2a(X[, j])
    }
    H.matsq <- H.mat \%*\% H.mat
    Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * diag(N)
    Var.Y.inv <- chol2inv(chol(Var.Y))
    log.lik0 <- dmvnorm(Y - alpha, rep(0, N), Var.Y, log = T)
    if (!silent) 
        cat("START iter", 0, log.lik0, "\n")
    log.lik1 <- log.lik0 + 2 * delt
    i <- 0
    while ((i != maxit) && (abs(log.lik0 - log.lik1) > delt)) {
        i <- i + 1
        log.lik0 <- log.lik1
        tmp.alpha <- crossprod(x0, Var.Y.inv)
        alpha <- as.vector(tcrossprod(Y, tmp.alpha)/tcrossprod(x0, 
            tmp.alpha))
        w.hat <- psi * lambda * H.mat \%*\% (Var.Y.inv \%*\% matrix(Y - 
            alpha, nc = 1))
        Var.w.hat <- Var.Y.inv + w.hat \%*\% t(w.hat)
        T1 <- sum(H.matsq * Var.w.hat)
        T2 <- crossprod(Y - alpha, crossprod(H.mat, w.hat))
        lambda <- as.vector(T2/T1)
        Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * 
            diag(N)
        Var.Y.inv <- chol2inv(chol(Var.Y))
        w.hat <- psi * lambda * H.mat \%*\% (Var.Y.inv \%*\% matrix(Y - 
            alpha, nc = 1))
        Var.w.hat <- Var.Y.inv + w.hat \%*\% t(w.hat)
        T3 <- crossprod(Y - alpha) + lambda^2 * sum(H.matsq * 
            Var.w.hat) - 2 * lambda * crossprod(Y - alpha, crossprod(H.mat, 
            w.hat))
        psi <- as.vector(N/T3)
        Var.Y <- lambda * lambda * psi * H.matsq + (1/psi) * 
            diag(N)
        log.lik1 <- dmvnorm(Y - alpha, rep(0, N), Var.Y, log = T)
        check <- i\%\%report.int
        if (!is.na(check) && check == 0 && !silent) {
            if (log.lik1 >= log.lik0) 
                cat("INCREASE iter", i, log.lik1, "\n")
            else cat("DECREASE iter", i, log.lik1, "\n")
        }
    }
    if (!silent) 
        cat("EM complete.\n", "\nNumber of iterations =", i, 
            "\n")
    if (!silent) 
        cat("Log-likelihood = ", log.lik1, "\n")
    list(alpha = alpha, lambda = lambda, psi = psi, log.lik = log.lik1, 
        no.iter = i, H.mat = H.mat, H.matsq = H.matsq)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
