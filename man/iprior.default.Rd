\name{iprior.default}
\alias{iprior.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
iprior.default(x, y, maxit = 50000, delt = 1e-05, report.int = 1000, silent = F, ...)
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
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
    silent = F, ...) 
{
    x <- as.matrix(x)
    y <- as.numeric(y)
    n <- length(y)
    est <- ipriorEM(x, y, maxit = maxit, delt = delt, report.int = report.int, 
        silent = silent)
    param <- c(est$alpha, est$lambda, est$psi)
    names(param) <- c("alpha", "lambda", "psi")
    Var.Y <- est$lambda * est$lambda * est$psi * est$H.matsq + 
        (1/est$psi) * diag(n)
    Var.Y.inv <- solve(Var.Y)
    w.hat <- (est$psi * est$lambda) * est$H.mat \%*\% (Var.Y.inv \%*\% 
        matrix(y - est$alpha, nc = 1))
    Var.w.hat <- Var.Y.inv + tcrossprod(w.hat)
    Y.hat <- est$alpha + est$lambda * as.vector(crossprod(est$H.mat, 
        w.hat))
    est$call <- match.call()
    est$coefficients <- param
    est$fitted.values <- Y.hat
    est$w.hat <- w.hat
    est$yval <- y
    est$xval <- x
    class(est) <- "iprior"
    est
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
