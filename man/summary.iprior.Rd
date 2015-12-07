\name{summary.iprior}
\alias{summary.iprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.iprior(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
function (object, ...) 
{
    ipriorloglik <- function(theta) {
        alpha <- theta[1]
        lambda <- theta[-c(1, length(theta))]
        psi <- theta[length(theta)]
        n <- length(fitted(object))
        y <- object$yval
        Var.Y <- lambda * lambda * psi * object$H.matsq + (1/psi) * 
            diag(n)
        loglik <- dmvnorm(y - alpha, rep(0, n), Var.Y, log = T)
        loglik
    }
    FisherInformation <- -hessian(ipriorloglik, coef(object))
    InverseFisher <- solve(FisherInformation)
    se <- sqrt(diag(InverseFisher))
    zval <- coef(object)/se
    tab <- cbind(Estimate = coef(object), S.E. = se, z = zval, 
        `P[|Z>z|]` = 2 * pnorm(-abs(zval)))
    res <- list(call = object$call, coefficients = tab)
    class(res) <- "summary.iprior"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
