\name{predict.iprior}
\alias{predict.iprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.iprior(object, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
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
function (object, newdata = NULL, ...) 
{
    if (is.null(newdata)) 
        ystar <- fitted(object)
    else {
        if (!is.null(object$formula)) {
            xstar <- model.matrix(object$formula, newdata)
            xcolnames <- colnames(xstar)
            xrownames <- rownames(xstar)
            wheres.int <- (colnames(xstar) == "(Intercept)")
            xstar <- xstar[, !wheres.int]
            x <- as.matrix(xstar)
        }
        else {
            xstar <- as.matrix(newdata)
        }
        Y <- object$y
        X <- object$x
        p <- ncol(X)
        H.mat <- 0
        for (j in 1:p) {
            H.mat <- H.mat + fn.H2a(X[, j], xstar)
        }
        ystar <- as.vector(object$alpha + object$lambda * (H.mat \%*\% 
            object$w.hat))
        names(ystar) <- xrownames
    }
    ystar
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
