\name{iprior}
\alias{iprior}
\alias{iprior.default}
\alias{iprior.formula}
\alias{print.iprior}
\alias{predict.iprior}
\alias{summary.iprior}
\alias{print.summary.iprior}
\title{Linear regression using I-priors}
\description{Fit a linear regression model with I-priors.}
\usage{
iprior(x, y, ...)

\method{iprior}{default}(x, y, ...)\method{iprior}{formula}(formula, data = list(), ...)\method{print}{iprior}(x, ...)\method{summary}{iprior}(object, ...)\method{predict}{iprior}(object, newdata=NULL, ...)
}
\arguments{
  \item{x}{ a numeric design matrix for the model. }  \item{y}{ a numeric vector of responses. }  \item{formula}{ a symbolic description of the model to be fit. }  \item{data}{ an optional data frame containing the variables in the model. }  \item{object}{ an object of class \code{"iprior"}, i.e., a fitted model. }  \item{\dots}{ not used. }}
%\details{ %%  ~~ If necessary, more details than the description above ~~ }
\value{  An object of class \code{"iprior"}, basically a list including elements  \item{coefficients}{  a named vector of coefficients }  \item{fitted.values}{ fitted values }}
\references{Wicher Bergsma [2014]}
\author{Haziq Jamil}
\note{Will add further notes in the future...}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(cats, package="MASS")
mod1 <- iprior(Hwt~Bwt, data=cats)
mod1
summary(mod1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }