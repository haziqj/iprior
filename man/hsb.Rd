% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\docType{data}
\name{hsb}
\alias{hsb}
\title{High school and beyond dataset}
\format{
A data frame of 7185 observations on 3 variables. \describe{
  \item{\code{mathach}}{Math achievement.} \item{\code{ses}}{Socio-Economic
  status.} \item{\code{schoolid}}{Categorical variable indicating the school
  the student went to. Treated as \code{\link{factor}}.} }
}
\source{
\href{http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/7896}{High
  School and Beyond, 1980: A Longitudinal Survey of Students in the United
  States (ICPSR 7896)}
}
\usage{
hsb
}
\description{
A national longitudinal survey of of students from public and private high
schools in the United States, with information such as students' cognitive
and non-cognitive skills, high school experiences, work experiences and
future plans collected.
}
\examples{
data(hsb)
str(hsb)
}
\references{
Rabe-Hesketh, S., & Skrondal, A. (2008). \emph{Multilevel and
  longitudinal modeling using Stata}. STATA press.

Raudenbush, S. W. (2004). \emph{HLM 6: Hierarchical linear and
  nonlinear modeling}. Scientific Software International.

Raudenbush, S. W., & Bryk, A. S. (2002). \emph{Hierarchical
  linear models: Applications and data analysis methods} (Vol. 1). Sage.
}
\keyword{datasets}
