################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2017  Haziq Jamil
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

#' Generate simulated data for smoothing models
#'
#' @param n Sample size.
#' @param seed (Optional) Random seed.
#'
#' @return A dataframe containing the response variable \code{y} and
#'   unidimensional explanatory variable \code{X}.
#'
#' @examples
#' gen_smooth(10)
#'
#' @export
gen_smooth <- function(n = 150, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  f <- function(x, truth = FALSE) {
    35 * dnorm(x, mean = 1, sd = 0.8) +
      65 * dnorm(x, mean = 4, sd = 1.5) +
      (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
      3 * dnorm(x, mean = 2.5, sd = 0.3)
  }
  x <- c(seq(0.2, 1.9, length = n * 5 / 8), seq(3.7, 4.6, length = n * 3 / 8))
  x <- sample(x, size = n)
  x <- x + rnorm(n, sd = 0.65)  # adding random fluctuation to the x
  x <- sort(x)
  y.err <- rt(n, df = 1)
  y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(n, mean = 4.1))  # adding random
  data.frame(y = y, X = x)
}

#' Random slopes model simulated data
#'
#' A simulated dataset to illustrate multilevel modelling using I-priors.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{\code{y}}{Response variable}
#'   \item{\code{x}}{Explanatory variable}
#'   \item{\code{grp}}{Factor, indicating the group of that observation.}
#' }
#'
#' @examples
#' data(simdat)
#' str(simdat)
"simdat"

#' High school and beyond dataset
#'
#' A national longitudinal survey of of students from public and private high
#' schools in the United States, with information such as students' cognitive
#' and non-cognitive skills, high school experiences, work experiences and
#' future plans collected.
#'
#' @format A data frame of 7185 observations on 3 variables. \describe{
#'   \item{\code{mathach}}{Math achievement.} \item{\code{ses}}{Socio-Economic
#'   status.} \item{\code{schoolid}}{Categorical variable indicating the school
#'   the student went to. Treated as \code{\link{factor}}.} }
#'
#' @source \href{http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/7896}{High
#'   School and Beyond, 1980: A Longitudinal Survey of Students in the United
#'   States (ICPSR 7896)}
#'
#' @references Rabe-Hesketh, S., & Skrondal, A. (2008). \emph{Multilevel and
#'   longitudinal modeling using Stata}. STATA press.
#' @references Raudenbush, S. W. (2004). \emph{HLM 6: Hierarchical linear and
#'   nonlinear modeling}. Scientific Software International.
#' @references Raudenbush, S. W., & Bryk, A. S. (2002). \emph{Hierarchical
#'   linear models: Applications and data analysis methods} (Vol. 1). Sage.
#'
#' @examples
#' data(hsb)
#' str(hsb)
"hsb"

#' High school and beyond dataset
#'
#' Smaller subset of \code{hsb}.
#'
#' A random subset of size 16 out of the original 160 groups.
#'
#' @format A data frame of 661 observations on 3 variables. \describe{
#'   \item{\code{mathach}}{Math achievement.} \item{\code{ses}}{Socio-Economic
#'   status.} \item{\code{schoolid}}{Categorical variable indicating the school
#'   the student went to. Treated as \code{factor}.} }
#'
#' @examples
#' data(hsbsmall)
#' str(hsbsmall)
"hsbsmall"

#' Air pollution and mortality
#'
#' Data on the relation between weather, socioeconomic, and air pollution
#' variables and mortality rates in 60 Standard Metropolitan Statistical Areas
#' (SMSAs) of the USA, for the years 1959-1961.
#'
#' Details.
#'
#' @format A data frame of 16 observations on 16 variables.
#' \describe{
#'   \item{\code{Mortality}}{Total age-adjusted mortality rate per 100,000.}
#'   \item{\code{Rain}}{Mean annual precipitation in inches.}
#'   \item{\code{Humid}}{Mean annual precipitation in inches.}
#'   \item{\code{JanTemp}}{Mean annual precipitation in inches.}
#'   \item{\code{JulTemp}}{Mean annual precipitation in inches.}
#'   \item{\code{Over65}}{Mean annual precipitation in inches.}
#'   \item{\code{Popn}}{Mean annual precipitation in inches.}
#'   \item{\code{Educ}}{Mean annual precipitation in inches.}
#'   \item{\code{Hous}}{Mean annual precipitation in inches.}
#'   \item{\code{Dens}}{Mean annual precipitation in inches.}
#'   \item{\code{NonW}}{Mean annual precipitation in inches.}
#'   \item{\code{WhiteCol}}{Mean annual precipitation in inches.}
#'   \item{\code{Poor}}{Mean annual precipitation in inches.}
#'   \item{\code{HC}}{Mean annual precipitation in inches.}
#'   \item{\code{NOx}}{Mean annual precipitation in inches.}
#'   \item{\code{SO2}}{Mean annual precipitation in inches.}
#' }
#' @references McDonald, G. C. and Schwing, R. C. (1973). Instabilities of
#'   regression estimates relating air pollution to mortality.
#'   \emph{Technometrics}, 15(3):463-481.
#'
#' @examples
#' data(pollution)
#' str(pollution)
"pollution"
