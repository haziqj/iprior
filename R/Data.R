################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2018  Haziq Jamil
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
#' @param x.jitter A small amount of jitter is added to the \code{X} variables
#'   generated from a normal distribution with mean zero and standard deviation
#'   equal to \code{x.jitter}.
#' @param xlim Limits of the \code{X} variables to generate from.
#'
#' @return A dataframe containing the response variable \code{y} and
#'   unidimensional explanatory variable \code{X}.
#'
#' @examples
#' gen_smooth(10)
#'
#' @export
gen_smooth <- function(n = 150, xlim = c(0.2, 4.6), x.jitter = 0.65,
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  f <- function(x, truth = FALSE) {
    35 * dnorm(x, mean = 1, sd = 0.8) +
      65 * dnorm(x, mean = 4, sd = 1.5) +
      (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
      3 * dnorm(x, mean = 2.5, sd = 0.3)
  }
  x <- c(seq(xlim[1], 1.9, length = n * 5 / 8),
         seq(3.7, xlim[2], length = n * 3 / 8))
  x <- sample(x, size = n)
  x <- x + rnorm(n, sd = x.jitter)  # adding random fluctuation to the x
  x <- sort(x)
  y.err <- rt(n, df = 1)
  y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(n, mean = 4.1))  # adding random
  data.frame(y = y, X = x)
}

#' Generate simulated data for multilevel models
#'
#' @param n Sample size. Input either a single number for a balanced data set,
#'   or a vector of length \code{m} indicating the sample size in each group.
#' @param seed (Optional) Random seed.
#' @param m Number of groups/levels.
#' @param sigma_e The standard deviation of the errors.
#' @param sigma_u0 The standard deviation of the random intercept.
#' @param sigma_u1 The standard deviation of the random slopes.
#' @param sigma_u01 The covariance of between the random intercept and the
#'   random slope.
#' @param beta0 The mean of the random intercept.
#' @param beta1 The mean of the random slope.
#' @param x.jitter A small amount of jitter is added to the \code{X} variables
#'   generated from a normal distribution with mean zero and standard deviation
#'   equal to \code{x.jitter}.
#'
#' @return A dataframe containing the response variable \code{y}, the
#'   unidimensional explanatory variables \code{X}, and the levels/groups
#'   (factors).
#'
#' @examples
#' gen_multilevel()
#'
#' @export
gen_multilevel <- function(n = 25, m = 6, sigma_e = 2, sigma_u0 = 2,
                           sigma_u1 = 2, sigma_u01 = -2, beta0 = 0, beta1 = 2,
                           x.jitter = 0.5, seed = NULL) {
  # Generates a data set according to the model \deqn{y_{ij} = \beta_{0j} +
  # \beta{1j}X_{ij} + \epsilon_{ij}} \deqn{\beta_{0j} \sim \text{N}(0,
  # \sigma_{u0}^2)} \deqn{\beta_{1j} \sim \text{N}(0, \sigma_{u1}^2)}
  # \deqn{\text{Cov}(\beta_{0j}, \beta_{1j}) = \sigma_{u01}} with
  # \eqn{i=1,\dots,n_j}  samples and \eqn{j=1,\dots,m} groups.
  if (!is.null(seed)) set.seed(seed)
  beta <- mvtnorm::rmvnorm(m, c(beta0, beta1),
                           sigma = matrix(c(sigma_u0 ^ 2, sigma_u01,
                                            sigma_u01, sigma_u1 ^ 2), nrow = 2))
  if (length(n) == 1) {
    n <- rep(n, m)
  } else if (length(n) != m) {
    stop("n must have length equal to the number of groups (m).", call. = FALSE)
  }
  dat <- as.data.frame(matrix(NA, nrow = sum(n), ncol = 3))
  n.cum <- c(0, cumsum(n))
  for (j in seq_len(m)) {
    x <- seq(0, 5, length = n[j]) + rnorm(n[j], sd = x.jitter)
    y <- beta[j, 1] + beta[j, 2] * x + rnorm(n[j], sd = sigma_e)
    dat[(n.cum[j] + 1):n.cum[j + 1], ] <- cbind(y, x, j)
  }
  names(dat) <- c("y", "X", "grp")
  dat$grp <- factor(dat$grp)
  dat
}

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

#' Results of I-prior cross-validation experiment on Tecator data set
#'
#' Results of I-prior cross-validation experiment on Tecator data set
#'
#' For the fBm and SE kernels, it seems numerical issues arise when using a
#' direct optimisation approach. Terminating the algorithm early (say using a
#' relaxed stopping criterion) seems to help.
#'
#' @format Results from iprior_cv cross validation experiment. This is a list of
#'   seven, with each component bearing the results for the linear, quadratic,
#'   cubic, fBm-0.5, fBm-MLE and SE I-prior models. The seventh is a summarised
#'   table of the results.
#'
#' @examples
#' # Results from the six experiments
#' print(tecator.cv[[1]], "RMSE")
#' print(tecator.cv[[2]], "RMSE")
#' print(tecator.cv[[3]], "RMSE")
#' print(tecator.cv[[4]], "RMSE")
#' print(tecator.cv[[5]], "RMSE")
#' print(tecator.cv[[6]], "RMSE")
#'
#' # Summary of results
#' print(tecator.cv[[7]])
#'
#' \dontrun{
#'
#' # Prepare data set
#' data(tecator, package = "caret")
#' endpoints <- as.data.frame(endpoints)
#' colnames(endpoints) <- c("water", "fat", "protein")
#' absorp <- -t(diff(t(absorp)))  # this takes first differences using diff()
#' fat <- endpoints$fat
#'
#' # Here is the code to replicate the results
#' mod1.cv <- iprior_cv(fat, absorp, folds = Inf)
#' mod2.cv <- iprior_cv(fat, absorp, folds = Inf, kernel = "poly2",
#'                      est.offset = TRUE)
#' mod3.cv <- iprior_cv(fat, absorp, folds = Inf, kernel = "poly3",
#'                      est.offset = TRUE)
#' mod4.cv <- iprior_cv(fat, absorp, method = "em", folds = Inf, kernel = "fbm",
#'                      control = list(stop.crit = 1e-2))
#' mod5.cv <- iprior_cv(fat, absorp, folds = Inf, kernel = "fbm",
#'                      est.hurst = TRUE, control = list(stop.crit = 1e-2))
#' mod6.cv <- iprior_cv(fat, absorp, folds = Inf, kernel = "se",
#'                      est.lengthscale = TRUE, control = list(stop.crit = 1e-2))
#'
#' tecator_res_cv <- function(mod) {
#'   res <- as.numeric(apply(mod$res[, -1], 2, mean))  # Calculate RMSE
#'   c("Training RMSE" = res[1], "Test RMSE" = res[2])
#' }
#'
#' tecator_tab_cv <- function() {
#'   tab <- t(sapply(list(mod1.cv, mod2.cv, mod3.cv, mod4.cv, mod5.cv, mod6.cv),
#'                   tecator_res_cv))
#'   rownames(tab) <- c("Linear", "Quadratic", "Cubic", "fBm-0.5", "fBm-MLE",
#'                      "SE-MLE")
#'   tab
#' }
#'
#' tecator.cv <- list(
#'   "linear"   = mod1.cv,
#'   "qudratic" = mod2.cv,
#'   "cubic"    = mod3.cv,
#'   "fbm-0.5"  = mod4.cv,
#'   "fbm-MLE"  = mod5.cv,
#'   "SE"       = mod6.cv,
#'   "summary"  = tecator_tab_cv()
#' )
#' }
#'
#'
"tecator.cv"
