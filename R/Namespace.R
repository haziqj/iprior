# Rcpp and RcppEigen stuff -----------------------------------------------------
#' @useDynLib iprior, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
# iprior package imports -------------------------------------------------------
#' @importFrom stats coef delete.response dnorm fitted logLik model.extract
#'   model.frame model.response na.action na.omit optim optimHess pnorm predict
#'   printCoefmat qnorm resid residuals rnorm rt terms
#' @importFrom utils capture.output combn object.size setTxtProgressBar str
#'   txtProgressBar
#' @import ggplot2
NULL

# Hacky way to pass R CMD CHECK "no visible binding",note ----------------------
globalVariables(c("BlockBStuff", "estl", "ind1", "ind2", "interactions", "intr",
                  "intr.3plus", "ipriorEM.env", "Iteration", "kernels",
                  "lambda", "loglik", "lower", "m", "maxit", "mod", "n",
                  "niter", "no.int", "no.int.3plus", "p", "parsm", "Pl",
                  "probit", "psi", "Psql", "Sl","stop.crit", "tt", "u", "upper",
                  "V", "value", "Var1", "Var2", "variable", "Vy.inv.y", "X",
                  "Xl", "xrownames", "xstar", "y", "yname"))
