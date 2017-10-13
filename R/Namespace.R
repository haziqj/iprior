# Rcpp and RcppEigen stuff -----------------------------------------------------
#' @useDynLib iprior, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
# iprior package imports -------------------------------------------------------
# @importFrom grDevices dev.new hcl
# @importFrom methods is
# @importFrom stats coef delete.response deviance dnorm kernel lm logLik
#   model.extract model.frame model.response optim optimHess pnorm printCoefmat
#   qnorm rnorm  terms fitted
# @importFrom utils combn setTxtProgressBar str txtProgressBar
#' @import ggplot2
NULL

# Hacky way to pass R CMD CHECK "no visible binding" note ----------------------
# globalVariables(c("BlockB", "BlockBstuff", "Hl", "Hlam.mat", "Pl", "Psql", "Sl",
#                   "V", "Var.Y.inv", "VarY.inv", "W.hat", "Y", "alpha",
#                   "force.nlm", "force.regEM", "hlamFn", "ind1", "ind2", "intr",
#                   "intr.3plus", "ipriorEM.env", "l", "lambda", "maxit", "model",
#                   "n", "nlm", "no.int", "no.int.3plus", "one.lam", "p", "parsm",
#                   "psi", "r", "report", "s", "stop.crit", "theta", "u", "w.hat",
#                   "x", "x0", "intercept", "probit", "rootkern", "Nystrom",
#                   "Nys.m", "Nys.seed", "Nys.kern", "Nys.samp", "A", "B", "a", "b",
#                   "y", "y.hat", "ipriorKernel"))
