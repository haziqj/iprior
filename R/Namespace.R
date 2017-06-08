# Rcpp and RcppEigen stuff -----------------------------------------------------
#' @useDynLib iprior, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
# iprior package imports -------------------------------------------------------
#' @importFrom grDevices dev.new hcl
#' @importFrom graphics abline lines par plot text
#' @importFrom methods is
#' @importFrom stats logLik optim coef delete.response kernel lm model.frame
#'   model.response pnorm printCoefmat qqnorm rnorm terms optimHess nlm
#' @importFrom utils combn setTxtProgressBar str txtProgressBar
#' @import RColorBrewer
NULL
