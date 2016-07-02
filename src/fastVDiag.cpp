#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Lower;

// Computing a quadratic matrix form in C++.
//
// Returns XdiagyXT.
//
// A fast implementation of XdiagyXT. This helps speed up
// the I-prior EM algorithm.
//
// @param X A symmetric, square matrix of dimension \code{n} by \code{n}
// @param y A vector of length \code{n}
//
// [[Rcpp::export]]

NumericMatrix fastVDiag(NumericMatrix X, NumericVector y) {
    const Eigen::Map<Eigen::MatrixXd> XX(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X));
    Eigen::MatrixXd XT(XX.transpose());
    unsigned int ncol = X.ncol();
    unsigned int nrow = X.nrow();
    int counter = 0;
    for (unsigned int j=0; j<ncol; j++) {
        for (unsigned int i=0; i<nrow; i++)  {
            X[counter++] *= y[j];
        }
    }
    const Eigen::Map<Eigen::MatrixXd> X2(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X));
    const Eigen::MatrixXd S(X2*XT);
    return wrap(S);
}
