#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;

// Multiplying a symmetric matrix by itself in C++.
//
// Returns the square of a symmetric matrix X.
//
// A fast implementation of X^2 for symmetric matrices. This helps
// speed up the I-prior EM algorithm.
//
// @param X A symmetric matrix
//
// [[Rcpp::export]]

Eigen::MatrixXd fastSquare(SEXP X) {
    const Map<MatrixXd> S(as<Map<MatrixXd> >(X));
    const int m(S.rows());
    const MatrixXd SS(MatrixXd(m,m).setZero().
                      selfadjointView<Lower>().rankUpdate(S.adjoint()));
    return SS;
}
