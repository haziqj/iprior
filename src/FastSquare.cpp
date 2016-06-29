#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;

//' Multiplying a symmetric matrix by itself in C++
//'
//' Returns the square of a symmetric matrix \code{AA}.
//'
//' @param AA A symmetric matrix
//'
//'
//' @details A fast implementation of \code{AA %*% AA} for symmetric
//' matrices. This helps speed up the I-prior EM algorithm.
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd FastSquare(SEXP AA) {
    const Map<MatrixXd> S(as<Map<MatrixXd> >(AA));
    const int m(S.rows());
    const MatrixXd SS(MatrixXd(m,m).setZero().
                      selfadjointView<Lower>().rankUpdate(S.adjoint()));
    return SS;
}
