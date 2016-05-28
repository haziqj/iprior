#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;

// [[Rcpp::export]]
Eigen::MatrixXd FastSquare(SEXP AA) {
    const Map<MatrixXd> S(as<Map<MatrixXd> >(AA));
    const int m(S.rows());
    const MatrixXd SS(MatrixXd(m,m).setZero().
                      selfadjointView<Lower>().rankUpdate(S.adjoint()));
    return SS;
}