#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::List;
using Rcpp::Named;
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

//' Eigen decomposition of a matrix in C++
//'
//' Returns the eigenvalues and eigenvectors of a matrix \code{M}.
//'
//' @param M A symmetric, positive-definite matrix
//'
//'
//' @details A fast implementation of \code{eigen()} for symmetric,
//' positive-definite matrices. This helps speed up the I-prior EM algorithm.
//' @export
// [[Rcpp::export]]

Rcpp::List EigenCpp(Eigen::Map<Eigen::MatrixXd> M) {
    VectorXd values;
    MatrixXd vectors;
    SelfAdjointEigenSolver<MatrixXd> es(M);
    return List::create(
        Named("values") = es.eigenvalues(),
        Named("vectors") = es.eigenvectors()
    );
}
