#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::List;
using Rcpp::Named;
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

//' Eigen decomposition of a matrix in C++.
//'
//' Returns the eigenvalues and eigenvectors of a matrix X.
//'
//' A fast implementation of eigen for symmetric, positive-definite
//' matrices. This helps speed up the I-prior EM algorithm.
//'
//' @param X A symmetric, positive-definite matrix
//'
// [[Rcpp::export]]

Rcpp::List eigenCpp(Eigen::Map<Eigen::MatrixXd> X) {
    VectorXd values;
    MatrixXd vectors;
    SelfAdjointEigenSolver<MatrixXd> es(X);
    return List::create(
        Named("values") = es.eigenvalues(),
        Named("vectors") = es.eigenvectors()
    );
}
