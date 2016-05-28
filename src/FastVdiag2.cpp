#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Lower;

// [[Rcpp::export]]
NumericMatrix FastVdiag2(NumericMatrix X, NumericVector y) {
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