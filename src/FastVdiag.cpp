#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix FastVdiag(NumericMatrix X, NumericVector y) {
    unsigned int ncol = X.ncol();
    unsigned int nrow = X.nrow();
    int counter = 0;
    for (unsigned int j=0; j<ncol; j++) {
        for (unsigned int i=0; i<nrow; i++)  {
            X[counter++] *= y[j];
        }
    }
    return X;
}