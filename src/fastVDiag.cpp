////////////////////////////////////////////////////////////////////////////////
//
//   iprior: Linear Regression using I-priors
//   Copyright (C) 2018 Haziq Jamil
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Lower;

//' Computing a quadratic matrix form in C++.
//'
//' Returns XdiagyXT.
//'
//' A fast implementation of XdiagyXT. This helps speed up
//' the I-prior EM algorithm.
//'
//' @param X A symmetric, square matrix of dimension \code{n} by \code{n}
//' @param y A vector of length \code{n}
//'
//' @export
//'
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
