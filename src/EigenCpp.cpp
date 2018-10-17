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

using Rcpp::List;
using Rcpp::Named;
using Eigen::Map;                     // 'maps' rather than copies
using Eigen::MatrixXd;                // variable size matrix, double precision
using Eigen::VectorXd;                // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;  // one of the eigenvalue solvers

//' Eigen decomposition of a matrix in C++.
//'
//' Returns the eigenvalues and eigenvectors of a matrix X.
//'
//' A fast implementation of eigen for symmetric, positive-definite
//' matrices. This helps speed up the I-prior EM algorithm.
//'
//' @param X A symmetric, positive-definite matrix
//'
//' @export
//'
// [[Rcpp::export]]

List eigenCpp(Eigen::Map<Eigen::MatrixXd> X) {
    VectorXd values;
    MatrixXd vectors;
    SelfAdjointEigenSolver<MatrixXd> es(X);
    return List::create(
        Named("values") = es.eigenvalues(),
        Named("vectors") = es.eigenvectors()
    );
}
