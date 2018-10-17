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

using Rcpp::as;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;

//' Multiplying a symmetric matrix by itself in C++.
//'
//' Returns the square of a symmetric matrix X.
//'
//' A fast implementation of X^2 for symmetric matrices. This helps
//' speed up the I-prior EM algorithm.
//'
//' @param X A symmetric matrix
//'
//' @export
//'
// [[Rcpp::export]]

Eigen::MatrixXd fastSquare(SEXP X) {
  const Map<MatrixXd> S(as<Map<MatrixXd> >(X));
  const int m(S.rows());
  const MatrixXd SS(MatrixXd(m,m).setZero().
                      selfadjointView<Lower>().rankUpdate(S.adjoint()));
  return SS;
}
