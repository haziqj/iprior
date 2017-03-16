////////////////////////////////////////////////////////////////////////////////
//
//   iprior: Linear Regression using I-priors
//   Copyright (C) 2016  Haziq Jamil
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
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#include <RcppEigen.h>
// #include <unsupported/Eigen/MatrixFunctions>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// Obtaing matrix square root in C++.
//
// Returns the square root of a positive definite matrix X.
//
// An implementation of X^{1/2} for symmetric positive definite matrices.
//
// @param X A positive definite matrix
//
// [[Rcpp::export]]

Eigen::MatrixXd fastSquareRoot(SEXP X) {
  const Map<MatrixXd> A(as<Map<MatrixXd> >(X));
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
  Eigen::MatrixXd sqrtA = es.operatorSqrt();
  return sqrtA;
}
