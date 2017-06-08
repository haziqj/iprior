////////////////////////////////////////////////////////////////////////////////
//
//   iprior: Linear Regression using I-priors
//   Copyright (C) 2017 Haziq Jamil
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

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP iprior_eigenCpp(SEXP);
extern SEXP iprior_fastSquare(SEXP);
extern SEXP iprior_fastSquareRoot(SEXP);
extern SEXP iprior_fastVDiag(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"iprior_eigenCpp",       (DL_FUNC) &iprior_eigenCpp,       1},
  {"iprior_fastSquare",     (DL_FUNC) &iprior_fastSquare,     1},
  {"iprior_fastSquareRoot", (DL_FUNC) &iprior_fastSquareRoot, 1},
  {"iprior_fastVDiag",      (DL_FUNC) &iprior_fastVDiag,      2},
  {NULL, NULL, 0}
};

void R_init_iprior(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
