################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2017  Haziq Jamil
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

#' @export
print.ipriorKernel2 <- function(x, ...) {
  tmp <- expand_Hl_and_lambda(x$Hl, seq_along(x$Hl), x$intr, x$intr.3plus)

  # if (isTRUE(x$probit)) {
  #   cat("Categorical response variables\n")
  # } else if (is.ipriorKernel_nys(x)) {
  #   cat("Nystrom kernel approximation ()\n")
  # }

  cat("Sample size:", x$n, "\n")
  cat("No. of covariates:", length(x$Xl), "\n")
  cat("No. of interactions:", x$no.int + x$no.int.3plus, "\n")

  cat("\n")
  cat("Kernel matrices:\n")
  for (i in seq_along(tmp$Hl)) {
    cat("", i, print_kern(tmp$Hl[[i]], ...), "\n")
  }
  cat("\n")
  cat("Hyperparameters to estimate:\n")
  if (x$thetal$n.theta > 0)
    cat(paste(names(x$thetal$theta), collapse = ", "))
  else
    cat("none")
}

print_kern <- function(x, ...) {
  kern.type <- attr(x, "kernel")
  res <- capture.output(str(x, ...))[1]
  res <- gsub(" num", kern.type, res)
  res
}

#' @export
summary.ipriorKernel2 <- function(object, ...) {
  y <- object$y

  res <- list(call = object$call)
  class(res) <- "ipriorKernel_summary"
  res
}

print.ipriorKernel_summary <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

