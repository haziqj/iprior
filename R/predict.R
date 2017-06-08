################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2016  Haziq Jamil
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

#' Predict for I-prior models.
#'
#' Calculated predicted values of an I-prior model for a set of new data. If no
#' new data specified, then the fitted values are returned. When not using the
#' formula interface to fit the model, then the new data supplied for
#' \code{predict} must be coerced into a list.
#'
#' @param object Objects of class \code{ipriorMod}.
#' @param newdata (optional) A data frame in which to look for variables with
#'   which to predict. If omitted, the fitted values are used.
#'
#'   Note that when using non-formula to fit, the explanatory variables must be
#'   supplied in a list, such as \code{newdata = list(x1.new, x2.new, x3.new)}.
#' @param ... This is not used here.
#'
#' @examples
#' # Non-formula fit
#' \donttest{mod <- iprior(y = stack.loss,
#'                         air = stack.x[,1],
#'                         water = stack.x[,2],
#'                         acid = stack.x[,3])}
#' \donttest{predict(mod, newdata = list(air = 58, water = 20, acid = 87))}
#'
#' # Formula fit
#' \donttest{mod.orange <- iprior(circumference ~ . ^ 2, data = Orange[-1, ])}
#' \donttest{predict(mod.orange, Orange[1, ])}
#'
#' @name predict
#' @export
predict.ipriorMod <- function(object, newdata = list(), ...) {
  list2env(object$ipriorKernel, environment())
  list2env(model, environment())
  environment(.lambdaExpand) <- environment()

  if (length(newdata) == 0) {
    ystar <- object$fitted
  } else {
    if (!is.null(object$formula)) {
      # Model has been fitted using formula interface
      mf <- model.frame(formula = object$formula, data = newdata)
      tt <- terms(mf)
      Terms <- delete.response(tt)
      xstar <- model.frame(Terms, newdata)
      xrownames <- rownames(xstar)
      if (one.lam) {
        xstar <- list(as.matrix(xstar))
      }
    } else {
      if (any(sapply(newdata, is.vector))) {
        newdata <- lapply(newdata, function(x) t(as.matrix(x)))
      }
      xstar <- newdata
      xrownames <- rownames(do.call(cbind, newdata))
    }

    # Define new kernel matrix -------------------------------------------------
    Hl <- .hMatList(x, kernel, intr, no.int, model$Hurst, intr.3plus,
                    rootkern = FALSE, xstar)  # can't square root if
                                             # matrix not square
    .lambdaExpand(object$lambda, env = environment())
    if (rootkern) {
      Hlam.mat <- object$psi *
        Reduce("+", mapply("*", Hl, lambda ^ 2, SIMPLIFY = FALSE))
      w.hat <- varyinv(object) %*% (Y - object$alpha)
    } else {
      Hlam.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
      w.hat <- object$w.hat
    }

    # Calculate fitted values --------------------------------------------------
    ystar <- as.vector(object$alpha + (Hlam.mat %*% w.hat))
    names(ystar) <- xrownames
  }
  ystar
}
