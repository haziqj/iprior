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

#' Recover the betas (slopes) of the regression curves
#'
#' Since an I-prior model does not deal with the betas (slopes or coefficients
#' of an ordinary linear model), this function calculates the slopes of the
#' regression curves based on the fitted values. While it is meant for Canonical
#' RKHS functions, it still works for FBM RKHS, in that the slopes returned
#' would be the "average" slope, i.e. the slope of the fitted straight line
#' through the fitted I-prior values. Currently, only models with one
#' explanatory variable are supported.
#'
#'
#'
#' @param object Object of class "ipriorMod"
#'
#' @return The slopes of the fitted regression curves.
#'
#' @examples
#' mod.iprior <- iprior(stack.loss ~ Air.Flow, data = stackloss)
#' slope(mod.iprior)
#'
#' @export
slope <- function(object){
	if (!is.ipriorMod(object)) {
	  stop("Input objects of class ipriorMod only.", call. = FALSE)
	}

	y <- object$fitted.values
	x <- object$ipriorKernel$x
	whichPearson <- isPea(object$ipriorKernel$model$kernel)
	whichCan <- isCan(object$ipriorKernel$model$kernel)
	if (sum(whichPearson, whichCan) > 1) {
	  stop("Functionality not supported yet for more than 2 variables.",
	       call. = FALSE)
	}
	# if (any(isFBM(object$ipriorKernel$kernel))) {
	#   stop("Functionality is only meant for Canonical kernels.", call. = FALSE)
	# }
	dat <- cbind(y, as.data.frame(x))
	if (any(whichPearson)) {
		grp <- dat[[which(whichPearson) + 1]]
		dat <- split(dat, grp)
		res <- sapply(dat, function(x) coef(lm(x[, 1] ~ x[, 2]))[2])
		names(res) <- levels(grp)
	} else {
		res <- coef(lm(dat[, 1] ~ dat[, 2]))[2]
		names(res) <- "slope"
	}
	res
}
