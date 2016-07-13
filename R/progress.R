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

#' EM algorithm progression results for fitted \code{ipriorMod} objects
#'
#' A table showing, for each EM iteration, the log-likelihood values, predicted
#' log-likelihood value, the change in log-likelihood value from the previous
#' iteration, and the progression, or trace, of the parameters. This table can
#' be called even if \code{silent = TRUE} or \code{progress = "none"} was called
#' when fitting the \code{ipriorMod} object.
#'
#' This is useful for diagnosing the EM algorithm, and for example, to obtain
#' the "traceplot" of the parameters. Note that the zeroth and final iterations
#' will always be shown in the table.
#'
#' @param object An object of class \code{ipriorMod}.
#' @param interval (optional) One of \code{"auto"}, \code{"all"}, or any number.
#'   This is an option to control how many EM iterations are displayed. Defaults
#'   to \code{"auto"}.
#'
#' @examples
#' mod <- iprior(len ~ . ^ 2, ToothGrowth, control = list(silent = TRUE))
#' progress(mod)
#'
#' # Works even when progress = "none" option called
#' mod <- iprior(len ~ . ^ 2, ToothGrowth, control = list(progress = "none"))
#' progress(mod, 50)
#' prog <- progress(mod, "all")
#' plot(prog$lambda1, type = "l")
#'
#' @export
progress <- function(object, interval = c("auto", "all", "input any number")) {
	if (class(object) != "ipriorMod") {
	  stop("Input ipriorMod class objectsonly.", call. = FALSE)
	}
  if (!is.null(object$optim.converged)) stop("I-prior model estimated using ipriorOptim - no progress report generated.", call. = FALSE)
	if (!object$converged) warning("The EM has not yet converged.", call. = FALSE)
  rn <- rownames(object$res.loglik)

	# Parameters -----------------------------------------------------------------
	res <- object$res.loglik
	cn <- c("Log-likelihood", "Pred.log-lik.", "Delta(i,i-1)",
	        colnames(object$res.param)[-1])
	res <- cbind(res, object$res.param[, -1])
	res <- as.data.frame(res, row.names = rn)
	colnames(res) <- cn

	# Trim the table -------------------------------------------------------------
	no.iter <- object$no.iter
	if (!is.numeric(interval)) interval <- match.arg(interval)
	if (interval == "auto") interval <- max(no.iter %/% 8, 1)
	if (interval == "all") interval <- 1
	if (interval == "input any number") {
	  stop("No, what I meant was that interval should be numeric! Try interval=10.",
	       call. = FALSE)
	}
	trim <- no.iter %/% interval
	first <- res[1, ]; last <- res[no.iter + 1, ]
	res.ind <- seq(from = interval + 1, to = (trim * interval) + 1 , by = interval)
	res <- res[res.ind, ]
	if (no.iter %% interval != 0) res <- rbind(first, res, last)
	else res <- rbind(first, res)
	res
}
