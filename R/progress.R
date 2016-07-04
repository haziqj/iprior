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
	if (!x$converged) warning("The EM has not yet converged.", call. = FALSE)
	rn <- rownames(x$res.loglik)

	# Parameters -----------------------------------------------------------------
	res <- x$res.loglik
	cn <- c("Log-likelihood", "Pred.log-lik.", "Delta(i,i-1)",
	        colnames(x$res.param)[-1])
	res <- cbind(res, x$res.param[, -1])
	res <- as.data.frame(res, row.names = rn)
	colnames(res) <- cn

	# Trim the table -------------------------------------------------------------
	no.iter <- x$no.iter
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
