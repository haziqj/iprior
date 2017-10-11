#' #' @export
#' print.ipriorMod_Nystrom <- function(x, ...) {
#'   class(x) <- "ipriorMod"
#'   print(x)
#' }
#'
#' #' @export
#' print.ipriorKernel_Nystrom <- function(x, ...) {
#'   class(x) <- "ipriorKernel"
#'   print(x)
#' }
#'
#' #' @export
#' summary.ipriorMod_Nystrom <- function(object, ...) {
#'   # Create table for summary output --------------------------------------------
#'   se <- object$se
#'   zval <- coef(object) / se
#'   tab <- cbind(Estimate   = round(coef(object), digits = 4),
#'                S.E.       = round(se, digits = 4),
#'                z          = round(zval, digits = 3),
#'                `P[|Z>z|]` = round(2 * pnorm(-abs(zval)), digits = 3))
#'   xname <- object$ipriorKernel$model$xname
#'   if (object$ipriorKernel$l == 1) {
#'     # only rename rows when using multiple lambdas
#'     lamnames <- c("(Intercept)", "lambda", "psi")
#'     rownames(tab) <- lamnames
#'   } else {
#'     lamnames <- paste0("lam", 1:(length(coef(object)) - 2))
#'     lamnames <- c("(Intercept)", paste(lamnames,
#'                                        object$ipriorKernel$model$lamnamesx,
#'                                        sep = "."), "psi")
#'     rownames(tab) <- lamnames
#'   }
#'   tab <- tab[-length(coef(object)), ]  # removes the psi from the table
#'
#'   res <- list(call = object$call, coefficients = tab,
#'               whichPearson = object$ipriorKernel$whichPearson,
#'               kernel = object$ipriorKernel$model$kernel,
#'               resid = object$ipriorKernel$Y[order(object$ipriorKernel$model$Nys.samp)] - fitted(object),
#'               log.lik = object$loglik,
#'               no.iter = object$no.iter, converged = FALSE,
#'               stop.crit = object$control$stop.crit,
#'               one.lam = object$ipriorKernel$model$one.lam, T2 = object$T2,
#'               l = object$ipriorKernel$l, p = object$ipriorKernel$p,
#'               Hurst = object$ipriorKernel$model$Hurst,formula = object$formula,
#'               psi.and.se = c(coef(object)[length(se)], se[length(se)]),
#'               xname = xname, no.int = object$ipriorKernel$no.int,
#'               optim.converged = object$optim.converged,
#'               rootkern = object$ipriorKernel$model$rootkern,
#'               Nystrom = list(m = object$ipriorKernel$model$Nys.kern),
#'               Nystrom.check = TRUE)
#'   class(res) <- "ipriorSummary"
#'   res
#' }
