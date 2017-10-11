#' #' @export
#' fitted.ipriorMod_Nystrom <- function(object, ...) {
#'   this.env <- environment()
#'   list2env(object, this.env)
#'   list2env(Nystrom_eigen(ipriorKernel, lambda, psi), envir = this.env)
#'   y.hat <- alpha + psi * (V * rep(u ^ 2, each = nrow(V))) %*% (t(V) %*% a)
#'   # Effectively this is SR estimates...
#'   as.numeric(y.hat)[order(object$ipriorKernel$model$Nys.samp)]
#' }
#'
#' #' @export
#' predict.ipriorMod_Nystrom <- function(object, newdata = list(), ...) {
#'   list2env(object, environment())
#'   list2env(ipriorKernel, environment())
#'   list2env(model, environment())
#'   environment(.lambdaExpand) <- environment()
#'
#'   if (length(newdata) == 0) {
#'     return(cat("No new data supplied. Use fitted() instead."))
#'   } else {
#'     if (!is.null(object$formula)) {
#'       # Model has been fitted using formula interface
#'       mf <- model.frame(formula = object$formula, data = newdata)
#'       tt <- terms(mf)
#'       Terms <- delete.response(tt)
#'       xstar <- model.frame(Terms, newdata)
#'       xrownames <- rownames(xstar)
#'       if (one.lam) {
#'         xstar <- list(as.matrix(xstar))
#'       }
#'     } else {
#'       if (any(sapply(newdata, is.vector))) {
#'         newdata <- lapply(newdata, function(x) (as.matrix(x)))
#'       }
#'       xstar <- newdata
#'       xrownames <- rownames(do.call(cbind, newdata))
#'     }
#'
#'     # Define new kernel matrix -------------------------------------------------
#'     tmp <- .reorder_ipriorKernel(object$ipriorKernel, seq_len(Nys.kern))
#'     list2env(tmp, environment())
#'     Hl <- .hMatList(x, kernel, intr, no.int, model$Hurst, intr.3plus,
#'                     rootkern = FALSE, xstar)
#'     .lambdaExpand(object$lambda, env = environment())
#'     Hlam.new.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))
#'
#'     # Calculate fitted values --------------------------------------------------
#'     list2env(Nystrom_eigen(object$ipriorKernel, lambda, psi), environment())
#'     D <- psi * (V * rep(u, each = nrow(V))) %*% (t(V) %*% a)
#'     E <- Hlam.new.mat %*% solve(A, cbind(A, B) %*% D)
#'     ystar <- as.numeric(alpha + E)
#'     names(ystar) <- xrownames
#'   }
#'   ystar
#' }
