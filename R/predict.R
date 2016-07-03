#' Predict for iprior models.
#'
#' Predict for iprior models.
#'
#' Predict  for iprior models.
#'
#' @param object Objects of class \code{iprior}.
#' @param newdata (Optional) A data frame in which to look for variables with
#'   which to predict. If omitted, the fitted values are used.
#' @param ... This is not used here.
#'
#' @examples
#' \donttest{mod <- iprior(stack.loss ~ ., data=stackloss)}
#' \donttest{predict(mod, newdata = data.frame(c(1,2,3)))}
#'
#' @export
predict.ipriorMod <- function(object, newdata = list(), ...) {
  list2env(object$ipriorKernel, environment())
  list2env(model, environment())
  environment(lambdaExpand) <- environment()

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
      xstar <- newdata
      xrownames <- rownames(do.call(cbind, newdata))
    }

    # Define new kernel matrix -------------------------------------------------
    Hl <- hMatList(x, kernel, intr, no.int, model$Hurst, xstar)
    lambdaExpand(object$lambda, env = environment())
    Hlam.mat <- Reduce("+", mapply("*", Hl, lambda, SIMPLIFY = FALSE))

    # Calculate fitted values --------------------------------------------------
    ystar <- as.vector(object$alpha + (Hlam.mat %*% object$w.hat))
    names(ystar) <- xrownames
  }
  ystar
}
