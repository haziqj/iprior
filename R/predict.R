#' Predict for iprior models.
#'
#' Predict for iprior models.
#'
#' Predict  for iprior models.
#'
#' @param object Objects of class \code{iprior}.
#' @param newdata (Optional) A data frame in which to look for variables with
#'   which to predict. If omitted, the fitted values are used.
#'
#' @examples
#' \donttest{mod <- iprior(stack.loss ~ ., data=stackloss)}
#' \donttest{predict(mod, newdata = data.frame(c(1,2,3)))}
#'
#' @export
predict.ipriorMod <- function(object, newdata = list(), ...) {
  list2env(object$ipriorKernel, environment())
  list2env(model, environment())

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
      xstar <- unlist(list(xstar), recursive = F)
    } else {
      xstar <- newdata
      xrownames <- rownames(do.call(cbind, newdata))
    }

    # Define new kernel matrix -------------------------------------------------
    H.mat <- hMatList(x, kernel, whichPearson, intr, no.int, gamma, xstar)
    lambda <- object$lambda
    if (parsm && no.int > 0) {
      for (j in 1:no.int) {
        lambda <- c(lambda, lambda[intr[1, j]] * lambda[intr[2, j]])
      }
    }
    H.mat.lam <- Reduce("+", mapply("*", H.mat[1:(p + no.int)],
                                    lambda[1:(p + no.int)], SIMPLIFY = F))

    # Calculate fitted values --------------------------------------------------
    ystar <- as.vector(object$alpha + (H.mat.lam %*% object$w.hat))
    names(ystar) <- xrownames
  }
  ystar
}
