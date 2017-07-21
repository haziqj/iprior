fitted.ipriorMod_Nystrom <- function(object, ...) {
  this.env <- environment()
  list2env(object, this.env)
  list2env(Nystrom_eigen(ipriorKernel, lambda, psi), envir = this.env)
  y.hat <- alpha + psi * (V * rep(u ^ 2, each = nrow(V))) %*% (t(V) %*% a)
  # Effectively this is SR estimates...
  as.numeric(y.hat)[order(object$ipriorKernel$model$Nys.samp)]
}

predict.ipriorMod_Nystrom <- function(object, newdata, ...) {
  cat("write me")
}
