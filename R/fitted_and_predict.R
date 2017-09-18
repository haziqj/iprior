predict_iprior <- function(y, Hlam, w) {
  y.hat <- Hlam %*% w
  resid <- y - y.hat
  train.error <- mean(resid ^ 2)
  list(y.hat = y.hat + attr(y, "scaled:center"), resid = resid, train.error = train.error)
}

se_yhat <- function(Hlam, psi) {
  list2env(eigen_Hlam(Hlam), environment())
  z <- psi * u ^ 2 + 1 / psi
  Vy.inv.Hlam <- vy_inv_a(1 / z, V, t(Hlam))
  diag(Hlam %*% Vy.inv.Hlam) + 1 / psi
}

fitted.ipriorMod2 <- function(object, intervals = FALSE, alpha = 0.05, ...) {
  y.hat <- object$fitted.values
  if (isTRUE(intervals)) {
    se <- se_yhat(get_Hlam(object$kernL, object$theta),
                  theta_to_psi(object$theta))
    lower <- y.hat + qnorm(alpha / 2) * se
    upper <- y.hat + qnorm(1 - alpha / 2) * se
    data.frame(y = y.hat, lower = lower, upper = upper)
  } else {
    return(y.hat)
  }
}
