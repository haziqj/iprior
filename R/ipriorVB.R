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
#
# A Variational Bayes implementation of iprior. Not currently used.
#
# ## ---- prelim ----
# library(ggplot2)
# library(progress)
# library(gridExtra)
#
# ## ---- variational.bayes ----
# ipriorVB <- function(y, X, kernel = "Canonical", maxit = 1000,
#                      stop.crit = 1e-7, silent = FALSE) {
#   if (kernel == "FBM") H <- iprior::fnH3(X)
#   else H <- iprior::fnH2(X)
#   H2 <- H %*% H
#   if (!silent) pb <- txtProgressBar(min = 0, max = maxit - 1, style = 3)
#   n <- length(y)
#
#   # Set up parameter results
#   lower.bound <- xi <- psi <- rep(NA, maxit)
#   u <- matrix(NA, ncol = n, nrow = maxit)
#
#   # Initialise
#   alpha <- mean(y)
#   xi[1] <- 1
#   xi2 <- 1
#   psi[1] <- 0.1
#   u[1, ] <- rep(0, n)
#   niter <- 1
#
#   for (t in 1:(maxit - 1)) {
#     # Update u
#     A <- xi2 * H2 + diag(1, n)
#     a <- as.numeric(xi[t] * crossprod(H, y - alpha))
#     tmp <- eigenCpp(A)
#     V <- tmp$vec
#     v <- abs(tmp$val)
#     u[t + 1, ] <- V %*% (diag(1 / v) %*% (t(V) %*% a) )
#     u.var <- fastVDiag(tmp$vec, 1 / tmp$val)
#     U <- u.var / psi[t] + tcrossprod(u[t + 1, ])
#
#     # Update lambda
#     ct <- sum(H2 * U); print(ct)
#     d <- as.numeric(crossprod(y - alpha, H) %*% u[t + 1, ])
#     xi[t + 1] <- d / ct
#     xi2 <- 1 / (psi[t] * ct) + (d / ct) ^ 2
#
#     # Update psi
#     A <- xi2 * H2 + diag(1, n)
#     r <- (sum((y - alpha) ^ 2) +  sum(A * U) - 2 * xi[t + 1] * d) / 2
#     psi[t + 1] <- (n + 1) / r
#
#     # Lower bound
#     lower.bound[t + 1] <- ((n + 1) / 2) * (1 - log(r) - log(n + 1)) -
#       ((n - 1) / 2) * log(2 * pi) + lgamma(n + 1) -
#       (1 / 2) * (determinant(A)$mod + log(ct))
#
#     lb.diff <- abs(lower.bound[t + 1] - lower.bound[t])
#     if (!is.na(lb.diff) && (lb.diff < stop.crit)) break
#     niter <- niter + 1
#     if (!silent) setTxtProgressBar(pb, t)
#   }
#   if (!silent) close(pb)
#
#   res <- list(w = u[niter, ] * psi[niter], alpha = alpha,
#               lambda = xi[niter] / psi[niter], psi = psi[niter],
#               lower.bound = lower.bound, kernel = kernel,
#               X = X, y = y, H = H, lambda.list = xi / psi, niter = niter)
#   class(res) <- "ipriorVB"
#   res
# }
#
# # I-prior probit fitted
# fitted.ipriorVB <- function(object) {
#   w <- object$w
#   lambda <- object$lambda
#   alpha <- object$alpha
#
#   as.numeric(alpha + lambda * object$H %*% w)
# }
#
# # I-prior probit predict
# predict.ipriorVB <- function(object, newdata, ...) {
#   w <- object$w
#   lambda <- object$lambda
#   alpha <- object$alpha
#
#   if (object$kernel == "Canonical") H.tilde <- fnH2(object$X, newdata)
#   if (object$kernel == "FBM") H.tilde <- fnH3(object$X, newdata)
#
#   ystar.hat <- as.numeric(alpha + lambda * H.tilde %*% w)
#   y.hat <- rep(0, nrow(newdata)); y.hat[ystar.hat >= 0] <- 1
#   p.hat <- pnorm(ystar.hat)
#
#   list(y = y.hat, prob = p.hat)
# }
#
# # I-prior probit plot
# plot.ipriorVB <- function(x, niter.plot = NULL, levels = NULL, ...) {
#   if (is.null(niter.plot)) niter.plot <- length(x$lower.bound) - 1
#   tmp <- as.factor(x$y)
#   if (!is.null(levels)) levels(tmp) <- levels
#   lb <- x$lower.bound[1:niter.plot]
#   error.rate <- x$error.rate[1:niter.plot]
#   maximin <- max(lb) - min(lb)
#   maximin.inv <- 1 / maximin
#   error.rate.scaled <- error.rate * maximin + min(lb)
#   plot.df1 <- data.frame(Iteration = 1:niter.plot,
#                          lower = lb,
#                          error = error.rate.scaled)
#   plot.df2 <- data.frame(Observation = 1:length(x$ystar),
#                          p.hat = fitted(x)$prob,
#                          Class = tmp)
#
#   p1 <- ggplot(plot.df1) +
#     geom_point(aes(x = Iteration, y = error, col = "Error rate")) +
#     geom_line(aes(x = Iteration, y = error, col = "Error rate",
#                   linetype = "Error rate")) +
#     geom_point(aes(x = Iteration, y = lower, col = "Lower bound")) +
#     geom_line(aes(x = Iteration, y = lower, col = "Lower bound",
#                   linetype = "Lower bound")) +
#     scale_linetype_manual(name = NULL,
#                           values = c("Lower bound" = "longdash",
#                                      "Error rate" = "solid")) +
#     scale_colour_manual(name = NULL,
#                         values = c("Lower bound" = "black",
#                                    "Error rate" = "lightgoldenrod4")) +
#     scale_y_continuous(
#       "Lower bound",
#       sec.axis = sec_axis(~ (. - min(lb)) * maximin.inv, name = "Error rate")
#     ) +
#     theme(legend.position = "top")
#
#   p2 <- ggplot(plot.df2, aes(x = Observation, y = p.hat, col = Class)) +
#     geom_point() +
#     labs(y = "Fitted probabilities")
#
#   grid.arrange(p1, p2, ncol = 1, nrow = 2, heights = c(6, 4))
# }
#
# print.ipriorVB <- function(x, newdata = NULL, testdata = NULL) {
#   cat("iterations = ", x$niter)
#   cat("\nlower bound = ", x$lower.bound[x$niter])
#   cat("\nalpha = ", x$alpha)
#   cat("\nlambda = ", x$lambda)
#   cat("\npsi = ", x$psi)
# }
#
