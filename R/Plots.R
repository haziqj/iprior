################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2018  Haziq Jamil
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

#' Plots for I-prior models
#'
#' There are three types of plots that are currently written in the package:
#' \describe{ \item{\code{plot_fitted}}{Plot the fitted regression line with
#' credibility bands.} \item{\code{plot_predict}}{Plot residuals against fitted
#' values.} \item{\code{plot_iter}}{Plot the progression of the log-likelihood
#' value over time.} } The S3 method \code{plot} for class \code{ipriorMod}
#' currently returns \code{plot_fitted}.
#'
#' See ggplot2 documentation for the plotting parameters.
#'
#' @param x An \code{ipriorMod} object.
#' @param X.var The index of the X variable to plot.
#' @param cred.bands Logical. Plot the confidence intervals? Defaults to
#'   \code{TRUE}.
#' @param niter.plot (Optional) Vector of length at most two, indicating the
#'   start and end points of the iterations to plot.
#' @param lab.pos Adjust the position of the log-likelihood label.
#' @param ... Not used
#' @param grp The index of the groups.
#' @param facet The index of the X variable in which to facet. This is a vector
#'   of maximum length 2.
#' @param show.legend Logical. Show legend?
#' @param show.points Logical. Show data points?
#' @param x.lab (Optional) X axis label.
#' @param y.lab (Optional) Y axis label.
#' @param grp.lab (Optional) The name for the groups, which is also the legend
#'   title.
#' @param grp.var Index of the grouping variable for multilevel plots.
#' @param extrapolate Logical. Extend the fitted regression line to fill the
#'   plot?
#' @param size Size of the fitted line
#' @param linetype Type of the fitted line
#' @param draws Number of draws for posterior predictive check.
#'
#' @export
plot.ipriorMod <- function(x, ...) {
  which.pearson <- x$ipriorKernel$which.pearson
  if (sum(which.pearson) > 0) {
    if (sum(which.pearson) == length(which.pearson)) {
      stop("Not available yet.", call. = FALSE)
    } else {
      plot_fitted_multilevel(x)
    }
  } else {
    plot_fitted(x)
  }
}

#' @rdname plot.ipriorMod
#' @export
plot_resid <- function(x) {
  # Args: x an ipriorMod object.
  plot.df <- as.data.frame(fitted(x)[1:2])
  ggplot(plot.df, aes(y, resid)) +
    geom_hline(yintercept = 0, col = "grey50", linetype = "dashed") +
    geom_point() +
    labs(x = "Fitted values", y = "Residuals") +
    theme_bw()
}

#' @rdname plot.ipriorMod
#' @export
plot_fitted_multilevel <- function(x, X.var = 1, grp.var = 1, facet = c(2, 3),
                                   cred.bands = TRUE, show.legend = TRUE,
                                   show.points = TRUE, x.lab = NULL,
                                   y.lab = NULL, grp.lab = NULL,
                                   extrapolate = FALSE) {
  which.pearson <- x$ipriorKernel$which.pearson
  cat.x <- which(which.pearson)
  cts.x <- which(!which.pearson)
  X      <- x$ipriorKernel$Xl[[cts.x[X.var]]]
  grp    <- x$ipriorKernel$Xl[[cat.x[grp.var]]]
  fit <- fitted(x, intervals = cred.bands)
  y.hat <- fit$y
  plot.df <- data.frame(y.hat = y.hat, x = X, grp = grp, y = get_y(x))
  plot.df
  if (isTRUE(extrapolate)) {
    x.min <- min(X)
    x.max <- max(X)
    X.ext <- seq(x.min, x.max, length = 20)
    grp.ext <- rep(unique(grp), each = 20)
    ext.df <- data.frame(X.ext, grp.ext)
    if (!is.null(x$ipriorKernel$formula)) {
      x.pos <- c(cts.x[X.var], cat.x[grp.var])
      colnames(ext.df)[1:2] <- attr(x$ipriorKernel$terms, "term.labels")[x.pos]
      fit.ext <- predict_iprior(x, ext.df, NULL, cred.bands, 0.05)
      y.hat.ext <- fit.ext$y
      plot.df.ext <- data.frame(y.hat = y.hat.ext, x = X.ext, grp = grp.ext,
                                y = NA)
    } else {
      stop("Not implemented yet. Sorry!")
    }
    plot.df <- rbind(plot.df, plot.df.ext)
  }

  if (length(cat.x) == 2) {
    plot.df <- cbind(plot.df, facet1 = x$ipriorKernel$Xl[[cat.x[facet[1]]]])
  }
  if (length(cat.x) == 3) {
    plot.df <- cbind(plot.df, facet2 = x$ipriorKernel$Xl[[cat.x[facet[2]]]])
  }

  if (is.null(x.lab)) x.lab <- x$ipriorKernel$xname[cts.x[X.var]]
  if (is.null(y.lab)) y.lab <- x$ipriorKernel$yname
  if (is.null(grp.lab)) grp.lab <- x$ipriorKernel$xname[cat.x[X.var]]
  nys.check <- is.ipriorKernel_nys(x$ipriorKernel)

  p <- ggplot(plot.df)

  if (isTRUE(cred.bands)) {
    plot.cred.bands.df <- data.frame(
      x = X, lower = fit$lower, upper = fit$upper, grp = grp
    )
    if (isTRUE(extrapolate)) {
      plot.cred.bands.df.ext <- data.frame(
        x = X.ext, lower = fit.ext$lower, upper = fit.ext$upper, grp = grp.ext
      )
      plot.cred.bands.df <- rbind(plot.cred.bands.df, plot.cred.bands.df.ext)
    }
    p <- p + geom_ribbon(data = plot.cred.bands.df, alpha = 0.15,
                         aes(x, ymin = lower, ymax = upper, fill = grp)) +
      scale_fill_discrete(name = grp.lab)
  }

  if (isTRUE(show.points)) {
    p <- p + geom_point(aes(x, y, col = grp), na.rm = TRUE)
  }

  p <- p +
    geom_line(aes(x, y.hat, col = grp)) +
    labs(x = x.lab, y = y.lab) +
    scale_colour_discrete(name = grp.lab) +
    theme_bw()

  if (length(cat.x) == 2) {
    p <- p + facet_grid(. ~ facet1)
  }
  if (length(cat.x) == 3) {
    p <- p + facet_grid(facet2 ~ facet1)
  }
  if (as.character(paste0(facet, collapse = "")) == "1") {
    p <- p + facet_wrap( ~ grp, nrow = 2)
  }

  if (!isTRUE(show.legend)) {
    p <- p + theme(legend.position = "none")
  }

  p
}

#' @rdname plot.ipriorMod
#' @export
plot_fitted <- function(x, X.var = 1, cred.bands = TRUE, size = 1,
                        linetype = "solid") {
  fit <- fitted(x, intervals = cred.bands)
  y.hat <- fit$y
  X <- x$ipriorKernel$Xl[[X.var]]
  if (!is.null(dim(X))) {
    if (ncol(X) > 1) X <- X[, X.var]
  }
  plot.df <- data.frame(y.hat = y.hat, x = X, y = get_y(x))
  x.lab <- x$ipriorKernel$xname[X.var]
  y.lab <- x$ipriorKernel$yname
  nys.check <- is.ipriorKernel_nys(x$ipriorKernel)

  p <- ggplot(plot.df)

  if (isTRUE(cred.bands)) {
    p <- p + geom_ribbon(aes(x = X, ymin = fit$lower, ymax = fit$upper),
                         fill = "grey", alpha = 0.5)
  }

  if (isTRUE(nys.check)) {
    p <- p + geom_point(aes(x, y), alpha = 0.15) +
    #   geom_point(
    #     data = plot.df[seq_len(x$ipriorKernel$nystroml$nys.size), ], aes(x, y),
    #     size = 2.5, col = "darkorange"
    # ) +
      geom_point(
        data = plot.df[seq_len(x$ipriorKernel$nystroml$nys.size), ], aes(x, y),
        size = 2, shape = 1, stroke = 1
      )
  } else {
    p <- p + geom_point(aes(x, y))
  }

  p + geom_line(aes(x, y.hat), col = "red3", size = size, linetype = linetype) +
    labs(x = x.lab, y = y.lab) +
    theme_bw()
}

#' @rdname plot.ipriorMod
#' @export
plot_iter <- function(x, niter.plot = NULL, lab.pos = c("up", "down")) {
  # Same code from iprobit, hence the lb references.

  if (x$niter < 2) stop("Nothing to plot.")

  lab.pos <- match.arg(lab.pos, c("up", "down"))
  if (lab.pos == "up") lab.pos <- -0.5
  else lab.pos <- 1.5

  lb.original <- x$loglik
  if (is.null(niter.plot)) niter.plot <- c(1, length(lb.original))
  else if (length(niter.plot) == 1) niter.plot <- c(1, niter.plot)
  niter.plot <- niter.plot[1]:niter.plot[2]
  lb <- lb.original[niter.plot]
  plot.df <- data.frame(Iteration = niter.plot, lb = lb)
  time.per.iter <- x$time$time / x$niter
  if (time.per.iter < 0.001) time.per.iter <- 0.001
  lb.lab <- rep("", length(lb))
  lb.lab[length(lb)] <- round(lb[length(lb)], 2)

  ggplot(plot.df, aes(x = Iteration, y = lb, label = max(lb))) +
    geom_line(col = "grey60") +
    geom_point() +
    geom_hline(yintercept = max(lb.original), linetype = 2, col = "red") +
    scale_x_continuous(
      sec.axis = sec_axis(~ . * time.per.iter, name = "Time (seconds)"),
      breaks = scales::pretty_breaks(n = min(5, ifelse(x$niter == 2, 1, x$niter)))
    ) +
    geom_text(aes(label = lb.lab), vjust = 1.5, size = 3.7) +
    annotate("text", col = "red3", x = niter.plot[1], y = max(lb.original),
             vjust = lab.pos, label = round(max(lb.original), 2), size = 3.7) +
    labs(y = "Log-likelihood") +
    theme_bw()
}

#' @rdname plot.ipriorMod
#' @export
plot_ppc <- function(x, draws = 100) {
  check_and_get_ipriorMod(x)
  Hlam <- get_Hlam(x$ipriorKernel, x$theta)
  psi <- theta_to_psi(x$theta, x$ipriorKernel)
  if (is.nystrom(x)) {
    eigen_Hlam_nys(Hlam, environment())  # assign u and V to environment
  } else {
    eigen_Hlam(Hlam, environment())  # assign u and V to environment
  }
  z <- psi * u ^ 2 + 1 / psi
  Vy.hat <- Hlam %*% A_times_a(1 / z, V, Hlam) + diag(1 / psi, x$ipriorKernel$n)
  y.hat <- fitted(x)$y
  ppc <- t(mvtnorm::rmvnorm(draws, mean = y.hat, sigma = Vy.hat))
  melted.ppc <- suppressMessages(reshape2::melt(data.frame(ppc = ppc)))

  p <- ggplot() +
    scale_x_continuous(breaks = NULL, name = expression(italic(y))) +
    scale_y_continuous(breaks = NULL) +
    geom_line(data = melted.ppc, stat = "density", alpha = 0.5,
              aes(x = value, group = variable, col = "yrep", size = "yrep")) +
    geom_line(aes(x = get_y(x), col = "y", size = "y"), stat = "density") +
    theme(legend.position = "bottom") +
    scale_colour_manual(
      name = NULL, labels = c("Observed", "Replications"),
      values = c("grey10", "steelblue3")
    ) +
    scale_size_manual(
      name = NULL, labels = c("Observed", "Replications"),
      values = c(1.1, 0.19)
    ) +
    labs(y = "Density", title = "Posterior predictive density check") +
    theme_bw() +
    theme(legend.position = c(0.9, 0.5))

  p
}

# plot_loglik <- function(x, xlim, ylim) {
#   check_and_get_ipriorMod(x)
#   if (missing(xlim)) xlim <- x$theta[1] + c(-1, 1) * 3
#   if (missing(ylim)) ylim <- x$theta[2] + c(-1, 1) * 3
#   theta.x <- seq(xlim[1], xlim[2], length = 50)
#   theta.y <- seq(ylim[1], ylim[2], length = 50)
#   plot.df <- expand.grid(theta.x, theta.y)
#   z <- seq_len(nrow(plot.df))
#   for (i in seq_along(z)) {
#     z[i] <- logLik(x$ipriorKernel, as.numeric(plot.df[i, ]))
#   }
#   plot.df <- cbind(plot.df, z)
#
#   ggplot(plot.df) +
#     geom_raster(aes(Var1, Var2, fill = z)) +
#     geom_contour(aes(Var1, Var2, z = z)) +
#     geom_vline(xintercept = x$theta[1]) +
#     geom_hline(yintercept = x$theta[2]) +
#     scale_fill_gradient(low = "white", high = "#00BFC4") +
#     theme_bw()
# }
