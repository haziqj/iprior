################################################################################
#
#   iprior: Linear Regression using I-priors
#   Copyright (C) 2017  Haziq Jamil
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
#' @param x An \code{ipriorMod} object.
#' @param X.var The index of the X variable to plot.
#' @param cred.bands Logical. Plot the confidence intervals? Defaults to \code{TRUE}.
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
plot_predict <- function(x) {
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
plot_fitted_multilevel <- function(x, X.var = 1, grp = 1, facet = c(2, 3),
                                   cred.bands = TRUE, show.legend = TRUE,
                                   show.points = TRUE, x.lab = NULL,
                                   y.lab = NULL, grp.lab = NULL) {
  fit <- fitted(x, intervals = cred.bands)
  y.hat <- fit$y

  which.pearson <- x$ipriorKernel$which.pearson
  cat.x <- which(which.pearson)
  cts.x <- which(!which.pearson)
  X      <- x$ipriorKernel$Xl[[cts.x[X.var]]]
  grp    <- x$ipriorKernel$Xl[[cat.x[grp]]]
  plot.df <- data.frame(y.hat = y.hat, x = X, grp = grp, y = get_y(x))
  if (length(cat.x) == 2) {
    plot.df <- cbind(plot.df, facet1 = x$ipriorKernel$Xl[[cat.x[facet[1]]]])
  }
  if (length(cat.x) == 3) {
    plot.df <- cbind(plot.df, facet2 = x$ipriorKernel$Xl[[cat.x[facet[2]]]])
  }

  if (is.null(x.lab)) x.lab <- x$ipriorKernel$xname[X.var]
  if (is.null(y.lab)) y.lab <- x$ipriorKernel$yname
  if (is.null(grp.lab)) grp.lab <- names(cat.x[grp])
  nys.check <- is.ipriorKernel_nys(x$ipriorKernel)

  p <- ggplot(plot.df)

  if (isTRUE(cred.bands)) {
    p <- p + geom_ribbon(aes(x = X, ymin = fit$lower, ymax = fit$upper,
                             fill = grp), alpha = 0.15) +
      scale_fill_discrete(name = grp.lab)
  }

  if (isTRUE(show.points)) {
    p <- p + geom_point(aes(x, y, col = grp))
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

  if (!isTRUE(show.legend)) {
    p <- p + theme(legend.position = "none")
  }

  p
}

#' @rdname plot.ipriorMod
#' @export
plot_fitted <- function(x, X.var = 1, cred.bands = TRUE) {
  fit <- fitted(x, intervals = cred.bands)
  y.hat <- fit$y
  X <- x$ipriorKernel$Xl[[X.var]]
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
        data = plot.df[seq_len(x$ipriorKernel$nystroml$nys.size), ], aes(x, y)
      )
  } else {
    p <- p + geom_point(aes(x, y))
  }

  p + geom_line(aes(x, y.hat), col = "red3") +
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
