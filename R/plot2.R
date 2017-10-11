plot_predict <- function(x) {
  # Args: x an ipriorMod2 object.
  plot.df <- as.data.frame(fitted(x)[1:2])
  ggplot(plot.df, aes(y, resid)) +
    geom_hline(yintercept = 0, col = "grey50", linetype = "dashed") +
    geom_point() +
    labs(x = "Fitted values", y = "Residuals") +
    theme_bw()
}

plot_fitted2 <- function(x, X.var = 1, ci = TRUE) {
  fit <- fitted(x, intervals = ci)
  y.hat <- fit$y
  X <- x$ipriorKernel$Xl[[X.var]]
  plot.df <- data.frame(y.hat = y.hat, x = X,
                        y = as.numeric(x$ipriorKernel$y) + x$intercept)
  x.lab <- x$ipriorKernel$xname[X.var]
  y.lab <- x$ipriorKernel$yname
  nys.check <- is.ipriorKernel_nys(x$ipriorKernel)

  p <- ggplot(plot.df)

  if (isTRUE(ci)) {
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
