plot_predict <- function(x) {
  # Args: x an ipriorMod2 object.
  plot.df <- as.data.frame(fitted(x)[1:2])
  ggplot(plot.df, aes(y, resid)) +
    geom_hline(yintercept = 0, col = "grey50", linetype = "dashed") +
    geom_point() +
    labs(x = "Fitted values", y = "Residuals") +
    theme_bw()
}

plot_fitted2 <- function(x, X.var = 1) {
  y.hat <- fitted(x)$y
  X <- x$kernL$Xl[[X.var]]
  plot.df <- data.frame(y.hat = y.hat, x = X,
                        y = as.numeric(x$kernL$y) + x$intercept)
  x.lab <- x$kernL$xname[X.var]
  y.lab <- x$kernL$yname

  ggplot(plot.df) +
    geom_point(aes(x, y)) +
    geom_line(aes(x, y.hat), col = "red3") +
    labs(x = x.lab, y = y.lab) +
    theme_bw()
}
