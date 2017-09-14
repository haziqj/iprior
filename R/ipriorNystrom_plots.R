#' @export
plot.ipriorMod_Nystrom <- function(x, ...) {
  plot_fitted(x)
}

#' New plots
#'
#' @param object an iprior object
#'
#' @name iprior_plot
#' @export
plot_fitted <- function(object) {
  list2env(object, environment())
  list2env(ipriorKernel, environment())
  list2env(model, environment())

  xval <- x[[1]][order(Nys.samp)]
  yval <- Y[order(Nys.samp)]
  x.order <- order(xval)
  y.hat <- fitted(object)
  df.plot <- data.frame(y = yval[x.order], x = xval[x.order], fitted = y.hat[x.order])

  ggplot(df.plot) +
    geom_point(aes(x = x, y = y)) +
    geom_line(aes(x = x, y = fitted), col = ggColPal(1), size = 1.3) +
    theme_bw()
}

# plot_multilevel_2 <- function(object) {
#   class(object) <- "ipriorMod"
#   object$ipriorKernel <- .reorder_ipriorKernel(object$ipriorKernel, order(object$ipriorKernel$model$Nys.samp))
#   object$fitted <- fitted(object)
#   object$residuals <- object$ipriorKernel$Y[order(object$ipriorKernel$model$Nys.kern)] - object$fitted
#   plot(object)
# }

#' @rdname iprior_plot
#' @export
plot_multilevel <- function(object) {
  list2env(object, environment())
  list2env(ipriorKernel, environment())
  list2env(model, environment())

  x <- object$ipriorKernel$x
  whichPearson <- isPea(object$ipriorKernel$model$kernel)
  which.cts <- which(!whichPearson)
  which.ctg <- which(whichPearson)
  x.cts <- unlist(x[which.cts], use.names = FALSE)
  x.ctg <- as.factor(x[[which.ctg]])
  y <- object$ipriorKernel$Y

  plot.df <- data.frame(x.cts, x.ctg, y)
  plot.df <- plot.df[order(Nys.samp), ]
  plot.df$y.hat <- fitted(object)

  ggplot(plot.df) +
    geom_text(aes(x = x.cts, y = y, col = x.ctg, label = x.ctg), size = 3) +
    geom_line(aes(x = x.cts, y = y.hat, col = x.ctg), size = 1.3) +
    theme_bw() +
    theme(legend.position = "none")

}
