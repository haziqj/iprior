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

#' Plots for \code{ipriorMod} objects
#'
#' Three plots are produced by default: Plot of fitted regression curve, plot of
#' fitted values against residuals, and a QQ-plot of the residuals. Note that
#' The plots of fitted regression line can be shown only if the explanatory
#' variable is of dimension one (i.e. only \code{p = 1} explanatory variable
#' used).
#'
#' @param x An object of class \code{ipriorMod}.
#' @param plots Option to control which plots to show. The options are:
#'   \describe{\item{\code{all}}{(default) All three plots are
#'   shown.}\item{\code{allinone}}{All three plots are shown one one
#'   screen.}\item{\code{fitted}}{Only the fitted regression curve is
#'   shown.}\item{\code{diagnostic}}{The two diagnostic plots are
#'   shown.}\item{\code{residuals}}{Only the plot of fitted against residuals
#'   shown.}\item{\code{qqplot}}{Only the QQ-plot of the residuals is shown.}}
#' @param own.labels Logical, useful when categorical variables has factor names
#'   which are long, because these are used as the points of the plots.
#' @param ... No further arguments are passed, so this is not used here.
#'
#' @examples
#' # Straight line regression (Canonical RKHS)
#' mod.orange <- iprior(circumference ~ age, Orange)
#' plot(mod.orange)
#'
#' # Multilevel model type plots (Canonical + Pearson RKHS)
#' mod.tooth1 <- iprior(len ~ ., ToothGrowth)
#' mod.tooth2 <- iprior(len ~ . ^ 2, ToothGrowth)
#' par(mfrow = c(1,2))
#' plot(mod.tooth1, plots = "fitted")  # random intercept
#' plot(mod.tooth2, plots = "fitted")  # random slopes & intercept
#'
#' # One-dimensional smoothing (FBM RKHS with Hurst coef. 0.5)
#' mod.cars <- iprior(dist ~ speed, cars, model = list(kernel = "FBM"))
#' plot(mod.cars, plots = "fitted")
#'
#' @name plot
#' @export
plot.ipriorMod <- function(x,
                           plots = c("all", "allinone", "fitted", "diagnostic",
                                     "residuals", "qqplot"),
                           own.labels = FALSE, ...) {
  object <- x
	x <- object$ipriorKernel$x
	y <- object$ipriorKernel$Y
	whichPearson <- isPea(object$ipriorKernel$model$kernel)
	xnames <- object$ipriorKernel$model$xname
	yname <- object$ipriorKernel$model$yname
	yhat <- object$fitted
	resid <- object$residuals
	top3 <- order(abs(resid), decreasing = TRUE)[1:3]
	# colx <- c(RColorBrewer::brewer.pal(9, "Set1")[-6],
	#           RColorBrewer::brewer.pal(12, "Paired")[c(2,4,6,8,10,12)],
	#           RColorBrewer::brewer.pal(8,"Dark2"))
	# colx <- c(RColorBrewer::brewer.pal(8, "Dark2"),
	#           RColorBrewer::brewer.pal(8, "Set2"))
	colx <- c(RColorBrewer::brewer.pal(9, "Set1")[-9],
	          RColorBrewer::brewer.pal(8, "Dark2"))
	colx[6] <- RColorBrewer::brewer.pal(8, "Set2")[6]

	if (!is.numeric(plots)) {
		thisplot <- match.arg(plots)
		if ((thisplot == "all") | (thisplot == "allinone")) thisplot <- 1:3
		else if (thisplot == "fitted") thisplot <- c(1)
		else if (thisplot == "diagnostic") thisplot <- c(2,3)
		else if (thisplot == "residuals") thisplot <- c(2)
		else if (thisplot == "qqplot") thisplot <- c(3)
	}
	else{
		thisplot <- plots
		if (any((thisplot > 3) | (thisplot < 1))) stop("Must be either 1, 2 or 3.")
	}
	whichplot <- (1:3) %in% thisplot

	# Plot 1: Fitted regression curve --------------------------------------------
	which.cts <- which(!whichPearson)
	which.ctg <- which(whichPearson)

	if ((length(which.cts) == 1) | (length(which.ctg) > 0)) {
	  x.cts <- unlist(x[which.cts], use.names = FALSE)
	  if (length(x.cts) != length(y)) stop("X variable more than 1 dimensions.")
	  if (length(which.ctg) == 0) {
		  # Plot only continuous variables -----------------------------------------
			plot1 <- function(z) {
				xorder <- order(x.cts)
				plot(x = x.cts, y = y, xlab = xnames[which.cts], ylab = yname,
				     main = "Fitted regression curve", cex = 0.55)
				lines(x = x.cts[xorder], y = yhat[xorder], col = colx[1], lwd = 1.55)
			}
		} else {
		  # Define categorical variables -------------------------------------------
		  if (length(which.ctg) == 1) x.ctg <- as.factor(x[[which.ctg]])
			else {
				x.ctg <- as.data.frame(x[which.ctg])
				x.ctg <- interaction(x.ctg)
				# x.ctg <- x.ctg[,2]
			}
			if (length(which.cts) == 0) {
			  # Plot only categorical variables --------------------------------------
				yhat.unq <- unique(yhat)
				plotlvl <- levels(x.ctg)
				grp <- as.numeric(x.ctg)
				plotlvl <- plotlvl[unique(grp)]
				grp <- as.numeric(factor(grp))
				if (own.labels) plotlvl <- unique(grp)
				plot1 <- function(z) {
					plot(x = grp, y = y, type = "n", xlab = xnames[which.ctg],
					     ylab = yname, main = "Fitted regression curve", xaxt = "n",
					     xlim = c(0.5, length(unique(grp)) + 0.5))
					for (i in unique(grp)) {
						text(x = grp[grp == i], y[grp == i], plotlvl[i],
						     col = colx[(i - 1) %% 16 + 1], cex = 0.8)
						abline(a = yhat.unq[i], b = 0, col = colx[(i - 1) %% 16 + 1])
					}
				}
			} else {
				# Multilevel plots -----------------------------------------------------
				plotlvl <- levels(x.ctg)
				grp <- as.numeric(x.ctg)
				plotlvl <- plotlvl[unique(grp)]
				grp <- as.numeric(factor(grp))
				if (own.labels) plotlvl <- unique(grp)
				plot1 <- function(z) {
					plot(x = x.cts, y = y, type = "n", xlab = xnames[which.cts],
					     ylab = yname, main = "Fitted regression curve")
					for (i in unique(grp)) {
						xorder <- order(x.cts[grp == i])
						text(x = x.cts[grp == i], y[grp == i], plotlvl[i],
						     col = colx[(i - 1) %% 16 + 1], cex = 0.55)
						lines(x = x.cts[grp == i][xorder], yhat[grp == i][xorder],
						      col = colx[(i - 1) %% 16 + 1], lwd = 1.55)
					}
				}
			}
		}
	} else {
		whichplot[1] <- FALSE
		message("Fitted plot not generated because dimensions of continuous x variables is greater than 1.")
	}

	# Plot 2: Fitted vs. Residuals -----------------------------------------------
	plot2 <- function(z){
		plot(x = yhat, y = resid, xlab = "Fitted values", ylab = "Residuals",
		     main = "Fitted vs. Residuals")
		text(yhat[top3], resid[top3], as.character(top3), pos = 4, cex = 0.8)
		abline(0, 0, col = 200, lty = 3, lwd = 1.5)
	}

	# Plot 3: Normal Q-Q plot ----------------------------------------------------
	plot3 <- function(z){
		tmp <- qqnorm(scale(resid), ylab = "Standardised Residuals")
		updown <- c(3,3,3); updown[which(tmp$y[top3,] > 0)] <- 1
		text(tmp$x[top3], tmp$y[top3, ], as.character(top3), pos = updown, cex = 0.8)
		abline(0, 1, col = 200, lty = 3, lwd = 1.5)
	}

	# Plot the plots -------------------------------------------------------------
	if (!is.numeric(plots) && match.arg(plots) == "allinone") {
		dev.new(width = 14, height = 5.3)
		par(mfrow = c(1,3))
		for (i in which(whichplot)) {
			do.call(paste0("plot", i), list(i))
		}
	}
	else{
		timetostop <- length(which(whichplot))
		tts.count <- 1
		for (i in which(whichplot)) {
			do.call(paste0("plot", i), list(i))
			if (!(tts.count == timetostop)) readline("Hit <Return> to see next plot")
			tts.count <- tts.count + 1
		}
	}
}
