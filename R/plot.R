#' @export
plot.ipriorMod <- function(x,
                           plots = c("all", "allinone", "fitted", "diagnostic",
                                     "residuals", "qqplot"),
                           UseOwnLabels = FALSE, ...) {
  object <- x
	x <- object$ipriorKernel$x
	y <- object$ipriorKernel$Y
	p <- object$ipriorKernel$p
	whichPearson <- isPea(object$ipriorKernel$model$kernel)
	no.plot <- sum(!whichPearson)
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

	if(!is.numeric(plots)){
		thisplot <- match.arg(plots)
		if ((thisplot == "all") | (thisplot == "allinone")) thisplot <- 1:3
		else if (thisplot == "fitted") thisplot <- c(1)
		else if (thisplot == "diagnostic") thisplot <- c(2,3)
		else if (thisplot == "residuals") thisplot <- c(2)
		else if (thisplot == "qqplot") thisplot <- c(3)
	}
	else{
		thisplot <- plots
		if (any((thisplot > 3) | (thisplot < 1))) stop("Must be between 0 and 1.")
	}
	whichplot <- (1:3) %in% thisplot

	### Plot 1: Fitted values
	cts.vars <- which(!whichPearson)
	ctg.vars <- which(whichPearson)

	if ((length(cts.vars) == 1) | (length(ctg.vars) > 0)) {
		if (length(ctg.vars) == 0) {
		  # Plot only continuous variables -----------------------------------------
			x.cts <- unlist(x[cts.vars], use.names=F)
			if (length(x.cts) != length(y)) stop("X variable not a vector.")
			plot1 <- function(z){
				xorder <- order(x.cts)
				plot(x=x.cts, y=y, xlab=xnames[cts.vars], ylab=yname, main="Fitted regression curve", cex=0.55)
				lines(x=x.cts[xorder], y=yhat[xorder], col=colx[1], lwd=1.55)
			}
		} else {
		  # Define categorical variables -------------------------------------------
		  if (length(ctg.vars) == 1) x.ctg <- as.factor(x[[ctg.vars]])
			else {
				x.ctg <- as.data.frame(x[ctg.vars])
				x.ctg <- interaction(x.ctg)
				# x.ctg <- x.ctg[,2]
			}
			if (length(cts.vars) == 0) {
			  # Plot only categorical variables --------------------------------------
				yhat.unq <- unique(yhat)
				plotlvl <- levels(x.ctg)
				grp <- as.numeric(x.ctg)
				plotlvl <- plotlvl[unique(grp)]
				grp <- as.numeric(factor(grp))
				if (UseOwnLabels) plotlvl <- unique(grp)
				plot1 <- function (z) {
					plot(x=grp, y=y, type="n", xlab=xnames[ctg.vars], ylab=yname, main="Fitted regression curve", xaxt="n", xlim=c(0.5, length(unique(grp))+0.5))
					for (i in unique(grp)) {
						text(x=grp[grp==i], y[grp==i], plotlvl[i], col=colx[(i-1)%%16+1], cex=0.8)
						abline(a=yhat.unq[i], b=0, col=colx[(i-1)%%16+1])
					}
				}
			} else {
				# Multilevel plots -----------------------------------------------------
			  x.cts <- as.numeric(x[[cts.vars]])
				if (length(x.cts) != length(y)) stop("X variable not a vector.")
				plotlvl <- levels(x.ctg)
				grp <- as.numeric(x.ctg)
				plotlvl <- plotlvl[unique(grp)]
				grp <- as.numeric(factor(grp))
				if (UseOwnLabels) plotlvl <- unique(grp)
				plot1 <- function(z) {
					plot(x = x.cts, y = y, type = "n", xlab = xnames[cts.vars],
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
		message("Fitted plot not generated because number of continuous x variables is greater than 1.")
	}

	### Plot 2: Fitted vs. Residuals
	plot2 <- function(z){
		plot(x=yhat, y=resid, xlab="Fitted values", ylab="Residuals", main="Fitted vs. Residuals")
		text(yhat[top3], resid[top3], as.character(top3), pos=4, cex=0.8)
		abline(0, 0, col=200, lty=3, lwd=1.5)
	}

	### Plot 3: Normal Q-Q plot
	plot3 <- function(z){
		tmp <- qqnorm(scale(resid), ylab="Standardised Residuals")
		updown <- c(3,3,3); updown[which(tmp$y[top3,] > 0)] <- 1
		text(tmp$x[top3], tmp$y[top3,], as.character(top3), pos=updown, cex=0.8)
		abline(0, 1, col=200, lty=3, lwd=1.5)
	}

	## Plot
	if(!is.numeric(plots) && match.arg(plots) == "allinone"){
		dev.new(width=14, height=5.3)
		par(mfrow=c(1,3))
		for(i in which(whichplot)){
			do.call(paste0("plot", i), list(i))
		}
	}
	else{
		timetostop <- length(which(whichplot))
		tts.count <- 1
		for(i in which(whichplot)){
			do.call(paste0("plot", i), list(i))
			if(!(tts.count == timetostop)) readline("Hit <Return> to see next plot")
			tts.count <- tts.count + 1
		}
	}

}
