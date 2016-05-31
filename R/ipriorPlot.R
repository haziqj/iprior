###
### Plot function
###

plot.iprior <- function(object, UseOwnLabels=F, plots=c("all", "allinone", "fitted", "diagnostic", "residuals", "qqplot"), ...){
	require(RColorBrewer, quietly=T)
	x <- object$xval
	y <- object$yval
	p <- object$p
	no.plot <- sum(!object$whichPearson)
	xnames <- colnames(x)
	yname <- object$yname
	yhat <- fitted(object)
	resid <- residuals(object)
	top3 <- order(abs(resid), decreasing=T)[1:3]
	colx <- c(brewer.pal(9, "Set1")[-6], brewer.pal(12, "Paired")[c(2,4,6,8,10,12)], brewer.pal(8,"Dark2"))

	if(!is.numeric(plots)){
		thisplot <- match.arg(plots)
		if((thisplot == "all") | (thisplot == "allinone")) thisplot <- 1:3
		else if(thisplot == "fitted") thisplot <- c(1)
		else if(thisplot == "diagnostic") thisplot <- c(2,3)
		else if(thisplot == "residuals") thisplot <- c(2)
		else if(thisplot == "qqplot") thisplot <- c(3)
	} 
	else{
		thisplot <- plots
		if(any((thisplot > 3) | (thisplot < 1))) stop("Must be between 0 and 1.")
	} 
	whichplot <- (1:3) %in% thisplot

	### Plot 1: Fitted values
	cts.vars <- which(!object$whichPearson)
	ctg.vars <- which(object$whichPearson)

	if(length(c(cts.vars, ctg.vars)) <= 2){
		x.cts <- x[,cts.vars]
		if(length(cts.vars) == 0){
			yhat.unq <- unique(yhat)
			x.ctg <- x[,ctg.vars]
			plotlvl <- levels(x.ctg)
			grp <- as.numeric(x.ctg)
			if(UseOwnLabels) plotlvl <- unique(grp)
			plot1 <- function(z){
				plot(x=grp, y=y, type="n", xlab=xnames[ctg.vars], ylab=yname, main="Fitted regression curve", xaxt="n", xlim=c(0.5, length(unique(grp))+0.5))
				for(i in unique(grp)){
					text(x=grp[grp==i], y[grp==i], plotlvl[i], col=colx[(i-1)%%22+1], cex=0.8)
					abline(a=yhat.unq[i], b=0, col=colx[(i-1)%%22+1])
				}
			}			
		}
		else if(length(ctg.vars) == 0){
			plot1 <- function(z){
				xorder <- order(x.cts)							
				plot(x=x.cts, y=y, xlab=xnames[cts.vars], ylab=yname, main="Fitted regression curve")
				lines(x=x.cts[xorder], y=yhat[xorder], col=colx[1])			
			}
		}
		else{
			x.ctg <- x[,ctg.vars]
			plotlvl <- levels(x.ctg)
			grp <- as.numeric(x.ctg)
			if(UseOwnLabels) plotlvl <- unique(grp)
			plot1 <- function(z){
				plot(x=x.cts, y=y, type="n", xlab=xnames[cts.vars], ylab=yname, main="Fitted regression curve")
				for(i in unique(grp)){
					xorder <- order(x.cts[grp==i])				
					text(x=x.cts[grp==i], y[grp==i], plotlvl[i], col=colx[(i-1)%%22+1], cex=0.8)
					lines(x=x.cts[grp==i][xorder], yhat[grp==i][xorder], col=colx[(i-1)%%22+1])
				}
			}
		}
	}
	else{
		whichplot[3] <- F
		message("Fitted plot not generated because number of x variables is greater than 1.")
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