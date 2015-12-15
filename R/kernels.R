## Canonical kernel function
fn.H2 <- function(x, y=NA){ #takes in vector of covariates
	if(is.na(sum(y))) y <- x
	z <- c(x,y)
	tmp <- tcrossprod(y, x)
	tmp
}

## (centred) Canonical kernel function
fn.H2a <- function(x, y=NA){ #takes in vector of covariates
	if(is.na(sum(y))) y <- x
	z <- c(x,y)
	xbar <- mean(z)
	tmp <- tcrossprod(y-xbar, x-xbar)
	tmp
}

## Pearson kernel
fn.H1 <- function(x, y=NA){	#takes the group vector as argument
	if(is.na(sum(y))) y <- x
	z <- c(x,y)
	N3 <- length(z)
	n1 <- as.data.frame(table(x))[,2]
	n2 <- as.data.frame(table(y))[,2]
	prop <- as.data.frame(table(z))[,2]/N3
	
	indexi <- c(0,cumsum(n2)); indexj <- c(0,cumsum(n1))
	tmp <- matrix(-1, nr=sum(n2), nc=sum(n1))
	for(k in 1:min(length(n1),length(n2))) tmp[(indexi[k]+1):indexi[k+1],(indexj[k]+1):indexj[k+1]] <- matrix(1/prop[k]-1, nr=n2[k], nc=n1[k])
	tmp
}