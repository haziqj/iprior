## Canonical kernel function
fn.H2 <- function(x, y=NULL){ #takes in vector of covariates
	if(is.null(y)) y <- x
	z <- c(x,y)
	tmp <- tcrossprod(y, x)
	tmp
}

## (centred) Canonical kernel function
fn.H2a <- function(x, y=NULL){ #takes in vector of covariates
	x <- as.numeric(x)
	if(is.null(y)) y <- x
	else y <- as.numeric(y)
	z <- c(x,y)
	xbar <- mean(z)
	tmp <- tcrossprod(y-xbar, x-xbar)
	tmp
}

## Pearson kernel
fn.H1 <- function(x, y=NULL){	#takes the group vector as argument
	x <- as.numeric(x)
	if(is.null(y)) y <- x
	else y <- as.numeric(y)
	z <- c(x,y)
	N3 <- length(z)
	n1 <- as.data.frame(table(x))[,2]
	n2 <- as.data.frame(table(y))[,2]
	prop <- table(z)/N3
	
	tmp <- matrix(NA, nr=sum(n2), nc=sum(n1))
	for(i in 1:sum(n2)){
		for(j in 1:sum(n1)){
			tmp[i,j] <- as.numeric(y[i] == x[j]) / prop[names(prop) == y[i]] - 1
		}
	}
	tmp
}