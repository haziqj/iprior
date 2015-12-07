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