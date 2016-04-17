## Canonical kernel function
fn.H2 <- function(x, y=NULL){ #takes in vector of covariates
	if(is.null(y)) y <- x
	tmp <- tcrossprod(y, x)
	tmp
}

## (centred) Canonical kernel function
fn.H2a <- function(x, y=NULL){ #takes in vector of covariates
	x <- as.numeric(x)
	if(is.null(y)) y <- x
	else y <- as.numeric(y)
	xbar <- mean(x)
	tmp <- tcrossprod(y-xbar, x-xbar)
	tmp
}

## Pearson kernel
fn.H1 <- function(x, y=NULL){ #takes in vectors of type factors
		
	fn.index <- function(k, tmp){
  		w <- tmp[[k]]
 		cbind(
    	row = rep(w, times=length(w):1 ) ,
    	col = unlist(lapply(1:length(w), function(x) c(NA,w)[-(0:x)])))
	}
	
	ytmp <- y
	if(is.null(ytmp)) y <- x
	## check if the inputs are non-factors, suggest not to use Pearson kernel
	if(any(!is.factor(x), !is.factor(y))) warning("Non-factor type vector used with Pearson kernel.", call.=F)
	
	z <- unlist(list(x,y)); z <- as.numeric(z)				#simply doing c(x,y) messes with the factors
	x <- z[1:length(x)]; y <- z[(length(x)+1):length(z)]
	if(any(is.na(match(y,x)))) stop("The vector y contains elements not belonging to x.")
	prop <- table(x)/length(x)

	
	tmpx <- unique(lapply(sort(x), function(k) which(x == k)))	
	indexx <- lapply(1:length(unique(x)), fn.index, tmp=tmpx)
	nx <- length(unique(x))
	#prop[as.numeric(names(prop)) > nx] <- Inf
	#indexx[[nx+1]] <- cbind(row=1:length(x), col=1:length(x))
	if(is.null(ytmp)) indexy <- indexx
	else{
		tmpy <- unique(lapply(sort(y), function(k) which(y == k)))
		indexy <- lapply(1:length(unique(y)), fn.index, tmp=tmpy)	
	}
	unqy <- sort(unique(y))

	tmp.mat <- matrix(-1, nr=length(y), nc=length(x))
	for(i in 1:length(unqy)){
		rowind <- indexy[[i]][,1]			#grabs row index from y
		colind <- indexx[[ unqy[i] ]][,2]	#and corresponding col index from x
		tmp.mat[rowind, colind] <- 1 / prop[ unqy[i] ] - 1
	}
	tmp.mat
}

## FBM kernel (to be developed)
fn.H3 <- function(x, y=NULL, gamma=NULL){ #takes in vector of covariates
	if(is.null(gamma)) gamma <- 0.5
	x <- as.numeric(x)
	n <- length(x)
	
	if(is.null(y)){
		w <- 1:n
		index <- cbind( row = rep(w, times=(length(w)-1):0 ) ,
						col = unlist(lapply(1:(length(w)-1), function(x) w[-(1:x)])) )
		index <- index[order(index[,2]),]				
		tmp2 <- abs(x[index[,1]])^(2*gamma) + abs(x[index[,2]])^(2*gamma) - abs(x[index[,1]]-x[index[,2]])^(2*gamma)
		mat0 <- diag(0,n)
		mat0[upper.tri(mat0)] <- tmp2
		tmp <- diag(2*abs(x)^(2*gamma) - abs(x - x)^(2*gamma)) + mat0 + t(mat0)
	}
	else{
		y <- as.numeric(y); m <- length(y)
		tmp <- matrix(NA, nc=n, nr=m)
		for(i in 1:m){
			for(j in 1:n){
				tmp[i,j] <- abs(y[i])^(2*gamma) + abs(x[j])^(2*gamma) - abs(y[i] - x[j])^(2*gamma)
			}
		}
	}
	
	tmp
}

