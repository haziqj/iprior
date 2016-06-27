## Canonical kernel function
fn.H2 <- function(x, y=NULL){ #takes in vector of covariates
	if(is.null(y)) y <- x
	tmp <- tcrossprod(y, x)
	tmp
}

## Centred Canonical kernel function
fn.H2a <- function(x, y=NULL){ #takes in vector of covariates
	x <- as.numeric(x)
	if(is.null(y)) y <- x
	else y <- as.numeric(y)
	xbar <- mean(x)
	tmp <- tcrossprod(y-xbar, x-xbar)
	class(tmp) <- "Canonical"
	tmp
}

### Pearson kernel
fn.H1 <- function(x, y=NULL){ #takes in vectors of type factors

	# fn.index <- function(k, tmp=list()){ #function to list row/col
  		# w <- tmp[[k]]
 		# cbind(
    	# row = rep(w, times=length(w):1 ) ,
    	# col = unlist(lapply(1:length(w), function(x) c(NA,w)[-(0:x)])))
	# }
	
	ytmp <- y
	if(is.null(ytmp)) y <- x
	## check if the inputs are non-factors, suggest not to use Pearson kernel
	if(any(!is.factor(x), !is.factor(y))) warning("Non-factor type vector used with Pearson kernel.", call.=F)
															
															#combine x and y, unfactorise them and work with numbers
	z <- unlist(list(x,y)); z <- as.numeric(z)				#simply doing c(x,y) messes with the factors
	x <- z[1:length(x)]; y <- z[(length(x)+1):length(z)]
	if(any(is.na(match(y,x)))) stop("The vector y contains elements not belonging to x.")
	prop <- table(x)/length(x)
	
	unqy <- sort(unique(y))
	unqx <- sort(unique(x))
	tmpx <- lapply(unqx, function (k) which(x == k))
	tmpy <- lapply(unqy, function (k) which(y == k))
	tmp <- lapply(1:length(unqy), function (k) expand.grid(tmpy[[k]], tmpx[[unqy[k]]]) )
	#side note: can avoid for loop below by combining the list tmp using
	#do.call(rbind, tmp) or the faster data.table option as.matrix(data.table::rbindlist(tmp))
	#but tests found that this is actually slower.
	
	mat <- matrix(-1, nr=length(y), nc=length(x))
	for (i in 1:length(unqy)) {
		mat[as.matrix(tmp[[i]])] <- 1/prop[unqy[i]] - 1
	}
	class(mat) <- "Pearson"
	mat
}

### FBM kernel
fn.H3 <- function(x, y=NULL, gamma=NULL){ #takes in vector of covariates
	if(is.null(gamma)) gamma <- 0.5
	x <- as.numeric(x)
	n <- length(x)
	
	if(is.null(y)){
		tmp <- matrix(0, n, n)
		index.mat <- upper.tri(tmp, diag=T)
		index <- which(index.mat, arr.ind=T)
		tmp2 <- abs(x[index[,1]])^(2*gamma) + abs(x[index[,2]])^(2*gamma) - abs(x[index[,1]]-x[index[,2]])^(2*gamma)
		tmp[index.mat] <- tmp2 
		tmp2 <- tmp; diag(tmp2) <- 0
		tmp <- tmp + t(tmp2) 
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

## Centred and scaled FBM kernel
fn.H3a <- function(x, y=NULL, gamma=NULL){ #takes in vector of covariates
	if(is.null(gamma)) gamma <- 0.5
	x <- as.numeric(x)
	n <- length(x)
	
	if(is.null(y)){
		A <- matrix(0, n, n)
		index.mat <- upper.tri(A)
		index <- which(index.mat, arr.ind=T)
		tmp2 <- abs(x[index[,1]]-x[index[,2]])^(2*gamma)
		A[index.mat] <- tmp2 
		A <- A + t(A)
		rvec <- apply(A, 1, sum)
		s <- sum(rvec)
		rvec1 <- tcrossprod(rvec, rep(1,n))
		tmp <- (A - rvec1/n - t(rvec1)/n + s/(n^2))/(-2)
	}
	else{
		y <- as.numeric(y); m <- length(y)	

		A <- matrix(0, n, n)
		index.mat <- upper.tri(A)
		index <- which(index.mat, arr.ind=T)
		tmp2 <- abs(x[index[,1]]-x[index[,2]])^(2*gamma)
		A[index.mat] <- tmp2 
		A <- A + t(A)		
		rvec <- apply(A, 1, sum)
		s <- sum(rvec)
		rvec1 <- tcrossprod(rep(1,m), rvec)
		
		B <- matrix(0, m, n)
		indexy <- expand.grid(1:m, 1:n)
		B[,] <- abs(y[indexy[,1]]-x[indexy[,2]])^(2*gamma)
		qvec <- apply(B, 1, sum)
		qvec1 <- tcrossprod(qvec, rep(1,n))
		
		tmp <- (B - qvec1/n - rvec1/n + s/(n^2))/(-2)
	}
	class(tmp) <- paste("FBM", gamma, sep=",")
	tmp
}

### 'vectorised' versions of kernel functions

fnH1 <- function(x, y=NULL){
	res <- 0
	if( (ncol(x) > 1) && !is.null(ncol(x)) ){
		if((ncol(x) != ncol(y)) && !is.null(y)) stop("New data is structurally unsimilar.")
		for(i in 1:ncol(x)) res <- res + fn.H1(x=x[,i], y=y[,i])
	}
	else res <- fn.H1(x, y)
	return(res)
}

fnH2 <- function(x, y=NULL){
	res <- 0
	if( (ncol(x) > 1) && !is.null(ncol(x)) ){
		if((ncol(x) != ncol(y)) && !is.null(y)) stop("New data is structurally unsimilar.")
		for(i in 1:ncol(x)) res <- res + fn.H2a(x=x[,i], y=y[,i])
	}
	else res <- fn.H2a(x, y)
	return(res)
}

fnH3 <- function(x, y=NULL, gamma){
	res <- 0
	if( (ncol(x) > 1) && !is.null(ncol(x)) ){
		if((ncol(x) != ncol(y)) && !is.null(y)) stop("New data is structurally unsimilar.")
		for(i in 1:ncol(x)) res <- res + fn.H3a(x=x[,i], y=y[,i], gamma)
	}
	else res <- fn.H3a(x, y, gamma)
	return(res)
}