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

### PEARSON KERNEL TEST
### n=100,1000,2000,...,10000
### with roughly 10% different factors. i.e. n=2000 has 200 unique factors
# N <- c(seq(100, 900, 100), seq(1000, 10000, 1000))
# res <- NULL
# for(n in N){
	# x <- sample(1:(n/10), n, T)
	# t1 <- system.time(fn.H1(x))[3]
	# t2 <- system.time(fn.H1a(x))[3]
	# res <- rbind(res, c(t1,t2))
# }
# res <- cbind(N, res)
# colnames(res) <- c("N", "old", "new")
# res
# save.image("kernel")

# # first plot
# plot(N[1:10], res[1:10,2], type="b", ylab="time (seconds)", xlab="sample size n", col=4, main="Timing comparison of old and new Pearson kernel functions")
# points(N[1:10], res[1:10,3], type="b", col=2)
# text(1000, 0.3, labels="new", col=2)
# text(965, 9.2, labels="old", col=4)

# # second plot
# plot(N[10:19], res[10:19,2], type="b", ylab="time (seconds)", xlab="sample size n", col=4, main="Timing comparison of old and new Pearson kernel functions")
# points(N[10:19], res[10:19,3], type="b", col=2)
# text(10000, 100, labels="new", col=2)
# text(9650, 3500, labels="old", col=4)


