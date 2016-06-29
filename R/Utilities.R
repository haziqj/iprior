### Creation of H.mat list
Hmat_list <- function (x, kernel, whichPearson, intr, no.int, gamma, xstar=list(NULL)) {
	p <- length(x)
	if (kernel == "FBM") H.mat <- mapply(fnH3, x=x[!whichPearson], y=xstar[!whichPearson], gamma=gamma, SIMPLIFY=F)
	else H.mat <- mapply(fnH2, x=x[!whichPearson], y=xstar[!whichPearson], SIMPLIFY=F)
	tmp <- mapply(fnH1, x=x[whichPearson], y=xstar[whichPearson], SIMPLIFY=F)
	H.mat <- c(H.mat, tmp)
	H.mat[c(which(!whichPearson), which(whichPearson))] <- H.mat
	if (!is.null(intr)) { #add in interactions, if any
		for(j in 1:no.int) {
			H.mat[[p+j]] <- H.mat[[ intr[1,j] ]] * H.mat[[ intr[2,j] ]]
			class(H.mat[[p+j]]) <- paste(class(H.mat[[ intr[1,j] ]]), class(H.mat[[ intr[2,j] ]]), sep=" x ")
		}
	}
	H.mat
}

### Indexer helper function
indx.fn <- function(k) { #indices for H.mat2
				 		 #note: intr, ind1 and ind2 are created in kernL()
	ind.int1 <- intr[1,]==k; ind.int2 <- intr[2,]==k	#locating var/kernel matrix
	ind.int <- which(ind.int1 | ind.int2)				#of interactions (out of 1:no.int)
	k.int <- ind.int+p	#which kernel matrix has interactions involves k
	k.int.lam <- c(intr[1,][ind.int2], intr[2,][ind.int1])	#which lambdas has interaction with k
	nok <- (1:p)[-k]	#all variables excluding k
	k.noint <- which(!(ind.int1 | ind.int2)) + p	#the opposite of k.int

	find.H2 <- function(z){ #this function finds position of H2
		x <- z[1]; y <- z[2]
		which( (ind1==x & ind2==y) | (ind2==x & ind1==y) )
	}

	### P.mat %*% R.mat + R.mat %*% P.mat indices
	za <- which((ind1 %in% k & ind2 %in% nok) | (ind2 %in% k & ind1 %in% nok))
	grid.PR <- expand.grid(k.int, nok)
	zb <- which(	(ind1 %in% grid.PR[,1] & ind2 %in% grid.PR[,2]) |
					(ind2 %in% grid.PR[,1] & ind1 %in% grid.PR[,2])
	)
	grid.PR.lam <- expand.grid(k.int.lam, nok)

	### P.mat %*% U.mat + U.mat %*% P.mat indices
	grid.PU1 <- expand.grid(k, k.noint)
	zc <- which(	(ind1 %in% grid.PU1[,1] & ind2 %in% grid.PU1[,2]) |
					(ind2 %in% grid.PU1[,1] & ind1 %in% grid.PU1[,2])
	)
	grid.PU2 <- expand.grid(k.int, k.noint)
	zd <- apply(grid.PU2, 1, find.H2)
	grid.PU.lam <- expand.grid(k.int.lam, k.noint)

	### P.mat %*% P.mat indices
	grid.Psq <- t(combn(c(k, k.int), 2))
	ze <- apply(grid.Psq, 1, find.H2)
	grid.Psq.lam <- NULL
	if(length(k.int.lam) > 0) grid.Psq.lam <- t(combn(c(0, k.int.lam), 2))

	list(	k.int=k.int, k.int.lam=k.int.lam,
			PRU=c(za,zc,zb,zd),
			PRU.lam1=c(	rep(0, length(nok)+length(k.noint)),
						grid.PR.lam[,1],
						grid.PU.lam[,1]	),
			PRU.lam2=c(nok, k.noint, grid.PR.lam[,2], grid.PU.lam[,2]),
			Psq=c(k, k.int), Psq.lam=k.int.lam,
			P2=ze, P2.lam1=grid.Psq.lam[,1], P2.lam2=grid.Psq.lam[,2]
	)
}

### Flatten function
# flatten <- function(x) {
	# len <- sum(rapply(x, function(x) 1L))
	# y <- vector("list", len)
	# i <- 0L
	# rapply(x, function(x) { i <<- i+1L; y[[i]] <<- x })
	# y
# }
#
###
### An internal function of ipriorEM() which does the pretty formatting in the report
###

decimalplaces <- function(x){	#this function calculates how many decimal places there are
  x <- as.numeric(x)
  if((x %% 1) != 0) nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  else return(0)
}

ipriorEMprettyIter <- function(x, iter){
  #first need to figure out the correct significant figures
  tmp <- prettyNum(x, digits=7, width=11)
  dig <- 7 + 10 - nchar(gsub(" ", "", tmp))
  for(i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11)


  #then trim down excess digits, because sometimes it overshoots to 11 width
  dig1 <- 7 + 10 - nchar(gsub(" ", "", tmp))
  ind <- dig1 < 7
  dig[ind] <- (dig + dig1 - 7)[ind]
  for(i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11)

  #sometimes it works, sometimes it doesn't. so run one more time.
  dig1 <- 7 + 10 - nchar(gsub(" ", "", tmp))
  ind <- dig1 < 7
  dig[ind] <- (dig + dig1 - 7)[ind]
  for(i in 1:length(dig)) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11)

  #after trimming to correct signif, just add the trailing zeroes.
  #makes use of function decimalplaces()
  miss <- 10 - nchar(gsub(" ", "", tmp))
  ind <- which(miss > 0)
  for(i in ind){
    if(!is.na(x[i])) tmp[i] <- prettyNum(x[i], digits=dig[i], width=11, nsmall=ifelse(decimalplaces(tmp[i])==0, miss[i]-1, decimalplaces(tmp[i])+miss[i]))
  }

  #what remains must be scientific numbers
  miss <- 10 - nchar(gsub(" ", "", tmp))
  ind <- which(miss > 0)
  for(i in ind){
    if(!is.na(x[i])) tmp[i] <- formatC(as.numeric(tmp[i]), format="e", digits=decimalplaces(tmp[i])-4+miss[i], width=11)
  }

  #then print result
  Iter <- format(paste0("Iteration " , iter, ":"), width=16, just="left")
  cat(Iter, tmp, "")
}

ipriorEMprettyLoglik <- function(x){
  tmp <- prettyNum(x, digits=8, width=10)
  miss <- 10 - nchar(gsub(" ", "", tmp))
  tmp <- format(x, digits=8, width=10, nsmall=decimalplaces(tmp)+miss)
  #one more time
  #dig <- 10 - nchar(gsub(" ", "", tmp))#; print(dig); print(decimalplaces(tmp))
  #tmp <- format(tmp, digits=8, width=10, nsmall=5)
  tmp
}

.onUnload <- function(libpath) {
  # Whenever you use C++ code in your package, you need to clean up after
  # yourself when your package is unloaded.
  library.dynam.unload("iprior", libpath)
}
