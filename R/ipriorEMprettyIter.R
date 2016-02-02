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