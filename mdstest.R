#!/usr/bin/Rscript

library("methods")
library("doParallel")
library("foreach")
library("flowCore")
library("MASS")

args <- commandArgs(TRUE)
filenames <- list.files("fcsFiles", pattern="*.fcs", full.names = TRUE)

print(filenames)
  
cl <- makeCluster(30)
registerDoParallel(cl)

print(cl)
#for(i in 1:2){
out <- foreach(i=1:30, .packages=c('flowCore', 'MASS'), .combine='c') %dopar% {

	tryCatch({
	
	#extract point data from file 
  	data <- read.FCS(filenames[i])
  	data <- exprs(data)
#  	data[data<1] = 1
	
 	data <- unique(na.omit(data))
  	data <-log10(data)
  	data <-  unique(na.omit(data))
#  	idx <- apply(data, 1, function(x) all(is.finite(x)))
#  	data <- data[idx, ]
  
  	nrows <- dim(data)[1]
  	ncols <- dim(data)[2]
	
	print(nrows)
	print(ncols)
  	data <- array(data, dim = c(nrows, ncols))

  	# normalize data to 0-1000 range in columns 3+
  	for (k in 3:dim(data)[2]){
    		m <- min(data[,k])
    		if (max(data[,k]) == m){
      			print(sprintf("col %d max: %.2f min : %.2f",k, m, max(data[,k]))) 
      			data[,k] <- data[,k] / max(data[,k])
    		}
    		else{
      			data[,k] <- (data[,k] - m) / (max(data[,k]) - m)
			data[,k] <- data[,k] * 1000
    		}
  	}

	#generate dist matrix from first 100 points
	d <- dist(head(data, 100), method = "manhattan")
	#run cmdscale
	fit <- cmdscale(d, k = 2, eig = FALSE, x.ret = FALSE, add = FALSE)	
	write.matrix(fit, file=sprintf("log/%s", i), sep = " ")
	#generate dist matrix from projected points
	p <- dist(fit, method = "manhattan")
	#return max distortion
	distortion <- max(abs(d - p))	
	distortion
	},
	error = function(e) {print(e);}
	)

}
	print(out)
	write(out, file = "outfile", ncolumns = 1, append = FALSE)
