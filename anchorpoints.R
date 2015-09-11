#!/usr/bin/Rscript

library("methods")
#library("doParallel")
#library("foreach")
library("flowCore")
#library("MASS")


args <- commandArgs(TRUE)
filenames <- list.files("fcsFiles", pattern="*.fcs", full.names = TRUE)

print(filenames)
  
#cl <- makeCluster(30)
#registerDoParallel(cl)
#print(cl)

options(width = 200)

#points 
npoints <- args[1]
if(is.na(npoints)){
	npoints <- 15
}

for(i in 26:26){
#foreach(i=1:30, .packages=c('flowCore', 'MASS'), .combine='c') %dopar% {

	tryCatch({
	
	#extract point data from file 
  	data <- read.FCS(filenames[i])
  	data <- exprs(data)

	#set values < 1 to 1
  	data[data<1] = 1
	
	#apply log to columns 
	data[, 3:12] = log10(data[, 3:12]) 
 	data <- unique(na.omit(data))


	#normalize columns
	for(j in 1:ncol(data)){
		col_max = max(data[,j])
		print(col_max)
		if(col_max != 0){
			data[,j] <- data[,j] * (1000 / col_max)
		}
	}
# 	data <- unique(na.omit(data))

	nrows <- dim(data)[1]
	ncols <- dim(data)[2]
	data <- array(data, dim = c(nrows, ncols))
	print(nrows)
	
	#pick three random points to be anchor points
	anchors <- data[sample(1:nrows, 3, replace = FALSE), ]

	#write anchor points to file so we know what they are
	write.table(as.matrix(anchors), file="anchorpoints", sep=" ", col.names = FALSE, row.names = FALSE)

	#remove some points from the end so all matrices have 'npoints' rows
	#this should probably not be done
	if( (nrows %% npoints) != 0 ){
		data <- head(data, -(nrows %% npoints))
	}

	#split the data matrix into submatrices with # of rows being 'npoints'
	subsets <- split(as.data.frame(data), rep(1:(nrows / npoints), each = npoints))

	#prepend anchor points to each submatrix
	subsets <- lapply(subsets, function(x) { rbind(anchors, x)})

	#compute distance matrices and write them to file
	for(j in 1:length(subsets)){
	
		d <- dist(subsets[[i]], method = "manhattan")

		cat(npoints, file=sprintf("apoints/%s.%s", i, j), sep="\n")	
		write.table(as.matrix(d), file=sprintf("apoints/%s.%s", i, j), append = TRUE, sep = " ", col.names = FALSE, row.names = FALSE)

	}

	},
	error = function(e) {print(e);}
	)

}
