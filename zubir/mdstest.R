#!/usr/bin/Rscript

library("methods")
#library("doParallel")
#library("foreach")
library("flowCore")
library("MASS")


args <- commandArgs(TRUE)
filenames <- list.files("fcsFiles", pattern="*.fcs", full.names = TRUE)

print(filenames)
  
#cl <- makeCluster(30)
#registerDoParallel(cl)
#print(cl)

options(width = 200)


for(i in 1:1){
#foreach(i=1:30, .packages=c('flowCore', 'MASS'), .combine='c') %dopar% {

	tryCatch({
	
	#extract point data from file 
  	data <- read.FCS(filenames[i])
  	data <- exprs(data)
  	data[data<1] = 1
 	data <- unique(na.omit(data))
	nrows <- dim(data)[1]
	ncols <- dim(data)[2]
	data <- array(data, dim = c(nrows, ncols))
	print(nrows)

	#sample 'size' points from data
	size <- args[1]
	if(is.na(size)){
		size <- 15
	}
	data <- data[sample(1:nrows, size, replace = FALSE), ]

	print(sprintf("Randomly selecting %s points", size))
	print("Printing first 15 from sample")
	print(head(data,15))
	
	#compute distance matrix of sampled points
	d <- dist(data, method = "manhattan")

	#write d to file
	cat(size, file=sprintf("dist/%s", i), sep="\n")	
	write.table(as.matrix(d), file=sprintf("dist/%s", i), append = TRUE, sep = " ", col.names = FALSE, row.names = FALSE)

	#run MDS
	fit <- cmdscale(d, eig = FALSE, x.ret = FALSE, add = FALSE, k = 2) 

	#generate dist matrix from MDS projected points
	p <- dist(fit, method = "manhattan")

	#calculate average and max distortion
	diff <- abs(d - p)
	max_distortion <- max(abs(d - p))	
	avg_distortion <- mean(abs(d - p))	

	print("Max distortion")
	print(max_distortion)
	print("Avg distortion")
	print(avg_distortion)

	#write results to file
	cat("Max Distortion: ", file=sprintf("distortion_mds/%s", i), sep="\n")	
	cat(max_distortion, file=sprintf("distortion_mds/%s", i), sep="\n\n", append = TRUE)	
	cat("Avg Distortion: ", file=sprintf("distortion_mds/%s", i), sep="\n", append = TRUE)	
	cat(avg_distortion, file=sprintf("distortion_mds/%s", i), sep="\n\n", append = TRUE)
	cat("Projected Points", file=sprintf("distortion_mds/%s", i), sep="\n", append = TRUE) 
	cat(fit, file=sprintf("distortion_mds/%s", i), sep="\n\n", append = TRUE)
	write.table(data, file=sprintf("distortion_mds/%s", i), sep=" ", append = TRUE, col.names = FALSE, row.names = FALSE)
	cat("", file=sprintf("distortion_mds/%s", i), sep="\n", append = TRUE)
	cat("Distance Matrix of Original Points", file=sprintf("distortion_mds/%s", i), sep="\n", append = TRUE)
	cat(size, file=sprintf("distortion_mds/%s", i), sep="\n", append = TRUE)
	write.table(as.matrix(d), file=sprintf("distortion_mds/%s", i), sep=" ", append = TRUE, col.names = FALSE, row.names = FALSE)


	},
	error = function(e) {print(e);}
	)

}
