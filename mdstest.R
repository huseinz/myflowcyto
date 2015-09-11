#!/usr/bin/Rscript

library("methods")
library("doParallel")
library("foreach")
library("flowCore")
library("MASS")


args <- commandArgs(TRUE)
filenames <- list.files("fcsFiles", pattern="*.fcs", full.names = TRUE)

print(filenames)
  
#cl <- makeCluster(30)
#registerDoParallel(cl)
#print(cl)

options(width = 200)

#sample 'size' points from data
size <- args[1]
if(is.na(size)){
	size <- 15
}

cat(size, file="outfile.txt", sep = "\n")

for(i in 26:26){
#foreach(i=1:30, .packages=c('flowCore', 'MASS'), .combine='c') %dopar% {

	tryCatch({
	
	#extract point data from file 
  	data <- read.FCS(filenames[i])
  	data <- exprs(data)

	#set values < 1 to 1
  	data[data<1] = 1
	
#	print("orig")
#	print(head(data, 10))
	#apply log to columns 
	data[, 3:12] = log10(data[, 3:12]) 
 	data <- unique(na.omit(data))

#	print("log")
#	print(head(data, 10))

	#normalize columns

#	apply(data, 2, function(x){x * 2})
	for(j in 1:ncol(data)){
		col_max = max(data[,j])
		print(col_max)
		data[,j] <- data[,j] * (1000 / col_max)
	}
# 	data <- unique(na.omit(data))

#	print("normalize")
#	print(head(data, 10))

	nrows <- dim(data)[1]
	ncols <- dim(data)[2]
	data <- array(data, dim = c(nrows, ncols))
	print(nrows)


	data <- data[sample(1:nrows, size, replace = FALSE), ]

	print(sprintf("Randomly selecting %s points", size))
	print("Printing first 15 from sample")
	print(head(data,15))
	
	#compute distance matrix of sampled points
	d <- dist(data, method = "manhattan")

	#write d to file
	cat(size, file=sprintf("dist/%s.2", i), sep="\n")	
	write.table(as.matrix(d), file=sprintf("dist/%s.2", i), append = TRUE, sep = " ", col.names = FALSE, row.names = FALSE)

if(FALSE){
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
	distortion_file = sprintf("distortion_mds/%s", i);

	cat("FCS File: ", file=distortion_file, sep = " ")
	cat(filenames[i], file=distortion_file, sep = "\n", append = TRUE)
	cat("", file=distortion_file, sep="\n", append = TRUE)
	cat("Number of points: ", file=distortion_file, sep = " ", append = TRUE)
	cat(size, file=distortion_file, sep = "\n", append = TRUE)
	cat("", file=distortion_file, sep="\n", append = TRUE)
	cat("Max Distortion: ", file=distortion_file, sep=" ", append = TRUE)	
	cat(max_distortion, file=distortion_file, sep="\n\n", append = TRUE)	
	cat("", file=distortion_file, sep="\n", append = TRUE)
	cat("Avg Distortion: ", file=distortion_file, sep=" ", append = TRUE)	
	cat(avg_distortion, file=distortion_file, sep="\n\n", append = TRUE)
	cat("", file=distortion_file, sep="\n", append = TRUE)
	cat("Projected Points", file=distortion_file, sep="\n", append = TRUE) 
	write.table(as.matrix(fit), file=distortion_file, sep=" ", append = TRUE, col.names = FALSE, row.names = FALSE)
	cat("", file=distortion_file, sep="\n", append = TRUE)
	cat("Original Points", file=distortion_file, sep="\n", append = TRUE) 
	write.table(data, file=distortion_file, sep=" ", append = TRUE, col.names = FALSE, row.names = FALSE)
	cat("", file=distortion_file, sep="\n", append = TRUE)
	cat("Distance Matrix of Original Points", file=distortion_file, sep="\n", append = TRUE)
	cat(size, file=distortion_file, sep="\n", append = TRUE)
	write.table(as.matrix(d), file=distortion_file, sep=" ", append = TRUE, col.names = FALSE, row.names = FALSE)

	#write table

	cat(max_distortion, file="outfile.txt", sep = "\n", append = TRUE)
	}
	},
	error = function(e) {print(e);}
	)

}
