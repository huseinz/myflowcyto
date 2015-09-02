#!/usr/bin/Rscript
#
# This R file:
# 1) asks for input file name, just input the name as "XXX.fcs", which should exist under folder inputFiles/ 
# 2) computes column-wise normalized centroids distances normalized to 0-1000
#    for inputFile named inputFiles/XXX.fcs, and writes distances to controlDist/XXX.csv
#    in the format
#       #rows #columns
#       followed by space seperated control points or control distances

library("methods")
library("flowCore")
library("cluster")

args <- commandArgs(TRUE)
#setwd("../")
wd <- getwd()
inputFilePath <- capture.output(cat(wd,'/fcsFiles'))
inputFilePath <- gsub(pattern = " ", replacement = "", x = inputFilePath, ignore.case = T)
print(inputFilePath)

csvFilePath <- capture.output(cat(wd,'/data'))
csvFilePath <- gsub(pattern = " ", replacement = "", x = csvFilePath, ignore.case = T)
print(csvFilePath)

  startTime <- proc.time()
  filename <- basename(args) #file.choose() #listOfFiles[l] #
  print(sprintf("Filename : %s",filename))
  dataFilename <- gsub(pattern = "fcs$", replacement = "1.txt", x = filename, ignore.case = T)  
  controlFilename <- gsub(pattern = "fcs$", replacement = "2.txt", x = filename, ignore.case = T)  
  
  print(dataFilename)  
  print(controlFilename)  

  setwd(inputFilePath)

  tryCatch({
  	sampdat <- read.FCS(filename)
	  setwd("../")
  	fdat <- exprs(sampdat)
  	fdat[fdat<1] = 1
 	fdat <- unique(na.omit(fdat))
  	logfdat <-log10(fdat)
  	logfdat <-  unique(na.omit(logfdat))
  	idx <- apply(logfdat, 1, function(x) all(is.finite(x)))
  	puredata <- logfdat[idx, ]
  	print(dim(puredata))
  
  	nrows <- dim(puredata)[1]
  	ncols <- dim(puredata)[2]
  	data <- array(puredata, dim = c(nrows, ncols))

  	setwd(csvFilePath)
#  	write.table(data, file=controlFilename,  sep=' ', append=FALSE, row.names=FALSE, col.names=FALSE)

  	# normalize data to 0-1 range in columns 3+
  	for (k in 3:dim(data)[2]){
    		m <- min(data[,k])
    		if (max(data[,k]) == m){
      			print(sprintf("col %d max: %.2f min : %.2f",k, m, max(data[,k]))) 
      			data[,k] <- data[,k] / max(data[,k])
    		}
    		else{
      			data[,k] <- (data[,k] - m) / (max(data[,k]) - m)
    		}
  	}
 
	d <- dist(data, method = "manhattan")
	fit <- cmdscale(d, eig = TRUE, k = 2)	
	#fit
#	numberOfControlPointsPerCluster <- 3
#  
#  	#kmeans cluster
#  	myKMeans <-  kmeans(data, optimumK, iter.max=30) #clara(data, optimumK)
# 	centroids <- myKMeans$centers
#  	clusters <-  myKMeans$cluster
#
#  	finalClusters <- clusters
#    
#  	numberAllControlPoints <- numberOfControlPointsPerCluster * optimumK
# 	outCentroids <- matrix(0, numberAllControlPoints, dim(data)[2])
#  
#  	# sub cluster each cluster to numberOfControlPointsPerCluster groups
#  	# use the first centroid of each subcluster as the out-centroid or all centroids
#  	#according to flag minimumControlPoints
#  	for (m in 1:dim(centroids)[1]){
#    
#    		clusterM <- which(clusters == m)
#    		rows <- length(clusterM)
#    		subData <- matrix(0, rows, dim(data)[2])
#    		for (xx in 1:rows){
#      			subData[xx, ] <- data[clusterM[xx],]
#    		}
#    
#    		subKMeans <-  kmeans(subData, numberOfControlPointsPerCluster, iter.max=30) #clara(subData, numberOfControlPointsPerCluster)
#    		subCenters <- subKMeans$centers
#    		subClusters <-  subKMeans$cluster
#
#    		finalClusters[clusters == m] <- subClusters + (m-1) * 3
#    
#    		index <- (m - 1) * numberOfControlPointsPerCluster
#    		for (p in 1:numberOfControlPointsPerCluster){
#       			outCentroids[index+p,] <- subCenters[p,]
#    		}
#  	}

  	setwd(csvFilePath)
  	write.table(data, file=dataFilename,  sep=' ', append=FALSE, row.names=FALSE, col.names=FALSE)

  	elapsedTime <- proc.time() - startTime
  	print(sprintf("%s %s %s %.2f %s","Processed ",filename," in ",elapsedTime["elapsed"]," seconds."))  
  },
error = function(e) { print(e);}
)

