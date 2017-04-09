commA <- c('spA', 'spB', 'spC', 'spD')

commB <- c('spC','spD','spE')

commC <- c('spA','spC','spF')

allComm <- list(commA, commB, commC)

x <- matrix(NA, nrow=3, ncol=6)
x[1,] <- c(1,1,1,1,0,0)
x[2,] <- c(0,0,1,1,1,0)
x[3,] <- c(1,0,1,0,0,1)


# get all pairwise minimums
nComm <- 3
minMat <- matrix(nrow = nComm, ncol = nComm)
for (i in 1:nComm) {
	for (j in 1:nComm) {
		if (i < j) {
			ij <- length(setdiff(allComm[[i]], allComm[[j]]))
			ji <- length(setdiff(allComm[[j]], allComm[[i]]))
			minMat[i,j] <- min(ij, ji)
		}
	}
}

maxMat <- matrix(nrow = nComm, ncol = nComm)
for (i in 1:nComm) {
	for (j in 1:nComm) {
		if (i < j) {
			ij <- length(setdiff(allComm[[i]], allComm[[j]]))
			ji <- length(setdiff(allComm[[j]], allComm[[i]]))
			maxMat[i,j] <- max(ij, ji)
		}
	}
}


# unique species
S_T <- length(unique(unlist(allComm)))
sumSi <- length(unlist(allComm))
Sdiff <- sumSi - S_T

# Jaccard dissimilarity for multiple sites
(sum(minMat, na.rm=TRUE) + sum(maxMat, na.rm=TRUE)) / (sum(Sdiff) + sum(minMat, na.rm=TRUE) + sum(maxMat, na.rm=TRUE))

# Jaccard dissimilarity turnover component
(2 * sum(minMat, na.rm=TRUE)) / (sum(Sdiff) + 2 * sum(minMat, na.rm=TRUE))

# Jaccard dissimilarity nestedness component
((sum(maxMat, na.rm=TRUE) - sum(minMat, na.rm=TRUE)) / (sum(Sdiff) + sum(minMat, na.rm=TRUE) + sum(maxMat, na.rm=TRUE))) * (sum(Sdiff) / (sum(Sdiff) + 2 * sum(minMat, na.rm=TRUE)))



multiPrepR <- function(allComm) {
	
	# Si = total number of species in site i
	# S_T = total number of species in all sites
	# b_ij = number of species exclusive to site i relative to j
	# b_ji = number of species exclusive to site j relative to i

	
	nComm <- length(allComm)
	
	# get all pairwise minimums
	minMat <- matrix(nrow = nComm, ncol = nComm)
	for (i in 1:nComm) {
		for (j in 1:nComm) {
			if (i < j) {
				ij <- length(setdiff(allComm[[i]], allComm[[j]]))
				ji <- length(setdiff(allComm[[j]], allComm[[i]]))
				minMat[i,j] <- min(ij, ji)
			}
		}
	}
	
	# get all pairwise maximums
	maxMat <- matrix(nrow = nComm, ncol = nComm)
	for (i in 1:nComm) {
		for (j in 1:nComm) {
			if (i < j) {
				ij <- length(setdiff(allComm[[i]], allComm[[j]]))
				ji <- length(setdiff(allComm[[j]], allComm[[i]]))
				maxMat[i,j] <- max(ij, ji)
			}
		}
	}
	
	# unique species
	S_T <- length(unique(unlist(allComm)))
	sumSi <- length(unlist(allComm))
	Sdiff <- sumSi - S_T

	sumMin <- sum(minMat, na.rm=TRUE)
	sumMax <- sum(maxMat, na.rm=TRUE)
	
	res <- c(sumMin, sumMax, Sdiff)
	
	return(res)
}


# Jaccard dissimilarity for multiple sites
jaccardDissimilarityMulti <- function(vec) {
	return((vec[1] + vec[2]) / (vec[3] + vec[1] + vec[2]))
}

# Jaccard dissimilarity turnover component
jaccardTurnoverMulti <- function(vec) {
	return((2 * vec[1]) / (vec[3] + 2 * vec[1]))
}

# Jaccard dissimilarity nestedness component
jaccardNestednessMulti <- function(vec) {
	return(((vec[2] - vec[1]) / (vec[3] + vec[1] + vec[2])) * (vec[3] / (vec[3] + 2 * vec[1])))
}

jaccardDissimilarityMulti(multiPrep(allComm))
jaccardTurnoverMulti(multiPrep(allComm))
jaccardNestednessMulti(multiPrep(allComm))



require(Rcpp)
sourceCpp('~/Dropbox/speciesRaster/src/speciesRaster_cppFunctions.cpp')


x <- readRDS('~/spByCell.rds')
nbList <- readRDS('~/nbList.rds')



listToTable <- function(q) {
	
	nSites <- length(q)
	allsp <- unique(unlist(q))
	
	mat <- matrix(nrow=nSites, ncol=length(allsp))
	colnames(mat) <- allsp

	for (i in 1:nSites) {
		mat[i,] <- as.numeric(allsp %in% q[[i]])
	}	
	return(mat)
}

multiSite <- function(spByCell, nbList, metric) {
	
	out <- numeric(length(spByCell))
	
	emptyInd <- ListIsEmptyR(spByCell)
	emptyCells <- emptyInd[[1]] + 1
	nonEmptyCells <- emptyInd[[2]] + 1
	
	for (i in 1:length(nonEmptyCells)) {
		cat(i, '\n')
				
		nbCells <- nbList[[nonEmptyCells[i]]]
		allInd <- ListIsEmptyR(spByCell[nbCells])
		ind <- allInd[[2]] + 1
		goodCells <- nbCells[ind]
		
		subList <- spByCell[goodCells]
		
		prepVec <- multiPrep(subList)
		out[nonEmptyCells[i]] <- jaccardDissimilarityMulti(prepVec)
	
	}
	
	out[emptyCells] <- NA
	
	return(out)
}

## Rcpp took 186.184 elapsed

system.time(cellBeta <- calcBetaMultiSite(x, nbList, metric = 'jaccardDissimilarityMulti'))
system.time(out <- multiSite(spByCell, nbList))

beta.multi(listToTable(subList), 'jaccard')
