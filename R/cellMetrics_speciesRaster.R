##' @title Cell Metrics
##' 
##' @description Calculate various morphological and phylogenetic community
##' metrics for every cell in a \code{speciesRaster} object. 
##' 
##' @param x object of class \code{speciesRaster}
##' 
##' @param metric name of metric to use, see Details. 
##'
##' @param var If a univariate morphological metric is specified, and the 
##' 	data in \code{x} are multivariate, which trait should be used?
##' 	This can also be a vector of column indices used to subset the data.
##' 
##' @param nreps Number of repetitions for Foote metric distribution.
##'
##' @param verbose Intended primarily for debugging, prints progress to the console
##' 
##' @return object of class \code{speciesRaster} where the raster represents
##' 	calculations of the metric at every cell. 
##' 
##' @details 
##' 	Univariate morphological metrics
##' 	\itemize{
##' 		\item{mean}
##' 		\item{median}
##' 		\item{range}
##'			\item{mean_NN_dist:} {mean nearest neighbor distance}
##'			\item{min_NN_dist:} {minimum nearest neighbor distance}
##' 		\item{variance}
##' 		\item{arithmeticWeightedMean} (see below)
##' 		\item{geometricWeightedMean} (see below)
##' 	}
##' 	Multivariate morphological metrics
##'		\itemize{
##'			\item{disparity} 
##' 		\item{range}
##' 		\item{rangePCA}
##' 		\item{mean_NN_dist:} {mean nearest neighbor distance}
##'			\item{min_NN_dist:} {minimum nearest neighbor distance}
##' 	}
##' 	Phylogenetic metrics
##' 	\itemize{
##'			\item{meanPatristic}
##'			\item{patristicNN:} {mean nearest neighbor in patristic distance}
##'			\item{phyloDisparity:} {sum of squared deviations in patristic distance}
##' 	}
##' 	Range-weighted metrics
##' 	\itemize{
##'			\item{weightedEndemism}
##'			\item{correctedWeightedEndemism:} {Weighted endemism standardized by species richness}
##'			\item{phyloWeightedEndemism:}
##' 	}
##'
##'		If data slot contains a pairwise matrix, var is ignored.
##'		Weighted mean options are available where, for each cell, a weighting scheme 
##' 	(inverse of species range sizes) is applied such that small-ranged species are 
##' 	up-weighted, and broadly distributed species are down-weighted. 
##' 	This can be a useful way to lessen the influence of broadly distributed species 
##' 	in the geographic mapping of trait data. 
##'
##' @examples
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' tamiasSpRas <- addTraits_speciesRaster(tamiasSpRas, tamiasTraits)
##'
##' # univariate morphological example
##' x <- cellMetrics_speciesRaster(tamiasSpRas, metric='mean', var='V2')
##' plot(x) 
##'
##' # multivariate morphological
##' x <- cellMetrics_speciesRaster(tamiasSpRas, metric='disparity')
##' plot(x)
##' 
##' # phylogenetic metrics
##' x <- cellMetrics_speciesRaster(tamiasSpRas, metric='meanPatristic')
##' plot(x)
##' 
##' @export


cellMetrics_speciesRaster <- function(x, metric, var = NULL, nreps = 20, verbose = FALSE) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}

	if (length(metric) > 1) {
		stop('You can only specify one metric.')
	}
	
	metric <- match.arg(metric, choices = c('mean', 'median', 'range', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean', 'rangePCA', 'disparity', 'mean_NN_dist', 'min_NN_dist', 'meanPatristic', 'patristicNN', 'phyloDisparity', 'weightedEndemism', 'correctedWeightedEndemism', 'phyloWeightedEndemism'))
	
	pairwise <- FALSE
	
	# if data is pairwise matrix, then set flag appropriately
	if (class(x[['data']]) %in% c('matrix', 'data.frame')) {
		if (identical(rownames(x[['data']]), colnames(x[['data']]))) {
			if (verbose) cat('\t...detected pairwise distance matrix...\n') 
			var <- NULL
			pairwise <- TRUE
			# make the diagonal and lower triangle NA
			x[['data']][lower.tri(x[['data']], diag = TRUE)] <- NA
			if (all(is.na(x[['data']][upper.tri(x[['data']])]))) {
				stop('There are no values in the upper triangle of the pairwise matrix.')
			}
		}
	}
	
	# if var is defined but data are univariate, disable var. 
	if (is.vector(x[['data']]) & !is.null(var)) {
		var <- NULL
	}
		
	if (!is.null(var) & class(x[['data']]) %in% c('matrix', 'data.frame')) {
		if (is.character(var)) {
			if (!var %in% colnames(x[['data']])) {
				stop('var not a valid column name of the data.')
			}
		} else if (is.numeric(var)) {
			if (!all(var %in% 1:ncol(x[['data']]))) {
				stop('Requested data column indices are out of range.')
			}
		}
	}
	
	if (class(x[['data']]) %in% c('matrix', 'data.frame') | is.vector(data)) {
		if (class(x[['data']]) %in% c('matrix', 'data.frame')) {
			metricType <- 'multiVar'
			if (!is.null(var) & length(var) == 1) {
				metricType <- 'uniVar'
			}
		} else {
			metricType == 'uniVar'
		}
	} else {
		metricType <- 'none'
	}
	
	# if a subset of data columns are requested, subset the data table
	if (!is.null(var) & class(x[['data']]) %in% c('matrix', 'data.frame')) {
		if (length(var) > 1) {
			x[['data']] <- x[['data']][, var]
		} else {
			x[['data']] <- setNames(x[['data']][, var], rownames(x[['data']]))
		}
		var <- NULL
	}
	
	if (metric %in% c('mean', 'median', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean') & metricType == 'multiVar' & !pairwise) {
		stop('If a univariate metric is requested from a multivariate dataset, a column name must be provided as var.')
	}
	
	if (metric %in% c('rangePCA', 'disparity') & metricType == 'uniVar' & !pairwise) {
		stop('A multivariate metric cannot be applied to a univariate dataset.')
	}
	
	# Prune species list according to metric of interest
	if (metric %in% c('mean', 'median', 'disparity', 'range', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean', 'rangePCA', 'mean_NN_dist', 'min_NN_dist')) {
		
		if (verbose) cat('\t...dropping species that are not in trait data...\n')
		# check that there is data in speciesRaster object
		if (is.null(x[['data']])) {
			stop('speciesRaster object does not contain trait data!')
		}
		
		# prune speciesRaster object down to species shared with data
		if (is.vector(x[['data']])) {
			x[[2]] <- intersectList(x[[2]], names(x[['data']]))
		} else {
			x[[2]] <- intersectList(x[[2]], rownames(x[['data']]))
		}
		
	 } else if (all(metric %in% c('meanPatristic', 'patristicNN','phyloDisparity', 'phyloWeightedEndemism'))) {
	 	
	 	# check that there is a phylogeny in speciesRaster object
		if (is.null(x[['phylo']])) {
			stop('speciesRaster object does not contain a phylo object!')
		}
		
		if (verbose) cat('\t...dropping species that are not in phylo data...\n')
	 	# prune speciesRaster object down to species shared with phylogeny
		x[[2]] <- intersectList(x[[2]], x[['phylo']]$tip.label)
	
	} else if (!metric %in% c('weightedEndemism', 'correctedWeightedEndemism')) {
		stop('Metric not recognized!')
	} else {
		# set empty entries to NA
		emptyEntries <- which(lengths(x[[2]]) == 0)
		x[[2]][emptyEntries] <- rep(list(NA), length(emptyEntries))
	}
	
	
	# create a mapping of which set of species are found in each cell
	if (verbose) cat('\t...Creating mapping of species combinations to cells...\n')
	uniqueComm <- unique(x[[2]])
	allComm <- sapply(x[[2]], function(y) paste(y, collapse='|'))
	uniqueCommLabels <- sapply(uniqueComm, function(y) paste(y, collapse='|'))

	# cellMap <- lapply(uniqueCommLabels, function(x) which(allComm == x))
	cellMap <- mapComm(uniqueCommLabels, allComm)
	names(cellMap) <- uniqueCommLabels	

	## ----------------------------------
	## MORPHOLOGY-RELATED METRICS
	
	## UNIVARIATE
	
	if (metric %in% c('mean', 'median', 'variance', 'range', 'mean_NN_dist', 'min_NN_dist') & !pairwise & metricType == 'uniVar') {
		if (verbose) cat('\t...calculating univariate metric:', metric, '...\n')
		trait <- x[['data']]
		resVal <- cellAvg(uniqueComm, trait = trait, stat = metric)
	}
	
	if (metric %in% c('arithmeticWeightedMean', 'geometricWeightedMean') & metricType == 'uniVar') {
		if (verbose) cat('\t...calculating univariate metric:', metric, '...\n')
		trait <- x[['data']]
			
		# get inverse range sizes (1 / cell counts)
		inverseCellCount <- 1 / x[['cellCount']]
		
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 0 & sapply(uniqueComm, anyNA) == FALSE)

		if (metric == 'arithmeticWeightedMean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				(1 / sum(inverseCellCount[y])) * sum(inverseCellCount[y] * trait[y])
			})
		}
		if (metric == 'geometricWeightedMean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				exp(sum(inverseCellCount[y] * log(trait[y])) / sum(inverseCellCount[y]))
			})			
		}	
	}
	
	if (metric %in% c('mean', 'median', 'variance') & pairwise) {
	
		# if pairwise matrix
	
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)

		if (metric == 'mean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(unlist(x[['data']][y, y]), na.rm = TRUE))
		} else if (metric == 'median') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::median(unlist(x[['data']][y, y]), na.rm = TRUE))
		} else if (metric == 'variance') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::var(unlist(x[['data']][y, y]), na.rm = TRUE))
		}
	}

	## MULTIVARIATE
	
	if (metric == 'disparity' & metricType == 'multiVar') {
		if (verbose) cat('\t...calculating multivariate metric:', metric, '...\n')
		# sum of the diagonal of the covariance matrix
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum(diag(cov(x[['data']][y,]))))
	}
	
	if (metric == 'range' & metricType == 'multiVar') {
		if (verbose) cat('\t...calculating multivariate metric:', metric, '...\n')
		# maximum of the distance matrix (0 if one sp)
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) max(dist(x[['data']][y, ])))
	}
	
	if (metric == 'rangePCA' & metricType == 'multiVar') {
		if (verbose) cat('\t...calculating multivariate metric:', metric, '...\n')
		pc <- prcomp(x[['data']])
		# retain 99% of the variation
		keep <- 1:which(cumsum(((pc$sdev)^2) / sum(pc$sdev^2)) >= 0.99)[1]
		if (length(keep) == 1) { keep <- 1:ncol(pc$x) }
		pc <- pc$x[,keep]
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
			sum(apply(pc[y,], 2, function(z) diff(range(z))))
		})	
	}
	
	if (metric == 'mean_NN_dist' & metricType == 'multiVar') {
		if (verbose) cat('\t...calculating multivariate metric:', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(nnDist(x[['data']][y, ], Nrep = nreps)$mean_dist))
	}

	if (metric == 'min_NN_dist' & metricType == 'multiVar') {
		if (verbose) cat('\t...calculating multivariate metric:', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) min(nnDist(x[['data']][y, ], Nrep = nreps)$mean_dist))
	}
	
	## ----------------------------------
	## PHYLOGENY-RELATED METRICS
	
	if (metric %in% c('meanPatristic', 'patristicNN', 'phyloDisparity')) {
		
		# calculate pairwise patristic distance
		patdist <- cophenetic(x[['phylo']])
		diag(patdist) <- NA
		
		if (metric == 'meanPatristic') {
			if (verbose) cat('\t...calculating phylo metric:', metric, '...\n')
			# meanPatristic is 0 if 1 species, NA if no species
			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal <- pbapply::pbsapply(uniqueComm, function(y) mean(unlist(patdist[y, y]), na.rm = TRUE))
			resVal[sapply(uniqueComm, length) == 1] <- 0
			resVal[sapply(uniqueComm, anyNA)] <- NA
		}
		
		if (metric == 'patristicNN') {
			if (verbose) cat('\t...calculating phylo metric:', metric, '...\n')
			# the mean of the minimum patristic distance for each species present
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(sapply(uniqueComm, length) > 1)
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				return(mean(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
			})
		}

		if (metric == 'phyloDisparity') {
			if (verbose) cat('\t...calculating phylo metric:', metric, '...\n')
			# the sum of the squared deviations from the mean
			# value of 0 if 1 species, NA if no species
			
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(sapply(uniqueComm, length) > 1)

			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum((patdist[y, y] - mean(patdist[y,y], na.rm = TRUE)) ^ 2, na.rm = TRUE) / choose(length(y), 2))
		}		
	}
	
	## ----------------------------------
	## RANGE-WEIGHTED METRICS
	
	if (metric %in% c('weightedEndemism', 'correctedWeightedEndemism')) {
		if (verbose) cat('\t...calculating weighted endemism metric...\n')
		
		# get inverse range sizes (1 / cell counts)
		inverseCellCount <- 1 / x[['cellCount']]

		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)		
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum(inverseCellCount[y]))
	
		# if corrected metric, we will standardize by species richness
		if (metric == 'correctedWeightedEndemism') {
			resVal[ind] <- resVal[ind] / lengths(uniqueComm[ind])
		}
	}
	
	if (metric == 'phyloWeightedEndemism') {
		if (verbose) cat('\t...calculating phylogenetic weighted endemism...\n')
		spEdges <- getRootToTipEdges(x[['phylo']])
		if (!'edgeArea' %in% names(x)) {
			if (verbose) cat('\t...calculating branch-specific range sizes...\n')
			x[['edgeArea']] <- do.call(cbind, phyloBranchRanges(x[['phylo']], convertNAtoEmpty(x[['speciesList']]), spEdges))
		}
		tipIndVec <- sapply(x[['phylo']]$tip.label, function(z) which(x[['phylo']]$tip.label == z))
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)

		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
			commEdges <- unique(unlist(spEdges[tipIndVec[y]])) + 1
			sub <- x[['edgeArea']][commEdges,]
			if (class(sub) == 'numeric') {
				sub <- matrix(sub, nrow = 1)
			}
			sum(sub[,1] / sub[,2])
		})
	}
	
 	
	
	## ----------------------------------
	## REMAP RESULTS TO RELEVANT CELLS AND ASSIGN TO RASTERS
	
	if (!is.list(resVal)) {
		resVal[is.nan(resVal)] <- NA
		cellVec <- numeric(length = length(allComm))
		for (i in 1:length(resVal)) {
			cellVec[cellMap[[i]]] <- resVal[i]
		}
		resRas <- raster::raster(x[[1]])
		raster::values(resRas) <- cellVec
	} else {
		resRas <- replicate(length(resVal), raster::raster(x[[1]]))
		for (i in 1:length(resVal)) {
			cellVec <- numeric(length = length(allComm))
			for (j in 1:length(resVal[[i]])) {
				cellVec[cellMap[[j]]] <- resVal[[i]][j]
			}
			raster::values(resRas[[i]]) <- cellVec
		}
		resRas <- raster::stack(resRas)
		names(resRas) <- names(resVal)
	}
	names(resRas) <- metric
	
	# prepare output object
	obj <- x
	obj[[1]] <- resRas
	obj[['geogSpecies']] <- sort(unique(unlist(x[['geogSpecies']])))
	
	class(obj) <- 'speciesRaster'
	return(obj)	
}




# # function to calculate phylo branch specific range areas

# phyloBranchRanges <- function(x) {

	# if (!'speciesRaster' %in% class(x)) {
		# stop('x must be of class speciesRaster.')
	# }
	
	# # check that there is a phylogeny in speciesRaster object
	# if (is.null(x[['phylo']])) {
		# stop('speciesRaster object does not contain a phylo object!')
	# }
	
	# phylo <- x[['phylo']]
	
	# allLeaves <- phangorn::Descendants(phylo, phylo$edge[,2], type = 'tips')
	# allLeaves <- lapply(allLeaves, function(y) phylo$tip.label[y])
	
	# branchTable <- matrix(nrow = length(phylo$edge.length), ncol = 2)
	# colnames(branchTable) <- c('branchLength', 'branchArea')
	# branchTable[,1] <- phylo$edge.length
	
	# for (i in 1:length(phylo$edge.length)) {
		
		# inCell <- sapply(x[['speciesList']], function(y) any(allLeaves[[i]] %in% y))
		# branchTable[i, 2] <- sum(inCell)	
	# }

	# x[['edgeArea']] <- branchTable
	# return(x)	
# }



# getCommEdges2 <- function(phylo, comm, indList) {
	
	# tipInd <- sapply(comm, function(z) which(phylo$tip.label == z))
	# edges <- unique(unlist(indList[tipInd]))
	# edges <- edges + 1
	# # edgeColors <- ifelse(1:length(phylo$edge.length) %in% edges, 'red','black')
	# # edgeWidths <- ifelse(1:length(phylo$edge.length) %in% edges, 2, 1)
	# # plot(phylo, tip.color=ifelse(phylo$tip.label %in% comm, 'blue','gray'), root.edge=TRUE, edge.color = edgeColors, edge.width = edgeWidths)
	# # nodelabels(node=nodes, frame='circle', bg='red', cex=0.1)
	
	# return(edges)
# }




# # library(raster)
# library(maptools)
# library(Rcpp)
# # convert polygon ranges to raster
# ranges <- rasterStackFromPolyList(tamiasPolyList, resolution = 20000)
# x <- createSpeciesRaster(ranges = ranges)
# x <- addPhylo_speciesRaster(x, tamiasTree)


convertNAtoEmpty <- function(spCellList) {
	for (i in 1:length(spCellList)) {
		if (all(is.na(spCellList[[i]]))) {
			spCellList[[i]] <- 'empty'
		}
	}
	return(spCellList)
}


# # phylo2 <- x[['phylo']]
# class(phylo2) <- 'list'

# sourceCpp('~/Desktop/testing.cpp')

# tipEdges <- getRootToTipEdges(phylo2)

# q <- phyloBranchRanges(phylo2, spCellList, tipEdges)
# q <- do.call(cbind, q)



