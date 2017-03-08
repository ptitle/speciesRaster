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
##' 
##' @param nthreads number of threads for parallel processing
##' 
##' @param nreps Number of repetitions for Foote metric distribution.
##' 
##' @return object of class \code{speciesRaster} where the raster represents
##' 	calculations of the metric at every cell. 
##' 
##' @details 
##' 	Univariate morphological metrics
##' 	\itemize{
##' 		\item{mean}
##' 		\item{median}
##' 	}
##' 	Multivariate morphological metrics
##'		\itemize{
##'			\item{disparity} 
##' 		\item{range}
##' 		\item{rangePCA}
##' 		\item{Foote_pi}
##' 		\item{Foote_ri}
##' 		\item{Foote_di}
##' 		\item{Foote_meanDist}
##' 		\item{Foote_sdDist}
##' 		\item{Foote_maxDist}
##' 	}
##' 	Phylogenetic metrics
##' 	\itemize{
##'			\item{meanPatristic}
##'			\item{patristicNN} {nearest neighbor in patristic distance}
##'			\item{phyloDisparity} {sum of squared deviations in patristic distance}
##' 	}
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



cellMetrics_speciesRaster <- function(x, metric, var = NULL, nthreads = 1, nreps = 20) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	metric <- sapply(metric, match.arg, choices = c('mean', 'median', 'range', 'rangePCA', 'disparity', 'Foote_pi', 'Foote_ri', 'Foote_di', 'Foote_meanDist', 'Foote_sdDist', 'Foote_maxDist', 'meanPatristic', 'patristicNN', 'phyloDisparity'), USE.NAMES = FALSE)
	
	if (length(metric) > 1) {
		if (any(metric %in% c('mean','median','meanPatristic','patristicNN','phyloDisparity','disparity', 'range', 'rangePCA'))) {
			stop('You can only specify one metric, unless it is the Foote metrics.')
		}
	}
	
	# Prune species list according to metric of interest
	if (all(metric %in% c('mean', 'median', 'disparity', 'range', 'rangePCA', 'Foote_pi', 'Foote_ri', 'Foote_di', 'Foote_meanDist', 'Foote_sdDist', 'Foote_maxDist'))) {
		
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
		
	 } else if (all(metric %in% c('meanPatristic', 'patristicNN','phyloDisparity'))) {
	 	
	 	# check that there is a phylogeny in speciesRaster object
		if (is.null(x[['phylo']])) {
			stop('speciesRaster object does not contain a phylo object!')
		}

	 	# prune speciesRaster object down to species shared with phylogeny
		x[[2]] <- intersectList(x[[2]], x[['phylo']]$tip.label)
	} else {
		stop('Metric not recognized or bad combination!')
	}
	
	# create a mapping of which set of species are found in each cell
	uniqueComm <- unique(x[[2]])
	allComm <- sapply(x[[2]], function(y) paste(y, collapse='|'))
	uniqueCommLabels <- sapply(uniqueComm, function(y) paste(y, collapse='|'))

	cellMap <- lapply(uniqueCommLabels, function(x) which(allComm == x))
	names(cellMap) <- uniqueCommLabels	

	## ----------------------------------
	## MORPHOLOGY-RELATED METRICS
	
	if (all(metric %in% c('mean', 'median'))) {
		if (is.null(var) & class(x[['data']]) %in% c('matrix', 'data.frame')) {
			stop(paste0("For metric ", metric, ", a variable must be specified."))
		}
		if (!is.null(var) & class(x[['data']]) %in% c('matrix', 'data.frame')) {
			if (!var %in% colnames(x[['data']])) {
				stop("'var' not in data.")
			}
		}
		if (is.vector(x[[4]])) {
			trait <- x[[4]]
		} else {
			trait <- setNames(x[[4]][, var], rownames(x[[4]]))
		}
		resVal <- cellAvg(uniqueComm, trait = trait, stat = metric)
	}
	
	if (all(metric == 'disparity')) {
		# sum of the diagonal of the covariance matrix
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- sapply(uniqueComm[ind], function(y) sum(diag(cov(x[['data']][y,]))))
	}
	
	if (all(metric == 'range')) {
		# maximum of the distance matrix (0 if one sp)
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- sapply(uniqueComm[ind], function(y) max(dist(x[['data']][y, ])))
	}
	
	if (all(metric == 'rangePCA')) {
		pc <- prcomp(x[['data']])
		# retain 99% of the variation
		keep <- 1:which(cumsum(((pc$sdev)^2) / sum(pc$sdev^2)) >= 0.99)[1]
		if (length(keep) == 1) {keep <- 1:ncol(pc$x)}
		pc <- pc$x[,keep]
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- sapply(uniqueComm[ind], function(y) {
			sum(apply(pc[y,], 2, function(z) diff(range(z))))
		})	
	}
	
	if (all(metric %in% c('Foote_pi', 'Foote_ri', 'Foote_di', 'Foote_meanDist', 'Foote_sdDist', 'Foote_maxDist'))) {
		
		requestedMetrics <- c('pi','ri', 'di', 'mean_dist', 'sd_dist', 'max_dist')[c('Foote_pi', 'Foote_ri', 'Foote_di', 'Foote_meanDist', 'Foote_sdDist', 'Foote_maxDist') %in% metric]
		
		nnDistRes <- lapply(uniqueComm, function(y) nnDist(x[['data']][y, ], Nrep = nreps, requested = requestedMetrics))	
		
		resVal <- vector('list', length = length(requestedMetrics))
		
		for (i in 1:length(requestedMetrics)) {
			resVal[[i]] <- sapply(nnDistRes, function(y) {
				if (!anyNA(y)) {
					return(mean(y[, requestedMetrics[i]]))
				} else {
					return(NA)
				}
			})	
		}

		if (length(resVal) == 1) {
			resVal <- resVal[[1]]
		} else {
			names(resVal) <- requestedMetrics
		}
		
	}
	
	## ----------------------------------
	## PHYLOGENY-RELATED METRICS
	
	if (all(metric %in% c('meanPatristic', 'patristicNN', 'phyloDisparity'))) {
		
		# calculate pairwise patristic distance
		patdist <- cophenetic(x[['phylo']])
		diag(patdist) <- NA
		
		if (metric == 'meanPatristic') {
			# meanPatristic is 0 if 1 species, NA if no species
			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal <- sapply(uniqueComm, function(y) mean(unlist(patdist[y, y]), na.rm = TRUE))
			resVal[sapply(uniqueComm, length) == 1] <- 0
			resVal[sapply(uniqueComm, anyNA)] <- NA
		}
		
		if (all(metric == 'patristicNN')) {
			# the mean of the minimum patristic distance for each species present
			resVal <- sapply(uniqueComm, function(y) {
				if (!anyNA(y)) {
					if (length(y) > 1) {
						return(mean(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
					} else {
						return(0)
					}
				} else {
					return(NA)
				}
			})
		}

		if (all(metric == 'phyloDisparity')) {
			# the sum of the squared deviations from the mean
			# value of 0 if 1 species, NA if no species
			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal <- sapply(uniqueComm, function(y) sum((patdist[y, y] - mean(patdist[y,y], na.rm = TRUE)) ^ 2, na.rm = TRUE))
			resVal[sapply(uniqueComm, length) == 1] <- 0
			resVal[sapply(uniqueComm, anyNA)] <- NA
		}
		
	}
	
	## ----------------------------------
	## REMAP RESULTS TO RELEVANT CELLS AND ASSIGN TO RASTERS
	
	if (!is.list(resVal)) {
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
	obj <- vector('list', length = 5)
	names(obj) <- c('raster', 'speciesList', 'geogSpecies', 'data', 'phylo')
	obj[[1]] <- resRas
	obj[[2]] <- x[[2]]
	obj[[3]] <- sort(unique(unlist(x[[2]])))
	
	class(obj) <- 'speciesRaster'
	return(obj)	
}









