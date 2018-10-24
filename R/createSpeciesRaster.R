##' @title Create speciesRaster
##' 
##' @description This function takes a rasterStack and generates both a richness 
##' raster and an associated list of species per cell, creating an object of 
##' class \code{speciesRaster}.
##' 
##' 
##' @param ranges Either a RasterStack, RasterBrick, or species by cell matrix. If 
##' 	raster objects, cell values can either be binary presence/absence, or probabilities.
##'
##' @param rasterTemplate If input is a species x cell matrix, then a rasterTemplate 
##' must be provided. Cells with a value of 1 will be processed, cells with a value 
##' of 0 will be ommitted. Therefore, all cells must have a value of 0/1.
##'
##' @param verbose Primarily intended for debugging, print progress to the console.
##'
##' 
##'	@details 
##' This function generates an object of class \code{speciesRaster}, which is a 
##'	list containing the following elements:
##'	\itemize{
##'		\item{\code{raster:}} {A raster representing counts of species per cell.}
##'		\item{\code{speciesList:}} {A list of species found in each cell.}
##'		\item{\code{geogSpecies:}} {a vector of unique species in all cells.}
##' 	\item{\code{cellCount:}} {a vector of counts of presence cells for each species.}
##'		\item{\code{data:}} {An empty spot that morphological data can be added to.}
##'		\item{\code{phylo:}} {An empty spot that a phylogeny can be added to.}
##'	}
##'	
##'	If input is a RasterStack, then all parameters are taken from that, such as resolution,
##' extent and projection. Any non-NA and non-zero cell is considered a presence.	
##'	This function expects that all input rasters in the rasterStack have presence values
##'	(i.e., at least 1 non-NA value). If any rasters have exclusively NA cells, then the 
##'	function will stop with a warning, and the output will be the index in the rasterStack
##'	of those rasters.
##' 
##' @return an object of class \code{speciesRaster}. 
##' 
##' @author Pascal Title
##'
##' @examples
##' library(raster)
##' library(maptools)
##' # example dataset: a list of 24 chipmunk distributions as polygons
##' head(tamiasPolyList)
##' 
##' # convert polygon ranges to raster
##' ranges <- rasterStackFromPolyList(tamiasPolyList, resolution = 20000)
##' 
##' spRas <- createSpeciesRaster(ranges = ranges)
##' 
##' spRas
##' 
##' 
##'
##'
##' @export



createSpeciesRaster <- function(ranges, rasterTemplate = NULL, verbose = FALSE) {
	
	if (all(!class(ranges) %in% c('RasterStack', 'RasterBrick', 'matrix', 'data.frame'))) {
		stop('Input must be a list of SpatialPolygons or a RasterStack.')
	}
		
	# prepare output object
	obj <- vector('list', length = 7)
	names(obj) <- c('raster', 'speciesList', 'cellCommInd', 'geogSpecies', 'cellCount', 'data', 'phylo')
	
	
	# if rasterstack as input
	if (all(class(ranges) %in% c('RasterStack', 'RasterBrick'))) {
		
		#check that all rasters have values
		if (verbose) cat('\t...Checking for empty rasters...\n')
		valCheck <- raster::minValue(ranges)
		badEntries <- which(is.na(valCheck))
		badEntriesRet <- badEntries
		if (length(badEntries) > 0) {
			badEntries <- paste(which(is.na(valCheck)), collapse = ', ')
			warning(paste0('The following rasters have no non-NA cells: ', badEntries, '.'))
			return(badEntriesRet)
		}
		
		# rasterstack calculations only
		# create matrix of cells (rows) x raster (cols)
	
		# prepare result objects
		ras <- raster::raster(ranges[[1]])
		raster::values(ras) <- 0
		cellCommVec <- integer(length = raster::ncell(ranges))
		spByCell <- vector('list', length = raster::ncell(ranges))
		
		# determine the size of rasterStack that can be processed in memory	
		if (verbose) cat('\t...Determining if rasterstack can be processed in memory...')
		if (raster::canProcessInMemory(ranges)) {

			if (verbose) cat('yes\n')
			mat <- matrix(nrow=raster::ncell(ranges), ncol=raster::nlayers(ranges))
			colnames(mat) <- names(ranges)

			for (i in 1:raster::nlayers(ranges)) {
				mat[, i] <- ranges[[i]][]
			}
		
			# set all NA to 0
			mat[is.na(mat)] <- 0
			
			# check if is binary
			if (identical(unique(as.numeric(mat)), c(0,1))) {
				# if not binary and probRanking is false, convert all to 0/1
				mat[mat != 0] <- 1
			}
				
			# get count of species per cell
			cellSums <- rowSums(mat)
								
			# assign values to result raster
			raster::values(ras) <- cellSums
			
			# get list of which species are found in each cell
			spByCell <- spListPerCell(mat)		
			
		} else {
			
			# data too big. Split into subsets of rows
			if (verbose) cat('no\n')
			if (verbose) cat('\t...Determining how many rasters can be processed in memory...')	
			n <- 1
			while (raster::canProcessInMemory(ranges[[1:n]])) {
				n <- n + 1
			}
			
			if (verbose) cat(n, '\n')
					
			indList <- split(1:raster::nlayers(ranges), ceiling(1:raster::nlayers(ranges)/n))
			
			pb <- raster::pbCreate(length(indList), progress = 'text')	
	
			cellVals <- vector('list', length = length(indList))
			SpByCellList <- vector('list', length = length(indList))
			
			for (i in 1:length(indList)) {
									
				submat <- matrix(nrow=raster::ncell(ranges), ncol=length(indList[[i]]))
				colnames(submat) <- names(ranges)[indList[[i]]]

				for (j in 1:length(indList[[i]])) {
					submat[, j] <- ranges[[indList[[i]][j]]][]
				}
				
				# set all NA to 0
				submat[is.na(submat)] <- 0
				
				# check if is binary
				if (identical(unique(as.numeric(submat)), c(0,1))) {
					# if not binary and probRanking is false, convert all to 0/1
					submat[submat != 0] <- 1
				}
					
				# get count of species per cell
				cellSums <- rowSums(submat)
									
				# assign values to result raster
				cellVals[[i]] <- cellSums
				
				# get list of which species are found in each cell
				SpByCellList[[i]] <- spListPerCell(submat)
				
				raster::pbStep(pb, step = i)	
			}

			raster::pbClose(pb, timer = FALSE)
		
			# combine pieces
			if (verbose) cat('\t...Assembling speciesRaster...\n')	
			raster::values(ras) <- rowSums(do.call(cbind, cellVals))
			# for now, replace all NA with 'empty'
			for (i in 1:length(SpByCellList)) {
				for (j in 1:length(SpByCellList[[i]])) {
					if (all(is.na(SpByCellList[[i]][[j]]))) {
						SpByCellList[[i]][[j]] <- 'empty'
					}
				}
			}
			spByCell <- mergeLists(SpByCellList)
			spByCell <- lapply(spByCell, unique)
			spByCell[sapply(spByCell, length) == 0] <- NA

		}
		
		
		# reduce spByCell to unique communities and track
		if (verbose) cat('\t...Reducing species list to unique sets...')
		uniqueComm <- unique(spByCell)
		spByCell2 <- sapply(spByCell, function(y) paste(y, collapse = '|'))
		uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
		for (i in 1:length(uniqueComm2)) {
			cellCommVec[which(spByCell2 == uniqueComm2[i])] <- i
		}		
		if (verbose) cat('done\n')		
		
		#remove zero cells
		ras[ras == 0] <- NA
		names(ras) <- 'spRichness'
		obj[['raster']] <- ras		
		obj[['speciesList']] <- uniqueComm
		obj[['cellCommInd']] <- cellCommVec
		obj[['geogSpecies']] <- sort(unique(names(ranges)))
		
		# calculate range area for each species ( = number of cells)
		if (verbose) cat('\t...Calculating species cell counts...\n\n')
		obj[['cellCount']] <- countCells(convertNAtoEmpty(obj[['speciesList']]), obj[['geogSpecies']])
		names(obj[['cellCount']]) <- obj[['geogSpecies']]
	}
		
	# input ranges can be a binary presence/absence sp x cell matrix
	# where rownames are species and columns are cells
	if (any(class(ranges) %in%  c('matrix', 'data.frame'))) {
		if (any(class(ranges) == 'data.frame')) {
			ranges <- as.matrix(ranges)
		}
		if (length(unique(rownames(ranges))) != nrow(ranges)) {
			stop('rownames in species x cell matrix must be unique.')
		}
		if (mode(ranges) != 'numeric') {
			stop('matrix data does not appear to be numeric.')	
		}
		if (!identical(as.numeric(range(as.vector(ranges))), c(0, 1))) {
			mode(ranges) <- 'logical'
			mode(ranges) <- 'numeric'
		}
		if (is.null(rasterTemplate)) {
			stop('If input is a species x cell matrix, then a raster template must be provided.')
		}
		if (raster::ncell(rasterTemplate) != ncol(ranges)) {
			stop('If input is species x cell matrix, then number of columns must equal the number of raster cells.')
		}
		if (!identical(as.numeric(range(raster::values(rasterTemplate))), c(0, 1))) {
			stop('rasterTemplate can only have values of 0 or 1.')
		}
		if (verbose) cat('\t...Using species by cell matrix...\n')
		if (verbose) cat('\t...Calculating species richness...\n')
		dropCells <- which(raster::values(rasterTemplate) == 0)
		raster::values(rasterTemplate) <- colSums(ranges)
		if (length(dropCells) > 0) {
			rasterTemplate[dropCells] <- 0
		}
		rasterTemplate[rasterTemplate == 0] <- NA
		
		names(rasterTemplate) <- 'spRichness'
		
		if (verbose) cat('\t...Indexing species in cells...\n')
		obj[['raster']] <- rasterTemplate
		
		spByCell <- apply(ranges, 2, function(x) names(x[which(x == 1)]))
		emptyInd <- which(sapply(obj[['speciesList']], length) == 0)
		if (length(dropCells) > 0) {
			emptyInd <- union(emptyInd, dropCells)
		}
		emptyList <- rep(list(NA), length(emptyInd))
		spByCell[emptyInd] <- emptyList

		# reduce spByCell to unique communities and track
		uniqueComm <- unique(spByCell)
		spByCell2 <- sapply(spByCell, function(y) paste(y, collapse = '|'))
		uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
		
		cellCommVec <- integer(length = length(spByCell))
		
		for (i in 1:length(uniqueComm2)) {
			cellCommVec[which(spByCell2 == uniqueComm2[i])] <- i
		}
		
		obj[['speciesList']] <- uniqueComm
		
		obj[['cellCommInd']] <- cellCommVec
		
		obj[['geogSpecies']] <- sort(rownames(ranges))
		
		# calculate range area for each species ( = number of cells)
		if (verbose) cat('\t...Calculating species cell counts...\n\n')
		obj[['cellCount']] <- rowSums(ranges)
		 
	}
	
	if (class(obj[[1]]) != 'RasterLayer') {
		stop('Input type not supported.')
	}
	
	class(obj) <- 'speciesRaster'
	return(obj)	
}


