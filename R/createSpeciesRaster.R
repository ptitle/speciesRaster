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
##' @param rasterTemplate If input is a species x cell matrix, then a rasterTemplate must be
##' 	provided.
##' 
##' @param probRanking If \code{ranges} is a rasterStack of probabilities rather than
##' 	presence/absence, then apply a probability ranking filter to identify species.
##'		See details. 
##'
##' @param verbose Primarily intended for debugging, print progress to the console.
##'
##' @param chunkSize number of raster rows to handle at one time. Inferred automatically 
##' 	if null.
##'
##' @param force temporary parameter that forces the function to take the more memory-
##' 	conservative approach, even if the dataset is small enough to handle otherwise.
##' 
##'	@details 
##' This function generates an object of class \code{speciesRaster}, which is a 
##'	list containing the following elements:
##'	\itemize{
##'		\item{\code{raster:}} {A raster representing counts of species per cell.}
##'		\item{\code{speciesList:}} {A list of species found in each cell.}
##'		\item{\code{geogSpecies:}} {a vector of unique species in all cells.}
##'		\item{\code{data:}} {An empty spot that morphological data can be added to.}
##'		\item{\code{phylo:}} {An empty spot that a phylogeny can be added to.}
##'	}
##'	
##'	If input is a RasterStack, then all parameters are taken from that, such as resolution,
##' extent and projection. If \code{probRanking = FALSE}, then any non-NA
##'	and non-zero cell is considered a presence. If \code{probRanking = TRUE}, then 
##'	cell values are evaluated as probabilities from 0 to 1, for example as output from 
##'	species distribution models (SDM). The probability ranking filter was described in 
##'	the SESAM framework (Guisan and Rahbek 2011, D'Amen et al. 2015) as a way to prevent
##'	overprediction of species richness from SDMs. If \code{probRanking = TRUE}, then the
##'	richness of a cell is calculated as the rounded sum of the probabilities of that cell,
##'	and then that number of species is preserved, selected via probability ranking. 
##'	
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
##' @references
##' Guisan A., Rahbek C. 2011. SESAM - a new framework integrating macroecological and
##' species distribution models for predicting spatio-temporal patterns of species 
##' assemblages. J. Biogeog.
##' 38:1433-1444.
##'
##' D'Amen M., Dubuis A., Fernandes R.F., Pottier J., Pellissier L., Guisan A. 2015. Using 
##' species richness and functional traits predictions to constrain assemblage predictions
##' from stacked species distribution models. J. Biogeog. 42:1255-1266.
##'
##' @export



createSpeciesRaster <- function(ranges, rasterTemplate = NULL, probRanking = FALSE, verbose = FALSE, chunkSize = NULL, force=FALSE) {
	
	if (!class(ranges) %in% c('RasterStack', 'RasterBrick', 'matrix')) {
		stop('Input must be a list of SpatialPolygons or a RasterStack.')
	}
		
	# prepare output object
	obj <- vector('list', length = 5)
	names(obj) <- c('raster', 'speciesList', 'geogSpecies', 'data', 'phylo')
	
	
	# if rasterstack as input
	if (class(ranges) %in% c('RasterStack', 'RasterBrick')) {
		
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
		spByCell <- vector('list', length = raster::ncell(ranges))
		
		if (verbose) cat('\t...Determining if rasterstack can be processed in memory...')
		if (raster::canProcessInMemory(ranges) & !force) {

			if (verbose) cat('yes\n')
			pBar <- FALSE
			mat <- matrix(nrow=raster::ncell(ranges), ncol=raster::nlayers(ranges))
			colnames(mat) <- names(ranges)

			for (i in 1:raster::nlayers(ranges)) {
				mat[, i] <- ranges[[i]][]
			}
		
			# set all NA to 0
			mat[is.na(mat)] <- 0
			
			# check if is binary
			if (identical(unique(as.numeric(mat)), c(0,1))) {
	 			probRanking <- FALSE
			} else {
				if (!probRanking) {
					# if not binary and probRanking is false, convert all to 0/1
					mat[mat != 0] <- 1
				}
			}
				
			# get count of species per cell
			cellSums <- rowSums(mat)
				
			if (probRanking) {
				cellSums <- round(cellSums)
			}
				
			# assign values to result raster
			raster::values(ras) <- cellSums
			
			# get list of which species are found in each cell
			# because we can't pass the entire data.table to rcpp code, break into chunks
			# of 100 records
			
			if (probRanking) {
				spByCell <- returnTopIndices(mat, cellSums)
				
				# convert indices to species names
				spByCell <- lapply(spByCell, function(x) colnames(mat)[x])
				
			} else {
				
				spByCell <- spListPerCell(mat)		
	
			}
			
		} else {
			
			# data too big. Split into subsets of rows
			if (verbose) cat('no\n')

			if (verbose) cat('\t...Determining appropriate chunk size...\n')
			if (is.null(chunkSize)) {			
				# testing chunksize
				# using code that is run for canProcessInMemory()
				test <- seq(5, nrow(ranges), by = 5)
				n <- 4 + raster::nlayers(ranges) - 1
				for (i in 1:length(test)) {
					check <- round(1.1 * ncol(ranges) * test[i]) * n  < 1e+06
					if (!check) {
						chunkSize <- test[i - 1]
						break
					}
	 			}
	 			if (length(chunkSize) < 1) {
	 				chunkSize <- min(test)
	 			}
	 		}		
			
			pBar <- TRUE
						
			nChunks <- round(nrow(ranges) / chunkSize)
			counter <- 1
			cellCounter <- 1
			
			if (verbose) cat('\t...Chunk size:', chunkSize, 'nChunks:', nChunks, '...\n')
			
			pb <- raster::pbCreate(nChunks, progress = 'text')

			for (i in 1:nChunks) {

				raster::pbStep(pb, step = i)
				
				start <- counter
				if (i == nChunks) {
					end <- nrow(ranges)
				} else {
					end <- start + chunkSize - 1
				}
				
				cellInd <- raster::cellFromRow(ranges[[1]], seq(start, end))
				
				mat <- matrix(nrow = length(cellInd), ncol = raster::nlayers(ranges))
				colnames(mat) <- names(ranges)
				
				for (j in 1:raster::nlayers(ranges)) {
	
					mat[, j] <- ranges[[j]][cellInd]
				}
							
				# set all NA to 0
				mat[is.na(mat)] <- 0
				
				if (!probRanking) {
					# if not binary and probRanking is false, convert all to 0/1
					mat[mat != 0] <- 1
				}
				
				# get count of species per cell
				cellSums <- rowSums(mat)
					
				if (probRanking) {
					cellSums <- round(cellSums)
				}
					
				# assign values to result raster
				raster::values(ras)[cellCounter : (cellCounter + nrow(mat) - 1)] <- cellSums
				
				# get list of which species are found in each cell
				# because we can't pass the entire data.table to rcpp code, break into chunks
				# of 100 records
				
				if (probRanking) {
					spByCell[cellCounter : (cellCounter + nrow(mat) - 1)] <- returnTopIndices(mat, cellSums)
					
					# convert indices to species names
					spByCell[cellCounter : (cellCounter + nrow(mat) - 1)] <- lapply(spByCell[cellCounter : (cellCounter + nrow(mat) - 1)], function(x) colnames(mat)[x])
					
				} else {
					
					spByCell[cellCounter : (cellCounter + nrow(mat) - 1)] <- spListPerCell(mat)		
		
				}
	
				counter <- counter + chunkSize
				cellCounter <- cellCounter + nrow(mat)
			}
		}
		if (pBar) {
			raster::pbClose(pb, timer = FALSE)
			cat('\n')
		}
			
		#remove zero cells
		ras[ras == 0] <- NA
		names(ras) <- 'spRichness'
		obj[[1]] <- ras		
		obj[[2]] <- spByCell
	}
	
	# input ranges can be a binary presence/absence sp x cell matrix
	# where rownames are species and columns are cells
	if (class(ranges) == 'matrix') {
		if (mode(ranges) == 'numeric') {
			if (identical(unique(as.vector(ranges)), c(0, 1))) {
				if (!is.null(rasterTemplate)) {
					if (raster::ncell(rasterTemplate) == ncol(ranges)) {
							
						raster::values(rasterTemplate) <- colSums(ranges)
						rasterTemplate[rasterTemplate == 0] <- NA
						
						obj[[1]] <- rasterTemplate
						obj[[2]] <- apply(ranges, 2, function(x) names(x[which(x == 1)]))
						
					} else {
						stop('If input is a species by cell matrix, then it must be numeric with 0 or 1 and a template raster must be provided with the same number of cells as there are columns in the species by cell matrix.')					
					}
				}
			}
		}		
	}
	
	if (class(obj[[1]]) != 'RasterLayer') {
		stop('Input type not supported.')
	}
	
	obj[[3]] <- sort(unique(names(ranges)))
	
	class(obj) <- 'speciesRaster'
	return(obj)	
}


