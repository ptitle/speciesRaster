##' @title Create speciesRaster
##' 
##' @description This function takes a list of SpatialPolygons or a rasterStack and generates
##' both a richness raster and an associated list of species per cell, creating an object of class \code{speciesRaster}.
##' 
##' 

##' @param ranges Either a list of SpatialPolygons, SpatialPolygonsDataFrames, a rasterstack,
##' 	rasterbrick, or species by cell matrix. If raster objects, cell values can either be
##'		binary presence/absence, or probabilities.

##' @param resolution the size of the cells of the resulting raster (ignored if input is
##' 	raster)

##' @param resUnits if 'degrees', then raster will be unprojected, if 'meters' then raster
##' 	will be projected to equal area Behrmann projection (ignored if input is raster)

##' @param extent if 'auto', then the maximal extent of the polygons will be used, if 
##' 	input is raster, this will be ignored as the extent of the RasterStack will 
##' 	be used. If not auto, must be a numeric vector of length 4 with minLong, 
##'		maxLong, minLat, maxLat. 

##' @param coverCutoff only applies if input is SpatialPolygons. In rasterization of the
##' 	polygons, the percent of a raster cell that must be covered by the polygon for it to
##' 	count as a presence. 

##' @param nthreads number of threads to use for parallelization of the function. The R
##' 	package	\code{parallel} must be loaded for \code{nthreads > 1}.

##' @param rasterTemplate If input is a species x cell matrix, then a rasterTemplate must be
##' 	provided.
##' 
##' @param probRanking If \code{ranges} is a rasterStack of probabilities rather than
##' 	presence/absence, then apply a probability ranking filter to identify species.
##'		See details. 
##' 

##'	@details 
##' This function generates an object of class \code{speciesRaster}, which is a 
##'	list containing the following elements:
##'	\itemize{
##'		\item{raster} {A raster representing counts of species per cell.}
##'		\item{speciesList} {A list of species found in each cell.}
##'		\item{geogSpecies} {a vector of unique species in all cells.}
##'		\item{data} {An empty spot that morphological data can be added to.}
##'		\item{phylo} {An empty spot that a phylogeny can be added to.}
##'	}
##'	If input is a list of SpatialPolygons, then resolution must be specified.
##'	The extent will be inferred as the minimum extent required to encompass all 
##'	polygons. All cells within a polygon are considered as a presence for that species.
##'	
##'	If input is a rasterStack, then all parameters are taken from that, and resolution
##'	and extent arguments are ignored. If \code{probRanking = FALSE}, then any non-NA
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
##' # We will define our raster template where will specify extent, resolution and projection
##' listExtent <- getExtentOfList(tamiasPolyList)
##' template <- raster(ext = listExtent, res = c(20000, 20000), crs=proj4string(tamiasPolyList[[1]]))
##' 
##' # rasterize
##' ranges <- list()
##' for (i in 1:length(tamiasPolyList)) {
##' 	cat(i, '\n')
##' 	tmp <- rasterize(tamiasPolyList[[i]], template)
##' 	if (is.na(minValue(tmp))) {
##' 		# species' range is too small and gets dropped
##' 		# add back in
##' 		presenceCells <- unique(cellFromXY(tmp, spsample(tamiasPolyList[[i]], 10, type='random')))
##' 		tmp[presenceCells] <- 1
##' 	}
##' 	values(tmp)[!is.na(values(tmp))] <- 1
##' 	ranges[[i]] <- tmp
##' }
##' names(ranges) <- names(tamiasPolyList)
##' 
##' ranges <- stack(ranges)
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
##' assemblages predicting spatio-temporal patterns of species assemblages. J. Biogeog.
##' 38:1433-1444.
##'
##' D'Amen M., Dubuis A., Fernandes R.F., Pottier J., Pellissier L., Guisan A. 2015. Using 
##' species richness and functional traits predictions to constrain assemblage predictions
##' from stacked species distribution models. J. Biogeog. 42:1255-1266.
##'
##' @export



createSpeciesRaster <- function(ranges, resolution = 1, resUnits = 'degrees', extent = 'auto', coverCutoff = 0.5, nthreads = 1, rasterTemplate = NULL, probRanking = TRUE) {
	
	resUnits <- match.arg(resUnits, c('degrees', 'meters'))
	
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}

	if (class(ranges) == 'list') {
		if (class(ranges[[1]]) != 'SpatialPolygons' & class(ranges[[1]]) != 'SpatialPolygonsDataFrame') {
			stop('Input must be a list of SpatialPolygons or a RasterStack.')
		}
	}
	
	if (class(ranges) != 'matrix' & is.null(names(ranges))) {
		stop('List must be named.')
	}
	
	# prepare output object
	obj <- vector('list', length = 5)
	names(obj) <- c('raster', 'speciesList', 'geogSpecies', 'data', 'phylo')
	
	# if spatialpolygons
	if (class(ranges) == 'list') {
		if (class(ranges[[1]]) == 'RasterLayer') {
			stop('Rasters must be provided as a RasterStack.')
		}
		if (class(ranges[[1]]) == 'SpatialPolygons' | class(ranges[[1]]) == 'SpatialPolygonsDataFrame') {
			
			# test that all have same CRS
			if (length(unique(sapply(ranges, sp::proj4string))) != 1) {
				stop('proj4string of all polygons must match.')
			}

			if ('auto' %in% extent) {
				#get overall extent
				masterExtent <- getExtentOfList(ranges)
				masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
			} else if (is.numeric(extent) & length(extent) == 4) {
				masterExtent <- list(minLong = extent[1], maxLong = extent[2], minLat = extent[3], maxLat = extent[4])
			} else {
				stop("extent must be 'auto' or a vector with minLong, maxLong, minLat, maxLat.")
			}
			
			if (resUnits == 'degrees') {
				proj <- '+proj=longlat +datum=WGS84'
			} else {
				#Behrmann equal area projection
				proj <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
			}
			
			# if proj does not match CRS of polygons, transform
			if (proj != sp::proj4string(ranges[[1]])) {
				ranges <- lapply(ranges, function(x) sp::spTransform(x, sp::CRS(proj)))
			}
			
			#create template raster
			ras <- raster::raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(resolution, 2), crs = proj)
			
			# get the percent cover for each polygon, for each cell
			if (nthreads > 1) {
				cl <- parallel::makePSOCKcluster(nthreads)
				parallel::clusterExport(cl = cl, varlist = c('ranges', 'ras', 'rasterize'), envir = environment())
				polycover <- parallel::parLapply(cl, ranges, function(x) raster::rasterize(x, ras, getCover = TRUE))
				parallel::stopCluster(cl)
			} else {
				polycover <- lapply(ranges, function(x) raster::rasterize(x, ras, getCover = TRUE))
			}
			
			polycover <- lapply(polycover, function(x) x / 100)
			
			#identify cells for each polygon that satisfy cutoff
			cellInd <- lapply(polycover, function(x) which(raster::values(x) > coverCutoff))
			
			# create richness raster
			raster::values(ras) <- 0
			pb <- txtProgressBar(min = 0, max = length(cellInd), style = 3)
			for (i in 1:length(cellInd)) {
				setTxtProgressBar(pb, i)
				raster::values(ras)[cellInd[[i]]] <- raster::values(ras)[cellInd[[i]]] + 1
			}
			ras[ras == 0] <- NA
			names(ras) <- 'spRichness'
			obj[[1]] <- ras

			#create list of cells with species names
			spByCell <- vector('list', length = raster::ncell(ras))
			for (i in 1:length(cellInd)) {
				if (length(cellInd[[i]]) > 0) {
					for (j in 1:length(cellInd[[i]])) {
						spByCell[[cellInd[[i]][j]]] <- c(spByCell[[cellInd[[i]][j]]], names(cellInd)[i])
					}
				}
			}
			ind <- which(sapply(spByCell, is.null) == TRUE)
			spByCell[ind] <- lapply(spByCell[ind], function(x) vector('character', length = 0))
			obj[[2]] <- spByCell

		}
	}
	
	
	# if rasterstack as input
	if (class(ranges) %in% c('RasterStack', 'RasterBrick')) {
		
		#check that all rasters have values
		valCheck <- raster::minValue(ranges)
		badEntries <- which(is.na(valCheck))
		badEntriesRet <- badEntries
		if (length(badEntries) > 0) {
			badEntries <- paste(which(is.na(valCheck)), collapse = ', ')
			warning(paste0('The following rasters have no non-NA cells: ', badEntries, '.'))
			return(badEntriesRet)
		}
		
		# rasterstack calculations only
		# convert all rasters to presence/absence: create matrix of cells (rows) x raster (cols)
		# extract values for all cells, all species by presence coordinate
		presCoords <- raster::rasterToPoints(ranges)
		presCells <- raster::cellFromXY(ranges[[1]], presCoords[,1:2])
		presCoords <- presCoords[,3:ncol(presCoords)]
		
		# place in full matrix of cell x species
		fullmat <- matrix(nrow = raster::ncell(ranges), ncol = raster::nlayers(ranges))
		fullmat[presCells, ] <- presCoords
		colnames(fullmat) <- colnames(presCoords)
				
		# get count of species per cell
		fullmat <- ifelse(is.na(fullmat), 0, 1)
		if (identical(unique(as.numeric(fullmat)), c(0,1))) {
			probRanking <- FALSE
		}
		if (probRanking) {
			cellSums <- round(rowSums(fullmat, na.rm = TRUE))
		} else {
			cellSums <- rowSums(fullmat)
		}
				
		ras <- ranges[[1]]
		raster::values(ras) <- cellSums
		
		#remove zero cells
		ras[ras == 0] <- NA
		names(ras) <- 'spRichness'
		obj[[1]] <- ras		

		# get list of which species are found in each cell
		if (probRanking) {
			# fullmat[which(is.na(fullmat), arr.ind = TRUE)] <- 0
			spByCell <- returnTopIndices(fullmat, cellSums)
			spByCell <- lapply(spByCell, function(x) colnames(fullmat)[x])
		} else {
			spByCell <- spListPerCell(fullmat)		
		}
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


