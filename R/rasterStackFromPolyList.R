##' @title Polygon List to rasterStack
##' 
##' @description Takes a list of polygons and creates a rasterStack.
##'
##'
##' @param polyList a list of SpatialPolygon objects, named with taxon names
##'
##' @param resolution vertical and horizontal size of raster cell, in units 
##'		of the polygons' projection
##'
##' @param retainSmallRanges boolean; should small ranged species be dropped or preserved.
##'		See details.
##'
##' @param extent if 'auto', then the maximal extent of the polygons will be used. 
##' 	If not auto, must be a numeric vector of length 4 with minLong, maxLong, minLat, maxLat.
##'
##' @param dropEmptyRasters if \code{TRUE}, then species that have no presence cells will be dropped.
##'		If \code{FALSE}, then rasters will remain, filled entirely with \code{NA}. 
##'
##' @param nthreads number of threads to use for parallelization of the function. 
##'
##' @details 
##' 	In the rasterization process, all cells for which the polygon covers the midpoint are
##' 	considered as present and receive a value of 1. If \code{retainSmallRanges = FALSE}, 
##' 	then species whose ranges are so small that no cell registers as present will be 
##' 	dropped. If \code{retainSmallRanges = TRUE}, then the cells that the small polygon
##' 	is found in will be considered as present.
##' 
##'		If \code{dropEmptyRasters = TRUE} and \code{retainSmallRanges = TRUE}, then the species that 
##' 	will be dropped are those that are outside of the requested extent (which in that case 
##' 	would be specified explicitly).  
##' 
##' @return an object of class \code{RasterStack} where all rasters contain values of 
##' either NA or 1. 
##' 
##' @author Pascal Title
##'
##' @examples
##' library(raster)
##' library(maptools)
##' # example dataset: a list of 24 chipmunk distributions as polygons
##' head(tamiasPolyList)
##' 
##' rangeStack <- rasterStackFromPolyList(tamiasPolyList, resolution = 50000)
##' rangeStack
##' 
##' @export

rasterStackFromPolyList <- function(polyList, resolution = 50000, retainSmallRanges = TRUE, extent = 'auto', dropEmptyRasters = TRUE, nthreads = 1) {
	
	if (class(polyList) == 'list') {
		if (!class(polyList[[1]]) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
			stop('Input must be a list of SpatialPolygons or a RasterStack.')
		}
	}

	if (is.null(names(polyList))) {
		stop('List must be named with species names.')
	}

	# test that all have same CRS
	if (length(unique(sapply(polyList, sp::proj4string))) != 1) {
		stop('proj4string of all polygons must match.')
	}
	proj <- sp::proj4string(polyList[[1]])

	if ('auto' %in% extent) {
		#get overall extent
		masterExtent <- getExtentOfList(polyList)
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
	} else if (is.numeric(extent) & length(extent) == 4) {
		masterExtent <- list(minLong = extent[1], maxLong = extent[2], minLat = extent[3], maxLat = extent[4])
	} else {
		stop("extent must be 'auto' or a vector with minLong, maxLong, minLat, maxLat.")
	}

	#create template raster
	ras <- raster::raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(resolution, 2), crs = proj)
	
	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads)
		parallel::clusterExport(cl = cl, varlist = c('polyList', 'ras', 'rasterize'), envir = environment())
		rasList <- pbapply::pblapply(polyList, function(x) raster::rasterize(x, ras), cl = cl)
		parallel::stopCluster(cl)
	} else {
		rasList <- pbapply::pblapply(polyList, function(x) raster::rasterize(x, ras))
	}
	
	# force non-NA values to be 1
	for (i in 1:length(rasList)) {
		raster::values(rasList[[i]])[!is.na(raster::values(rasList[[i]]))] <- 1
	}
	
	ret <- raster::stack(rasList)
	
	# if user wants to retain species that would otherwise be dropped
	# sample some random points in the range and identify cells
	valCheck <- raster::minValue(ret)
	smallSp <- which(is.na(valCheck))

	if (retainSmallRanges) {
		
		if (length(smallSp) > 0) {
			for (i in 1:length(smallSp)) {
				try(presenceCells <- unique(raster::cellFromXY(ret[[smallSp[i]]], sp::spsample(polyList[[smallSp[i]]], 10, type = 'random'))), silent = TRUE)
				if ('try-error' %in% class(presenceCells)) {
					counter <- 1
					while ('try-error' %in% class(presenceCells) & counter < 10) {
						try(presenceCells <- unique(raster::cellFromXY(ret[[smallSp[i]]], sp::spsample(polyList[[smallSp[i]]], 10, type = 'random'))), silent = TRUE)
					}
				}
				ret[[smallSp[i]]][presenceCells] <- 1
			}
		}
		
	} else {
		
		# drop those species with no presences (due to small range size)
		ret <- ret[[setdiff(1:raster::nlayers(ret), smallSp)]]
	}
	
	# if requested, drop rasters that are entirely NA (and print to screen for reference)
	if (dropEmptyRasters) {
		valCheck <- raster::minValue(ret)
		badEntries <- which(is.na(valCheck))
		badEntriesInd <- badEntries
		badEntries <- sort(names(ret)[badEntries])
		if (length(badEntries) > 0) {
			cat('The following rasters have no non-NA cells:\n\n')
			for (i in 1:length(badEntries)) {
				cat('\t', badEntries[i], '\n')
			}
			
			# drop empty rasters
			ret <- ret[[setdiff(1:raster::nlayers(ret), badEntriesInd)]]
		}		
	}

	return(ret)
}	

