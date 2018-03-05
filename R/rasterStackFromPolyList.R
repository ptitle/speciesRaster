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
##' 	If not auto, can be a SpatialPolygon object, in which case the resulting rasterStack
##'		will be cropped and masked with respect to the polygon, or a SpatialPoints object, 
##' 	from which an extent object will be generated, or a numeric vector of length 4 
##' 	with minLong, maxLong, minLat, maxLat. If 'interactive', then an interactive plot
##'		will appear in which the user can draw the desired polygon extent.
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
##'		If \code{dropEmptyRasters = TRUE} and \code{retainSmallRanges = TRUE}, then the species 
##' 	that will be dropped are those that are outside of the requested extent (which in that
##' 	case would be specified explicitly).  
##'
##'		In interactive mode for defining the extent, the user can draw a bounding polygon on a 
##'		map. A function call will then be printed to the console so that the user can hard-code 
##'		that bounding polygon in future calls to this function.
##' 
##' 	In interactive mode, the basemap is from \url{www.naturalearthdata.com}. 
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

	if (class(extent) == 'character') {
		
		#get overall extent
		masterExtent <- getExtentOfList(polyList)
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
		
	} else if (is.numeric(extent) & length(extent) == 4) {
		# use user-specified bounds
		masterExtent <- list(minLong = extent[1], maxLong = extent[2], minLat = extent[3], maxLat = extent[4])
		
	} else if (class(extent) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
		# get extent from points
		masterExtent <- raster::extent(extent)
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)

	} else if ("Extent" %in% class(extent)) {
		
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
		
	} else if (class(extent) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
		
		# use extent polygon
		if (!is.na(sp::proj4string(extent))) {
			if (!identical(sp::proj4string(extent), proj)) {
				stop('extent must have same projection as polyList.')
			}
		} else {
			sp::proj4string(extent) <- proj
		}
		
		masterExtent <- raster::extent(extent)
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
	} else {
		stop("extent must be 'auto', a SpatialPolygon or a vector with minLong, maxLong, minLat, maxLat.")
	}
	
	# interactive extent: if this option is selected, a coarse richness raster
	# will be plotted, so that the user can designate an extent polygon
	# at this point, masterExtent is a the maximal extent
	wkt <- NULL
	if (class(extent) == 'character') {
		if (extent == 'interactive') {
			
			# coarse template
			# if projected, use 100km, if not, use 20 degrees
			quickRes <- ifelse(sp::is.projected(polyList[[1]]), 100000, 20)
			quickTemplate <- raster::raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(quickRes, 2), crs = proj)
			quick <- lapply(polyList, function(x) raster::rasterize(x, quickTemplate))
			rich <- raster::calc(raster::stack(quick), fun=sum, na.rm = TRUE)
			
			# add map for context
			if (sp::is.projected(polyList[[1]])) {
				wrld <- sp::spTransform(worldmap, sp::CRS(proj))
			} else {
				wrld <- worldmap
			}
			
			plot(rich, legend = FALSE)
			cat('\n\tAn interactive coarse-grain map has been displayed.\n')
			cat('\n\tPlease wait until plot is completed......')
			
			plot(wrld, add = TRUE, lwd = 0.5)
			cat('done!\n')
			graphics::title(main = 'Define your extent polygon.')
		
			cat('\tClick on the map to create a polygon that will define the extent of the rasterStack.')
			cat('\tRight-clicking will close the polygon and terminate the interactive plot.\n\n')
			
			userPoly <- raster::drawPoly(sp = TRUE, col='red', xpd=NA)
			proj4string(userPoly) <- proj
			masterExtent <- raster::extent(userPoly)
			masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
			extent <- userPoly
			
			# display call so user can use this extent in the future
			wkt <- rgeos::writeWKT(userPoly)
			
			grDevices::dev.off()
		}
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
		
	# if extent was polygon, then mask rasterStack
	if (class(extent) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {

		ret <- raster::mask(ret, extent)
	}
	
	# if user wants to retain species that would otherwise be dropped
	# sample some random points in the range and identify cells
	valCheck <- raster::minValue(ret)
	smallSp <- which(is.na(valCheck))

	if (retainSmallRanges) {
		
		if (length(smallSp) > 0) {
			for (i in 1:length(smallSp)) {
				
				cover <- raster::rasterize(polyList[[smallSp[i]]], ras, getCover = TRUE)
				presenceCells <- which(raster::values(cover) > 0)

				if (length(presenceCells) > 0) {
					ret[[smallSp[i]]][presenceCells] <- 1
				}
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

	if (!is.null(wkt)) {
		cat('\n\tUse the same extent in the future by supplying the following to the extent argument:\n\n')
		cat(paste0('\trgeos::readWKT("', wkt, '")'), '\n\n')
	}

	return(ret)
}	

