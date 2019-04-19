##' @title Polygon List to rasterStack
##' 
##' @description Takes a list of polygons and creates a rasterStack.
##'
##'
##' @param polyList a list of polygons objects (sf or sp), named with taxon names
##'
##' @param resolution vertical and horizontal size of raster cell, in units 
##'		of the polygons' projection
##'
##' @param retainSmallRanges boolean; should small ranged species be dropped or preserved.
##'		See details.
##'
##' @param extent if 'auto', then the maximal extent of the polygons will be used. 
##' 	If not auto, can be a SpatialPolygon or sf object, in which case the resulting rasterStack
##'		will be cropped and masked with respect to the polygon, or a spatial coordinates object, 
##' 	from which an extent object will be generated, or a numeric vector of length 4 
##' 	with minLong, maxLong, minLat, maxLat. If 'interactive', then an interactive plot
##'		will appear in which the user can draw the desired polygon extent.
##'
##' @param dropEmptyRasters if \code{TRUE}, then species that have no presence cells will be dropped.
##'		If \code{FALSE}, then rasters will remain, filled entirely with \code{NA}. 
##'
##'
##' @details 
##' 	In the rasterization process, all cells for which the polygon covers the midpoint are
##' 	considered as present and receive a value of 1. If \code{retainSmallRanges = FALSE}, 
##' 	then species whose ranges are so small that no cell registers as present will be 
##' 	dropped. If \code{retainSmallRanges = TRUE}, then the cells that the small polygon
##' 	is found in will be considered as present, even if it's a small percent of the cell.
##' 
##'		If \code{dropEmptyRasters = TRUE} and \code{retainSmallRanges = TRUE}, then the species 
##' 	that will be dropped are those that are outside of the requested extent (which in that
##' 	case would be specified explicitly).  
##'
##'		In interactive mode for defining the extent, the user can draw a bounding polygon on a 
##'		map. The drawn polygon will then be printed to the console so that the user can hard-code 
##'		that bounding polygon in future calls to this function.
##'
##'		Any SpatialPolygon or SpatialPoints objects are converted to objects of class \code{sf}.
##'	
##'		This function uses the \code{fasterize} package for conversion from polygon to raster.
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
##' library(sf)
##' # example dataset: a list of 24 chipmunk distributions as polygons
##' head(tamiasPolyList)
##' 
##' rangeStack <- rasterStackFromPolyList(tamiasPolyList, resolution = 50000)
##' rangeStack
##' 
##' @export

rasterStackFromPolyList <- function(polyList, resolution = 50000, retainSmallRanges = TRUE, extent = 'auto', dropEmptyRasters = TRUE) {
	
	if (is.list(polyList)) {
		if (inherits(polyList[[1]], c('SpatialPolygons', 'SpatialPolygonsDataFrame'))) {
			# if class SpatialPolygons, convert to sf
			for (i in 1:length(polyList)) {
				polyList[[i]] <- sf::st_as_sf(polyList[[i]])
			}
		}
	} else if (!inherits(polyList[[1]], c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'sf', 'sfc'))) {
		stop('polyList must be a list of SpatialPolygons or Simple Features.')
	}

	if (is.null(names(polyList))) {
		stop('List must be named with species names.')
	}

	# test that all have same CRS
	if (!any(sapply(polyList, function(x) identical(sf::st_crs(polyList[[1]]), sf::st_crs(x))))) {
		stop('proj4string or EPSG of all polygons must match.')
	}
	proj <- sf::st_crs(polyList[[1]])
	
	# if WKT string, then convert to sf polygon
	if (inherits(extent, 'character')) {
		if (grepl('POLYGON \\(', extent)) {
			extent <- sf::st_sf(sf::st_as_sfc(extent))
		}
	}

	if (inherits(extent, 'character')) {
		
		if (!extent %in% c('auto', 'interactive')) {
			stop("If extent is a character vector, it can only be 'auto' or 'interactive'.")
		}
		
		if (extent == 'auto') {
		
			#get overall extent
			masterExtent <- getExtentOfList(polyList)
			masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
		}
	} else if (is.numeric(extent) & length(extent) == 4) {
		# use user-specified bounds
		masterExtent <- list(minLong = extent[1], maxLong = extent[2], minLat = extent[3], maxLat = extent[4])
		
	} else if (inherits(extent, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
		if (inherits(extent, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
			extent <- sf::st_as_sf(extent)
		}
			
		if (!is.na(sf::st_crs(extent))) {
			if (!identical(sf::st_crs(extent), proj)) {
				extent <- sf::st_transform(extent, crs = proj)
			}
		} else {
			sf::st_crs(extent) <- proj
		}
	  
	  if (inherits(extent, 'sfc')) {
	    extent <- sf::st_sf(extent)
	  }
		
		# get extent from spatial object
		masterExtent <- raster::extent(extent)
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)

	} else if (inherits(extent, 'Extent')) {
		
		masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
		
	} else {
		stop("extent must be 'auto', a spatial object or a vector with minLong, maxLong, minLat, maxLat.")
	}
	
	# interactive extent: if this option is selected, a coarse richness raster
	# will be plotted, so that the user can designate an extent polygon
	wkt <- NULL
	if (inherits(extent, 'character')) {
		if (extent == 'interactive') {

			interactive <- interactiveExtent(polyList)
			extent <- interactive$poly
			wkt <- interactive$wkt

			masterExtent <- raster::extent(extent)
			masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
		}
	}

	#create template raster
	ras <- raster::raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(resolution, 2), crs = proj$proj4string)
	
	rasList <- pbapply::pblapply(polyList, function(x) fasterize::fasterize(x, ras, fun = 'sum'))
		
	# force non-NA values to be 1
	for (i in 1:length(rasList)) {
		raster::values(rasList[[i]])[!is.na(raster::values(rasList[[i]]))] <- 1
	}
	
	ret <- raster::stack(rasList)
		
	# if extent was polygon, then mask rasterStack
	if (inherits(extent, c('sf', 'sfc'))) {
		ret <- raster::mask(ret, extent)
	}
	
	valCheck <- raster::minValue(ret)
	smallSp <- which(is.na(valCheck))

	# if user wants to retain species that would otherwise be dropped
	# we will use the getCover argument in raster::rasterize to find cells that are at all covered by polygon
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
		badEntries <- which(is.na(valCheck))
		badEntriesInd <- badEntries
		badEntries <- sort(names(ret)[badEntries])
		if (length(badEntries) > 0) {
			message('The following rasters have no non-NA cells:\n\n')
			for (i in 1:length(badEntries)) {
				message('\t', badEntries[i], '\n')
			}
			
			# drop empty rasters
			ret <- ret[[setdiff(1:raster::nlayers(ret), badEntriesInd)]]
		}		
	}

	if (!is.null(wkt)) {
		message('\n\tUse the same extent in the future by supplying the following string to the extent argument:\n\n')
		message(paste0('\t"', wkt, '"'), '\n\n')
	}

	return(ret)
}	
