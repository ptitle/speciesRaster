##' @title Occurrences to rasterStack
##' 
##' @description Takes a set of occurrences and rasterizes the coordinates for each taxon.
##'
##'
##' @param occ a table of species coordinates, a list of species-specific tables of coordinates, 
##'		a spatial coordinates object, or a species-specific list of spatial coordinate objects (sf or sp).
##'		If in table form, each coordinate pair must have an associated species name. If in list form, 
##' 	each element of the list must be named with the name of the species.
##'
##' @param resolution vertical and horizontal size of raster cell, in coordinate units.
##'
##' @param extent if 'auto', then the maximal extent of the coordinates will be used. 
##' 	If not auto, can be a SpatialPolygon or sf object, in which case the resulting rasterStack
##'		will be cropped and masked with respect to the polygon, or a spatial coordinates object, 
##' 	from which an extent object will be generated, or a numeric vector of length 4 
##' 	with minLong, maxLong, minLat, maxLat. 
##'
##' @param coordHeaders headers for longitude and latitude columns (x and y), only necessary if \code{occ} is a table.
##'
##' @param taxonHeader header for taxon labels, only necessary if \code{occ} is a table.
##'
##' @param crs coordinate reference system, only necessary if \code{occ} is a table. Otherwise, this
##' 	information is pulled from \code{occ}. EPSG:4326 indicates unprojected long/lat.
##'
##' @details 
##' 	If there are spaces in taxon names, those will be replaced with underscores.
##'
##'		Any SpatialPolygon or SpatialPoints objects are converted to objects of class \code{sf}.
##' 
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
##' # To illustrate usage, we will randomly sample some coordinates from each species polygon
##' # and generate a couple of alternative input formats
##' 
##' # list of sf spatial objects
##' spOccList <- lapply(tamiasPolyList, function(x) st_sample(x, size = 10, type= 'random'))
##' spStack <- rasterStackFromOccurrences(spOccList, resolution = 50000, crs = st_crs(tamiasPolyList[[1]]))
##'
##' # list of coordinate tables
##' spOccList2 <- lapply(spOccList, function(x) st_coordinates(x))
##' spStack <- rasterStackFromOccurrences(spOccList2, resolution = 50000, crs = st_crs(tamiasPolyList[[1]]))
##' 
##' # single table of coordinates
##' spOccList3 <- spOccList2
##' for (i in 1:length(spOccList3)) {
##' 	spOccList3[[i]] <- cbind.data.frame(taxon = names(spOccList3)[i], spOccList3[[i]])
##' 	colnames(spOccList3[[i]]) <- c('taxon', 'X', 'Y')
##' }
##' spOccList3 <- do.call(rbind, spOccList3)
##' rownames(spOccList3) <- NULL
##' spOccList3[, 'taxon'] <- as.character(spOccList3[, 'taxon'])
##' spStack <- rasterStackFromOccurrences(spOccList3, resolution = 50000, 
##' 	coordHeaders = c('X', 'Y'), crs = st_crs(spOccList[[1]]))
##'
##' # a single labeled spatial object
##' spOccList4 <- st_as_sf(spOccList3[, c('taxon', 'X', 'Y')], coords = c('X','Y'), 
##' 	crs = st_crs(spOccList[[1]])$proj4string)
##' spStack <- rasterStackFromOccurrences(spOccList4, resolution = 50000)
##' 	
##' 
##' @export


rasterStackFromOccurrences <- function(occ, resolution = 50000, extent = 'auto', coordHeaders = c('long','lat'), taxonHeader = 'taxon', crs = '+init=epsg:4326') {
	
	# detect format and convert. Target is a list of separate tables of coordinates, one per species.
	if (inherits(occ, 'list')) {
		
		# do all elements in the list have the same class?
		if (length(unique(sapply(occ, class, simplify = FALSE))) != 1) {
			stop('Not all elements in occ have the same class.')
		}
		
		if (inherits(occ[[1]], c('SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {

			if (inherits(occ[[1]], c('SpatialPoints', 'SpatialPointsDataFrame'))) {
				occ <- lapply(occ, sf::st_as_sf)
			}
			
			crs <- sf::st_crs(occ[[1]])$proj4string
			
			if (length(unique(sapply(occ, function(x) sf::st_crs(x)$proj4string, simplify = FALSE))) != 1) {
				stop('Not all elements in occ have the same projection.')
			}
			
			occ <- lapply(occ, sf::st_geometry)

			if (inherits(occ[[1]], 'sfc_POINT')) {
				occ <- lapply(occ, sf::st_coordinates)
			}
	
			if (inherits(occ[[1]], 'sfc_MULTIPOINT')) {
				occ <- lapply(occ, function(x) sf::st_coordinates(x)[, 1:2])
			}	
		}
		
		for (i in 1:length(occ)) {
			occ[[i]] <- cbind.data.frame(taxon = names(occ)[i], occ[[i]])
			colnames(occ[[i]]) <- c(taxonHeader, coordHeaders)
		}
	
		occ <- do.call(rbind, occ)
		rownames(occ) <- NULL
		occ[, taxonHeader] <- as.character(occ[, taxonHeader])			
			
	} else {
		
		if (inherits(occ, c('SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
			if (inherits(occ, c('SpatialPoints', 'SpatialPointsDataFrame'))) {
				occ <- sf::st_as_sf(occ)
			}
			
			crs <- sf::st_crs(occ)$proj4string
			
			headers <- colnames(occ)
			headers <- setdiff(headers, coordHeaders)
			occ <- as.data.frame(occ)		
			occ <- cbind.data.frame(occ[, setdiff(headers, 'geometry')], sf::st_coordinates(occ[, 'geometry']))
			colnames(occ) <- c(setdiff(headers, 'geometry'), coordHeaders)
		} else if (!inherits(occ, c('matrix', 'data.frame'))) {
			stop('Input format of occ not recognized.')
		}
	}
	

	# Now, regardless of input format, occ should now be a single table of coordinates and associated taxon names
	
	if (!all(coordHeaders %in% colnames(occ))) {
		stop('Coordinate headers not found in occ.')
	}
	if (!taxonHeader %in% colnames(occ)) {
		stop('taxon header not found in occ.')
	}
	
	occ <- as.data.frame(occ, stringsAsFactors = FALSE)
	occ[, coordHeaders[1]] <- as.numeric(occ[, coordHeaders[1]])
	occ[, coordHeaders[2]] <- as.numeric(occ[, coordHeaders[2]])

	occ <- occ[, c(taxonHeader, coordHeaders)]
	occ[, taxonHeader] <- gsub('^\\s+|\\s+$', '', occ[, taxonHeader])
	occ[, taxonHeader] <- gsub('\\s+', '_', occ[, taxonHeader])
	
	# are any records missing critical information?
	if (any(occ[, taxonHeader] == '' | is.na(occ[, taxonHeader]))) {
		ind <- which(occ[, taxonHeader] == '' | is.na(occ[, taxonHeader]))
		message('\t', length(ind), ' records were dropped because taxon was not specified.')
		occ <- occ[ - ind, ]
	}
	
	if (any(occ[, coordHeaders[1]] == '' | is.na(occ[, coordHeaders[1]]) | occ[, coordHeaders[2]] == '' | is.na(occ[, coordHeaders[2]]))) {
		ind <- which(occ[, coordHeaders[1]] == '' | is.na(occ[, coordHeaders[1]]) | occ[, coordHeaders[2]] == '' | is.na(occ[, coordHeaders[2]]))
		message('\t', length(ind), ' records were dropped because they lacked coordinates.')
		occ <- occ[ - ind, ]
	}
	
	# split into list of tables
	occList <- split(occ, occ[, taxonHeader])		
	message('\t', 'Detected ', length(occList), ' taxa.')
	
	proj <- sf::st_crs(crs)
		
		# extent
	if (inherits(extent, 'character')) {
		
		if (extent != 'auto') {
			stop("If extent is a character vector, it must be 'auto'.")
		}
		
		if (extent == 'auto') {
		
			#get overall extent
			masterExtent <- list(minLong = min(occ[, coordHeaders[1]]), maxLong = max(occ[, coordHeaders[1]]), minLat = min(occ[, coordHeaders[2]]), maxLat = max(occ[, coordHeaders[2]]))
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
			
	#create template raster
	ras <- raster::raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(resolution, 2), crs = proj$proj4string)

	# rasterize each taxon occurrences
	rasList <- vector('list', length(occList))
	names(rasList) <- names(occList)
	for (i in 1:length(occList)) {
		rasList[[i]] <- raster::rasterize(occList[[i]][, coordHeaders], ras, fun = 'count')
	}

	# force non-NA values to be 1
	for (i in 1:length(rasList)) {
		raster::values(rasList[[i]])[!is.na(raster::values(rasList[[i]]))] <- 1
	}
	
	ret <- raster::stack(rasList)
	
	if (raster::nlayers(ret) == 0) {
		stop('Empty raster layers.')
	}
		
	# if extent was polygon, then mask rasterStack
	if (inherits(extent, c('sf', 'sfc'))) {
		ret <- raster::mask(ret, extent)
	}
	
	valCheck <- raster::minValue(ret)
	if (any(is.na(valCheck) | valCheck != 1)) {
		stop('Rasters are not binary.')
	}

	return(ret)
}




