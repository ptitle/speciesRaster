##' @title Extract from speciesRaster
##'
##' @description Return species from intersection between spatial points or polygons
##' and a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}

##' @param spatial coordinates as either a spatial points object (sp or sf), 
##' 	a marix/dataframe with two columns, a numeric vector of c(long, lat), 
##'		or as a spatial polygon object (sp or sf).

##' @param collapse boolean; if \code{TRUE}, then a vector of unique species 
##' 	is returned, pooled from all cells, if \code{FALSE}, then list is 
##' 	returned with species from every cell as intersected by \code{spatial}.
##'
##' @details If \code{spatial} is a spatial object, it will be transformed to the same 
##' 	 projection as \code{x} if needed. If \code{spatial} is not a spatial object,
#' 		it is assumed to be in the same projection as \code{x}.
##'
##' @return A vector of species if \code{collapse = TRUE}, or a list of species by
##' 	cell if \code{collapse = FALSE}. 
##'
##' @author Pascal Title
##'
##' @examples
##' library(sf)
##' # get the projection of the speciesRaster
##' proj <- summary(tamiasSpRas)$crs
##' # define some points
##' pts <- rbind.data.frame(
##' 		c(-120.5, 38.82),
##' 		c(-84.02, 42.75),
##' 		c(-117.95, 55.53))
##' colnames(pts) <- c('x', 'y')
##' ptsSF <- st_as_sf(pts, coords = 1:2, crs = "+init=epsg:4326")
##' pts <- st_coordinates(st_transform(ptsSF, crs = proj))
##' 
##' # extract with table of coordinates
##' extractFromSpeciesRaster(tamiasSpRas, pts)
##' 
##' # extract with spatial points object
##' extractFromSpeciesRaster(tamiasSpRas, ptsSF)
##'
##' # extract with spatial polygon
##' hull <- st_convex_hull(st_union(ptsSF))
##' extractFromSpeciesRaster(tamiasSpRas, hull)
##' 
##' 
##' # returns each cell's contents
##' extractFromSpeciesRaster(tamiasSpRas, hull, collapse=FALSE)
##' 
##' # collapses results to unique set of species
##' extractFromSpeciesRaster(tamiasSpRas, hull, collapse=TRUE)
##' 
##' @export


extractFromSpeciesRaster <- function(x, spatial, collapse=TRUE) {
	
	if (!inherits(x, 'speciesRaster')) {
		stop('x must be of class speciesRaster.')
	}
	
	if (inherits(spatial, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
		if (inherits(spatial, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
			spatial <- sf::st_as_sf(spatial)
		}
		
		if (!is.na(sf::st_crs(spatial))) {
			if (!identical(sf::st_crs(spatial), sf::st_crs(x[[1]]))) {
				spatial <- sf::st_transform(spatial, crs = sf::st_crs(x[[1]]))
			}
		} else {
			stop('spatial must have projection information.')
		}
			
		spatial <- sf::st_geometry(spatial)
		
		if (inherits(spatial, 'sfc_POINT')) {
			spatial <- sf::st_coordinates(spatial)
		}

		if (inherits(spatial, 'sfc_MULTIPOINT')) {
			spatial <- sf::st_coordinates(spatial)[, 1:2]
		}

	}
	
	if (inherits(spatial, 'sfc_POLYGON')) {
		cells <- raster::cellFromPolygon(x[[1]], p = as(spatial, 'Spatial'))[[1]]
	} else if (inherits(spatial, c('matrix', 'data.frame', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
		cells <- raster::cellFromXY(x[[1]], spatial)
	} else if (is.numeric(spatial) & is.vector(spatial)) {
		cells <- raster::cellFromXY(x[[1]], matrix(spatial, ncol = 2))
	} else {
		stop('format of spatial not recognized.')
	}
	
	# extract species from speciesRaster
	res <- sapply(x[['cellCommInd']][cells], function(y) x[['speciesList']][[y]])
	names(res) <- paste0('cell', cells)
	
	if (collapse) {
		return(sort(unique(unlist(res))))
	} else {
		return(res)
	}
}