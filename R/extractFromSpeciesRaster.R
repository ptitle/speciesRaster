##' @title Extract from speciesRaster
##'
##' @description Return species from intersection between spatial points or polygons
##' and a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}

##' @param spatial coordinates as either a SpatialPoints object, or as a 
##' 	marix/dataframe with two columns, or as a SpatialPolygons object.

##' @param collapse boolean; if \code{TRUE}, then a vector of unique species 
##' 	is returned, pooled from all cells, if \code{FALSE}, then list is 
##' 	returned with species from every cell as intersected by \code{spatial}. 
##'
##' @return A vector of species if \code{collapse = TRUE}, or a list of species by
##' 	cell if \code{collapse = FALSE}. 
##'
##' @author Pascal Title
##'
##' @examples
##' library(rgeos)
##' library(sp)
##' # get the projection of the speciesRaster
##' proj <- summary(tamiasSpRas)$projection
##' # define some points
##' pts <- rbind(
##' 		c(-120.5, 38.82),
##' 		c(-84.02, 42.75),
##' 		c(-117.95, 55.53))
##' pts <- SpatialPoints(coords=pts, proj4string=CRS('+proj=longlat +datum=WGS84'))
##' pts <- spTransform(pts, CRS(proj))
##' pts <- coordinates(pts)
##' 	
##' extractFromSpeciesRaster(tamiasSpRas, pts)
##' 
##' box <- readWKT("POLYGON((-120 38, -119 39, -120 40, -120 38))", 
##'		p4s=CRS('+proj=longlat +datum=WGS84'))
##' # transform to the same coordinate system as the speciesRaster object
##' box <- spTransform(box, CRS(proj))
##' 
##' # returns each cell's contents
##' extractFromSpeciesRaster(tamiasSpRas, box, collapse=FALSE)
##' 
##' # collapses results to unique set of species
##' extractFromSpeciesRaster(tamiasSpRas, box, collapse=TRUE)
##' 
##' @export


extractFromSpeciesRaster <- function(x, spatial, collapse=TRUE) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	if (class(spatial) %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
		cells <- raster::cellFromPolygon(x[[1]], p = spatial)[[1]]
	} else if (class(spatial) %in% c('matrix', 'data.frame', 'SpatialPoints', 'SpatialPointsDataFrame')) {
		cells <- raster::cellFromXY(x[[1]], spatial)
	} else {
		stop('format of spatial not recognized.')
	}
	
	# extract species from speciesRaster
	res <- sapply(cells, function(y) x[[2]][[y]])
	names(res) <- paste0('cell', cells)
	
	if (collapse) {
		return(sort(unique(unlist(res))))
	} else {
		return(res)
	}
}