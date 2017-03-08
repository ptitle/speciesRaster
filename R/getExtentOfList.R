##' @title Get extent of list
##'
##' @description Given a list of SpatialPolygons, return an extent
##' object that encompasses all items. 
##'
##' @param shapes a list of SpatialPolygons
##'
##' @return An object of class \code{extent}. 
##'
##' @author Pascal Title
##' 
##' @examples
##' getExtentOfList(tamiasPolyList)
##' 
##' @export


getExtentOfList <- function(shapes) {	
	x <- lapply(shapes, function(x) sp::bbox(x))
	minLong <- min(sapply(x, function(x) x[1], simplify = TRUE))
	maxLong <- max(sapply(x, function(x) x[3], simplify = TRUE))
	minLat <- min(sapply(x, function(x) x[2], simplify = TRUE))
	maxLat <- max(sapply(x, function(x) x[4], simplify = TRUE))
	
	res <- raster::extent(shapes[[1]])
	res@xmin <- minLong
	res@xmax <- maxLong
	res@ymin <- minLat
	res@ymax <- maxLat
	
	return(res)
}	
