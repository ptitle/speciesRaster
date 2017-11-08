##' @title Extract min and max for multiple speciesRasters
##'
##' @description Extracts the range of values across a list of speciesRaster
##' 	objects for use in plotting
##'
##' @param x list of objects of class \code{speciesRaster} or \code{RasterLayer}.
##'
##' @details If the user would like to plot multiple speciesRaster objects
##' with a standardized color ramp, then the returned values from this function
##' can be supplied to \code{\link{plot.speciesRaster}}. Also works with RasterLayer
##' objects. 
##'
##' @return a numeric vector of length 2: overall min and max value. 
##'
##' @author Pascal Title
##'  
##' @export


getMultiMapRamp <- function(x) {
	
	if (class(x) == 'list') {
		if (!all(unique(sapply(x, class)) %in% c('RasterLayer','speciesRaster'))) {
			stop('Input must be list of speciesRaster or RasterLayer objects.')
		}
		
		for (i in 1:length(x)) {
			if (class(x[[i]]) == 'speciesRaster') {
				x[[i]] <- x[[i]][[1]]
			}
		}
		
		return(range(c(sapply(x, raster::minValue), sapply(x, raster::maxValue))))
	}
	else {
		if (class(x) == 'RasterLayer') {
			return(c(raster::minValue(x), raster::maxValue(x)))
		}
		if (class(x) == 'speciesRaster') {
			return(c(raster::minValue(x[[1]]), raster::maxValue(x[[1]])))
		}
	}	
}
