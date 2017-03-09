##' @title Plot speciesRaster
##'
##' @description Plot a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param log boolean; should the cell values be logged?
##' @param ... additional parameters will be passed to the raster package plot function.
##'
##' @return Nothing is returned. 
##'
##' @author Pascal Title
##' 
##' @examples
##' plot(tamiasSpRas)
##' 
##' plot(tamiasSpRas, legend=FALSE, axes=FALSE, box=FALSE)
##' addRasterLegend(tamiasSpRas[[1]], location = 'top', ramp=c('blue','yellow','red'))
##' 
##' @export

plot.speciesRaster <- function(x, log = FALSE, ...) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('Object must be of class speciesRaster')
	}
		
	colramp <- grDevices::colorRampPalette(c('blue','yellow','red'))
	if (!log) {
		raster::plot(x[[1]], col=colramp(100), ...)
	} else {
		raster::plot(log(x[[1]]), col=colramp(100), ...)
	}
}