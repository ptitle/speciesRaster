##' @title Plot speciesRaster
##'
##' @description Plot a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param log boolean; should the cell values be logged?
##' @param includeLegend boolean; should legend be included?
##' @param colorvec vector of color names that will be used for the color ramp
##' @param box boolean; should box be drawn around plot?
##' @param axes boolean; should axes be included?
##' @param location location of legend, if included. See \code{\link{addRasterLegend}}.
##' @param ... additional parameters will be passed to the \code{\link{addRasterLegend}} function.
##'
##' @return Nothing is returned. 
##'
##' @author Pascal Title
##' 
##' @examples
##' plot(tamiasSpRas)
##' 
##' plot(tamiasSpRas, legend=FALSE, axes=FALSE, box=FALSE)
##' addRasterLegend(tamiasSpRas, location = 'top', ramp=c('blue','yellow','red'))
##' 
##' @export

plot.speciesRaster <- function(x, log = FALSE, includeLegend = TRUE, colorvec = c('blue', 'yellow', 'red'), box=TRUE, axes=TRUE, location = 'right', ...) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('Object must be of class speciesRaster')
	}
	
	colramp <- grDevices::colorRampPalette(colorvec)
	
	if (grepl('beta', names(x[[1]]))) {
		raster::plot(x[[1]], col = colramp(100), breaks = seq(from = 0, to = 1, length.out = (100)), legend = FALSE, box = box, axes = axes)
		if (includeLegend) {
			addRasterLegend(x[[1]], location = location, ramp = colorvec, ncolors=100, minmax = c(0, 1), ...)
		}
	} else {
		if (!log) {
			raster::plot(x[[1]], col = colramp(100), box = box, axes = axes, legend = FALSE)
			addRasterLegend(x[[1]], location = location, ramp = colorvec, ncolors=100, ...)
		} else {
			raster::plot(log(x[[1]]), col = colramp(100), box = box, axes = axes, legend = FALSE)
			addRasterLegend(log(x[[1]]), location = location, ramp = colorvec, ncolors=100, ...)
		}
	}
}