##' @title Plot speciesRaster
##'
##' @description Plot a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param log boolean; should the cell values be logged?
##' @param includeLegend boolean; should legend be included?
##' @param col either a vector of color names that will be interpolated, or a color ramp
##' 	function that takes an integer (see for example \code{\link{colorRampPalette}})
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
##' plot(tamiasSpRas, includeLegend=FALSE, axes=FALSE, box=FALSE)
##' addRasterLegend(tamiasSpRas, location = 'top', ramp=c('blue','yellow','red'))
##' 
##' @export

plot.speciesRaster <- function(x, log = FALSE, includeLegend = TRUE, col = c('blue', 'yellow', 'red'), box=TRUE, axes=TRUE, location = 'right', ...) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('Object must be of class speciesRaster')
	}
	
	if (class(col) == 'function') {
		colramp <- col
	} else {
		if (class(col) == 'character') {
			colramp <- grDevices::colorRampPalette(col)
		}
	}
	
	# if (grepl('beta', names(x[[1]]))) {
		# raster::plot(x[[1]], col = colramp(100), breaks = seq(from = 0, to = 1, length.out = (100)), legend = FALSE, box = box, axes = axes)
		# if (includeLegend) {
			# addRasterLegend(x[[1]], location = location, ramp = colorvec, ncolors=100, minmax = c(0, 1), ...)
		# }
	# } else {
		if (!log) {
			raster::plot(x[[1]], col = colramp(100), box = box, axes = axes, legend = FALSE)
			if (includeLegend) {
				addRasterLegend(x[[1]], location = location, ramp = colramp, ncolors=100, ...)
			}
		} else {
			raster::plot(log(x[[1]]), col = colramp(100), box = box, axes = axes, legend = FALSE)
			if (includeLegend) {
				addRasterLegend(log(x[[1]]), location = location, ramp = colramp, ncolors=100, ...)
			}
		}
	# }
}