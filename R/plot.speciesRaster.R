##' @title Plot speciesRaster
##'
##' @description Plot a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param log boolean; should the cell values be logged?
##' @param colorRampRange numeric vector of min and max value for scaling the color
##' 	ramp. Automatically inferred if set to \code{NULL}. This is relevant if multiple
##' 	plots are desired on the same scale. See \code{\link{getMultiMapRamp}}. 
##' @param includeLegend boolean; should legend be included?
##' @param col either a vector of color names that will be interpolated, or a color ramp
##' 	function that takes an integer (see for example \code{\link{colorRampPalette}})
##'	@param includeWorldMap boolean; should a world map be plotted?
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
##' # Example for how to plot multiple speciesRasters on the same color scale
##' # for illustration purposes, we will compare weighted endemism to
##' # phylogenetic weighted endemism
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' spRas1 <- cellMetrics_speciesRaster(tamiasSpRas, metric='weightedEndemism')
##' spRas2 <- cellMetrics_speciesRaster(tamiasSpRas, metric='phyloWeightedEndemism')
##' # get global min and max values
##' minmax <- getMultiMapRamp(list(spRas1, spRas2))
##' par(mfrow = c(1,2))
##' plot(spRas1, colorRampRange = log(minmax), log = TRUE, location='right')
##' plot(spRas2, colorRampRange = log(minmax), log = TRUE, location='left')
##' 
##' @export

plot.speciesRaster <- function(x, log = FALSE, colorRampRange = NULL, includeLegend = TRUE, col = c('blue', 'yellow', 'red'), includeWorldMap = TRUE, box=TRUE, axes=TRUE, location = 'right', ...) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('Object must be of class speciesRaster')
	}
	
	if (!is.null(colorRampRange)) {
		if (class(colorRampRange) != 'numeric' | length(colorRampRange) != 2) {
			stop('colorRampRange must be a vector of length 2: min and max.')
		}
	} else {
		colorRampRange <- c(raster::minValue(x[[1]]), raster::maxValue(x[[1]]))
		if (log) {
			colorRampRange <- log(colorRampRange)
		}	
	}
	
	if (class(col) == 'function') {
		colramp <- col
	} else {
		if (class(col) == 'character') {
			colramp <- grDevices::colorRampPalette(col)
		}
	}
	if (!log) {
		raster::plot(x[[1]], col = colramp(100), box = box, axes = axes, legend = FALSE, zlim = colorRampRange)
		if (includeLegend) {
			addRasterLegend(x[[1]], location = location, ramp = colramp, ncolors = 100, minmax = colorRampRange, ...)
		}
	} else {
		raster::plot(log(x[[1]]), col = colramp(100), box = box, axes = axes, legend = FALSE, zlim = colorRampRange)
		if (includeLegend) {
			addRasterLegend(log(x[[1]]), location = location, ramp = colramp, ncolors=100, minmax = colorRampRange, ...)
		}
	}
	
	if (includeWorldMap) {
		# add map for context
		wrld <- sp::spTransform(worldmap, sp::CRS(raster::projection(x[[1]])))
		plot(wrld, add = TRUE, lwd = 0.5)		
	}
}