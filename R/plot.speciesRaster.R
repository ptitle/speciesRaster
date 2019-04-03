##' @title Plot speciesRaster
##'
##' @description Plot a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param log boolean; should the cell values be logged?
##' @param colorRampRange numeric vector of min and max value for scaling the color
##' 	ramp. Automatically inferred if set to \code{NULL}. This is relevant if multiple
##' 	plots are desired on the same scale. See \code{\link{getMultiMapRamp}}. 
##' @param legend boolean; should legend be included?
##' @param col either a vector of color names that will be interpolated, or a color ramp
##' 	function that takes an integer (see for example \code{\link{colorRampPalette}})
##'	@param includeWorldMap boolean; should a world map be plotted?
##' @param box boolean; should box be drawn around plot?
##' @param axes boolean; should axes be included?
##' @param location location of legend, if included. See \code{\link{addRasterLegend}}.
##' @param singleSpCol color for single-species cells. See details.
##' @param ... additional parameters will be passed to the \code{\link{addRasterLegend}} function.
##'
##'
##' @details If \code{x} is a metric as generated with \code{cellMetrics_speciesRaster} that returns 0 
##' 	for single-species cells, then those cells (that have a value of 0) will be plotted in gray (or any color
##' 	as specified with \code{singleSpCol}).
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
##' @rdname plot
##' @export

plot.speciesRaster <- function(x, log = FALSE, colorRampRange = NULL, legend = TRUE, col = c('blue', 'yellow', 'red'), includeWorldMap = TRUE, box=TRUE, axes=TRUE, location = 'right', singleSpCol = gray(0.9), ...) {
	
	if (!inherits(x, 'speciesRaster')) {
		stop('Object must be of class speciesRaster')
	}
	
	# if x is a speciesRaster that represents a metric that only makes sense for communities with multiple species, 
	# then single species cells have a value of zero, and we will plot those cells as gray.
	if (names(x[[1]]) %in% c('range', 'mean_NN_dist', 'min_NN_dist', 'variance', 'disparity', 'rangePCA', 'meanPatristic', 'patristicNN', 'phyloDisparity', 'PSV')) {
		# determine which cells have just 1 species
		singleSpCells <- singleSpCellIndex(x)
		x[[1]][singleSpCells] <- NA
		singleSpRas <- raster::raster(x[[1]])
		singleSpRas[] <- NA
		singleSpRas[singleSpCells] <- 1
		plotSingleCells <- TRUE
	} else {
		plotSingleCells <- FALSE
	}	
	
	if (!is.null(colorRampRange)) {
		if (class(colorRampRange) != 'numeric' | length(colorRampRange) != 2) {
			stop('colorRampRange must be a vector of length 2: min and max.')
		}
	} else {
		colorRampRange <- raster::cellStats(x[[1]], stat=range)
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
	} else {
		raster::plot(log(x[[1]]), col = colramp(100), box = box, axes = axes, legend = FALSE, zlim = colorRampRange)			
	}
	
	if (plotSingleCells) {
		raster::plot(singleSpRas, col = singleSpCol, box = FALSE, axes = FALSE, legend = FALSE, add = TRUE)
	}
	
	if (legend) {
		if (!log) {
			addRasterLegend(x[[1]], location = location, ramp = colramp, ncolors = 100, minmax = colorRampRange, ...)
		} else {
			addRasterLegend(log(x[[1]]), location = location, ramp = colramp, ncolors=100, minmax = colorRampRange, ...)
		}		
	}

	if (includeWorldMap) {
		# add map for context
		wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
		graphics::plot(wrld, add = TRUE, lwd = 0.5)
	}
}