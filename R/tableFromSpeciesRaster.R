##' @title Data table from speciesRaster
##'
##' @description Given one or several speciesRaster objects, create a table of 
##' 	values and associated coordinate data.
##'
##' @param ... objects of class \code{speciesRaster}, \code{RasterLayer} or \code{RasterStack}.
##' @param n number of cells to randomly subsample, no subsampling if \code{NULL}
##'	@param dropSingleSpCells logical; should cells with single species be excluded?
##'
##' @details A set of cells are identified in the speciesRaster objects. If \code{n=NULL},
##'		then all cells are used, otherwise cells are randomly subsampled. Values at those
##' 	cells are then returned. This table construction can be particularly useful for 
##' 	subsequent statistical analyses.
##'
##'		Only cells with data in all inputs are returned. If n is greater than the number
##' 	of cells with data, then fewer than n cells will be returned. 
##'
##' @return data.frame with input variables, as well as \code{"long"} and \code{"lat"}.
##'
##' @author Pascal Title
##'
##' @examples
##' 
##' tamiasSpRas
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' tamiasSpRas <- addTraits_speciesRaster(tamiasSpRas, tamiasTraits)
##' morphoDisp <- cellMetrics_speciesRaster(tamiasSpRas, metric='disparity')
##' meanPat <- cellMetrics_speciesRaster(tamiasSpRas, metric='meanPatristic')
##' 
##' tableFromSpeciesRaster(tamiasSpRas, morphoDisp, meanPat, n = 100, dropSingleSpCells = TRUE)
##' 
##' @export

tableFromSpeciesRaster <- function(..., n = NULL, dropSingleSpCells = TRUE) {

	x <- list(...)
	
	if (!all(sapply(x, inherits, c('speciesRaster', 'RasterLayer', 'RasterStack')))) {
		stop('All input objects must be either of class speciesRaster, RasterLayer or RasterStack')
	}
	
	# x is a list of speciesRaster or raster objects

	for (i in 1:length(x)) {
		if (inherits(x[[i]], 'speciesRaster')) {
			if (dropSingleSpCells) {
				singleCellInd <- singleSpCellIndex(x[[i]])
				x[[i]][[1]][singleCellInd] <- NA
			}
			x[[i]] <- x[[i]][[1]]
		} else if (inherits(x[[i]], 'RasterStack')) {
			x[[i]] <- raster::unstack(x[[i]])
		}
	}
	
	x <- unlist(x, recursive = FALSE)
	
	if (any(c('x', 'y') %in% sapply(x, names))) {
		stop('Column names x and y are reserved for the coordinate fields, but one of the inputs is named x or y.')
	}

	# now x is a list of rasters
	# to identify which cells are not NA for all supplied rasters, sum them.
	x <- raster::stack(x)
	rasterSum <- raster::calc(x, fun = sum)

	if (is.null(n)) {
		selectedCells <- 1:raster::ncell(rasterSum)
	} else {
		selectedCells <- raster::sampleRandom(rasterSum, na.rm = TRUE, size = n, cells = TRUE)[,1]
	}

	df <- as.data.frame(matrix(nrow = length(selectedCells), ncol = raster::nlayers(x) + 2))
	colnames(df) <- c(names(x), 'x', 'y')

	for (i in 1:raster::nlayers(x)) {
		df[,i] <- x[[i]][selectedCells]
	}

	df[, c('x', 'y')] <- raster::xyFromCell(rasterSum, cell = selectedCells)

	return(df)
}