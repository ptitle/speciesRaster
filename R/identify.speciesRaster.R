##' @title Identify species in speciesRaster
##'
##' @description Returns the species found in a cell that is clicked on
##' in a plotted speciesRaster object. 
##'
##' @param x object of class \code{speciesRaster}
##' @param ... further arguments passed to \code{\link{identify}}
##' @param returnCell boolean; if FALSE, then species names are returned,
##' 	if TRUE, then cell numbers are returned. 
##'
##' @details This is a wrapper function for the \code{identify} function
##' in base graphics. This is primarily intended as a useful function for 
##' data exploration and spot-checking. 
##' 
##' @return A vector of species names or cell numbers. 
##'
##' @author Pascal Title
##' 
##' @rdname identify
##' @export

identify.speciesRaster <- function(x, ..., returnCell = FALSE) {
	
	if (raster::ncell(x[[1]]) != length(x[['cellCommInd']])) {
		stop('There is a mismatch between the raster and the species cell list. Different number of elements!')
	}
	grid <- raster::xyFromCell(x[[1]], 1:raster::ncell(x[[1]]))
	cell <- identify(x = grid[,1], y = grid[,2], plot = FALSE, n = 1, ...)
	if (returnCell) {
		return(cell)
	} else {
		return(x[['speciesList']][[x[['cellCommInd']][cell]]])
	}
}