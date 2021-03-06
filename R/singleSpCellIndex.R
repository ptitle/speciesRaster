##' @title Identify single-species cells
##'
##' @description Given a speciesRaster object, return the raster cell indices
##' 	of those cells that have just one species. 
##'
##' @param x object of class \code{speciesRaster}
##'
##' @details This function can be useful when further analyzing speciesRaster
##' 	objects generated by \code{\link{cellMetrics_speciesRaster}}, as it might
##' 	make sense to exclude these single-species cells in further analyses. 
##'
##' 
##' @return numeric vector of raster cell numbers. 
##'
##' @author Pascal Title
##' 
##' @examples
##' singleSpCellIndex(tamiasSpRas)
##' 
##' @export


# identify single-species cells

singleSpCellIndex <- function(x) {

	if (!inherits(x, 'speciesRaster')) {
		stop('Object must be of class speciesRaster')
	}
	
	singleSp <- which(lengths(x[['speciesList']]) == 1 & sapply(x[['speciesList']], anyNA) == FALSE)
	singleSpCells <- which(x[['cellCommInd']] %in% singleSp)
	return(singleSpCells)
}