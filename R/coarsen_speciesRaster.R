##' @title Coarsen a speciesRaster
##'
##' @description Return a speciesRaster object with a coarser resolution.
##'
##' @param x object of class \code{speciesRaster}

##' @param fact the aggregation factor as the number of cells in the horizontal
##'		and vertical directions.
##'
##' @return \code{speciesRaster} with a coarser resolution. 
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasSpRas
##' 
##' # coarsen such that groups of 4 cells are aggregated to 1 cell
##' # (2 cells in horizontal and 2 cells in vertical directions)
##' coarsen_speciesRaster(tamiasSpRas, fact = 2)
##' 
##' @export


coarsen_speciesRaster <- function(x, fact) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	template <- raster::aggregate(x[[1]], fact)
	raster::values(template) <- 1:raster::ncell(template)
	
	# determine which large cell each small cell falls in
	cellPts <- raster::xyFromCell(x[[1]], 1:raster::ncell(x[[1]]))
	e <- raster::extract(template, cellPts)
	
	eList <- lapply(1:raster::ncell(template), function(y) which(e == y))
	
	# for each new cell, find the unioned set of species found in all included cells
	newSpList <- vector('list', length = raster::ncell(template))
	for (i in 1:length(eList)) {
		tmp <- Reduce(union, x[['speciesList']][eList[[i]]])
		tmp <- tmp[stats::complete.cases(tmp)]
		newSpList[[i]] <- tmp
	}
	
	newSpList[which(sapply(newSpList, length) == 0)] <- NA
		
	res <- x
	res[[1]] <- template
	raster::values(res[[1]]) <- sapply(newSpList, length)
	raster::values(res[[1]])[sapply(newSpList, anyNA)] <- NA
	res[[2]] <- newSpList
	
	return(res)

}