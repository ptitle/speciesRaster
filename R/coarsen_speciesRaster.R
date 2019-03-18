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
##' @details Species will be included in the coarser cell if they are found in
##' 	>= 50 percent of the original-resolution cells that are being aggregated. 
##' 	This should be thought of as more of a convenience function for experimentation.
##' 	Creating a coarser resolution speciesRaster from species ranges with 
##'		\code{\link{createSpeciesRaster}} is the more appropriate approach.
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
	
	if (!inherits(x, 'speciesRaster')) {
		stop('x must be of class speciesRaster.')
	}
	
	template <- raster::aggregate(x[[1]], fact)
	raster::values(template) <- 1:raster::ncell(template)
	
	# determine which large cell each small cell falls in
	cellPts <- raster::xyFromCell(x[[1]], 1:raster::ncell(x[[1]]))
	e <- raster::extract(template, cellPts)
	
	eList <- lapply(1:raster::ncell(template), function(y) which(e == y))
	
	# unroll condensed sp x cell list
	spCellList <- expandSpeciesCellList(x)
	
	# for each new cell, find the unioned set of species found in all included cells
	# for each set of cells, keep those species that are found in 50% or more of the cells
	newSpList <- vector('list', length = raster::ncell(template))
	for (i in 1:length(eList)) {
		newSpList[[i]] <- keepMajoritySpecies(spCellList[eList[[i]]])		
	}
	
	newSpRas <- rebuildSpeciesRaster(x, newSpList, template)
	
	return(newSpRas)

}

# function returns the species found in 50% or more of the list entries
keepMajoritySpecies <- function(z) {
	
	cutoff <- length(z) / 2
	uniqueSp <- Reduce(union, z)
	uniqueSp <- uniqueSp[stats::complete.cases(uniqueSp)]
	
	if (length(uniqueSp) > 0) {
		majoritySp <- sapply(uniqueSp, function(y) sum(sapply(z, function(yy) y %in% yy)))
		majoritySp <- names(majoritySp)[majoritySp >= cutoff]
		if (length(majoritySp) > 0) {
			return(majoritySp)
		} else {
			return(NA)
		}
		
	} else {
		return(NA)
	}
}