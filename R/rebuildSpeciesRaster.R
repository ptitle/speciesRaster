##' @title Rebuild speciesRaster
##'
##' @description Given a list of species in each cell, all components of the
##' 	speciesRaster object are reindexed and regenerated.
##'
##' @param x object of class \code{speciesRaster}
##' @param spCellList list in which each element is a character vector of species names.
##' @param rasterTemplate a separate raster can be provided if it has different 
##' 	characteristics than the one in \code{x}.
##'
##' @details This function is used internally by \code{\link{dropFromSpeciesRaster}} and 
##' 	by \code{\link{coarsen_speciesRaster}}. 
##'
##' @return new \code{speciesRaster} object.
##'
##' @author Pascal Title
##'	
##' 
##' @export


rebuildSpeciesRaster <- function(x, spCellList, rasterTemplate = NULL) {
	
	if (!is.null(rasterTemplate)) {
		if (raster::ncell(rasterTemplate) != length(spCellList)) {
			stop('ncells of rasterTemplate not equal to length of spCellList.')
		}
		
		x[[1]] <- rasterTemplate
	}
	
	# convert spCellList into condensed version: 
	# 	list of unique cell communities
	#	vector of indices to map these back to raster cells
	cellCommVec <- integer(length = raster::ncell(x[[1]]))
	uniqueComm <- unique(spCellList)
	fullList2 <- sapply(spCellList, function(y) paste(y, collapse = '|'))
	uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
	for (i in 1:length(uniqueComm2)) {
		cellCommVec[which(fullList2 == uniqueComm2[i])] <- i
	}
	
	uniqueSp <- sort(unique(unlist(uniqueComm)))
	
	spCellCount <- countCells(convertNAtoEmpty(spCellList), uniqueSp)
	names(spCellCount) <- uniqueSp
	
	rasterVals <- lengths(spCellList)
	rasterVals[which(sapply(spCellList, anyNA))] <- NA
	raster::values(x[[1]]) <- rasterVals
	
	x[['speciesList']] <- uniqueComm
	x[['cellCommInd']] <- cellCommVec
	x[['geogSpecies']] <- uniqueSp
	x[['cellCount']] <- spCellCount

	return(x)
}