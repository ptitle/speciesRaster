##' @title Drop species from speciesRaster
##'
##' @description Removes particular species from a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param sp a character vector of species names to be dropped.
##'
##' @details If species in \code{sp} are not in \code{x}, they will be ignored.
##'
##' @return new \code{speciesRaster} object.
##'
##' @author Pascal Title
##' 
##' @examples
##' tamiasSpRas
##'
##' newSpRas <- dropFromSpeciesRaster(tamiasSpRas, sp = c('Tamias_alpinus', 'Tamias_bulleri'))
##'
##' setdiff(tamiasSpRas[['geogSpecies']], newSpRas[['geogSpecies']])
##'	
##' 
##' @export



dropFromSpeciesRaster <- function(x, sp) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	if (any(sp %in% x[['geogSpecies']])) {
	
		fullList <- expandSpeciesCellList(x)
		
		for (i in 1:length(fullList)) {		
			if (any(sp %in% fullList[[i]])) {
				tmp <- sort(setdiff(fullList[[i]], sp))
				if (length(tmp) > 0) {
					fullList[[i]] <- tmp
				} else {
					fullList[[i]] <- NA
				}
			}		
		}
		
		raster::values(x[[1]]) <- lengths(fullList)
		raster::values(x[[1]])[which(sapply(fullList, anyNA) == TRUE)] <- NA
		
		# reduce spByCell to unique communities and track
		cellCommVec <- integer(length = raster::ncell(x[[1]]))
		uniqueComm <- unique(fullList)
		fullList2 <- sapply(fullList, function(y) paste(y, collapse = '|'))
		uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
		for (i in 1:length(uniqueComm2)) {
			cellCommVec[which(fullList2 == uniqueComm2[i])] <- i
		}
		
		x[['geogSpecies']] <- setdiff(x[['geogSpecies']], sp)
		x[['cellCount']] <- x[['cellCount']][x[['geogSpecies']]]
	
		x[['speciesList']] <- uniqueComm
		x[['cellCommInd']] <- cellCommVec
	}

	return(x)
}