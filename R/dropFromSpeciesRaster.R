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
	
	if (!inherits(x, 'speciesRaster')) {
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
		
		return(rebuildSpeciesRaster(x, fullList))
		
	} else {
		return(x)
	}
}