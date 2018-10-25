##' Function to expand condensed species list to full set of cells

##' @title Expand species list
##'
##' @description The speciesRaster object contains an accounting of species per
##' cell in a condensed format. This function returns a complete list of species
##' per cell. 
##'
##' @param x object of class \code{speciesRaster}
##'
##' @return list of species for each cell.
##'
##' @author Pascal Title
##' 
##' @examples
##' tamiasSpRas
##'	head(expandSpeciesCellList(tamiasSpRas))
##' 
##' @export



expandSpeciesCellList <- function(x) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}	
	
	return(sapply(x[['cellCommInd']], function(y) x[['speciesList']][[y]]))
	
}