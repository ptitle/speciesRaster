##' @title Get phylogeny from speciesRaster
##'
##' @description Return the phylo item from a speciesRaster object. 
##'
##' @param x object of class \code{speciesRaster}
##'
##' @return The phylogeny in the \code{phylo} location of the speciesRaster.
##'
##' @author Pascal Title
##' 
##' @examples
##' x <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' getPhylo_speciesRaster(x)
##' 
##' @export

getPhylo_speciesRaster <- function(x) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
		
	if (!class(x[['phylo']]) %in% 'phylo') {
		stop('Phylogeny not present.')
	}
	
	return(x[['phylo']])	
}