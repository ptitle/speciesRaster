##' @title Get data from speciesRaster
##'
##' @description Return the data item from a speciesRaster object. 
##'
##' @param x object of class \code{speciesRaster}
##'
##' @return The object found in the \code{data} location of the speciesRaster.
##' 	Could be either a vector or a matrix. 
##'
##' @author Pascal Title
##'
##' @examples
##' x <- addTraits_speciesRaster(tamiasSpRas, tamiasTraits)
##' getData_speciesRaster(x)
##' 
##' @export

getData_speciesRaster <- function(x) {
	
	if (!inherits(x, 'speciesRaster')) {
		stop('x must be of class speciesRaster.')
	}
		
	if (!inherits(x[['data']], c('numeric', 'matrix', 'data.frame'))) {
		stop('Data not present.')
	}
	
	return(x[['data']])	
}