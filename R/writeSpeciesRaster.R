##' @title Save speciesRaster object
##'
##' @description Write a speciesRaster object to disk. 
##'
##' @param x object of class \code{speciesRaster}
##' @param filename filename with no extension
##'
##' @details This function ensures that the raster data is stored in memory
##' and then writes a .rds file with xz compression. This file can be read back in with \code{\link{readRDS}}.
##'
##' @return Nothing is returned, but object is written to disk.
##'
##' @author Pascal Title
##'
##' @examples
##' \dontrun{
##' #save
##' save(tamiasSpRas, '~/tamiasSpeciesRaster')
##' 
##' # read back in
##' tamiasSpRas <- readRDS('~/tamiasSpeciesRaster.rds')
##' }
##' @export

writeSpeciesRaster <- function(x, filename) {
	
	if (!inherits(x, 'speciesRaster')) {
		stop('x must be of class speciesRaster.')
	}

	# check if raster values are in memory, and if not, move them to memory
	if (!raster::inMemory(x[[1]])) {
		x[[1]] <- raster::setValues(x[[1]], 1:raster::ncell(x[[1]]))
	}
	
	if (!grepl('\\.rds', filename)) {
		filename <- paste0(filename, '.rds')
	}
	
	saveRDS(x, file = filename, compress = 'xz')	
}

