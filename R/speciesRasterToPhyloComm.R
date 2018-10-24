##' @title Convert speciesRaster to community matrix
##'
##' @description Given specific sites, convert speciesRaster to 
##' 	phylocomm matrix, with sites as rows, and species as columns
##'
##' @param x object of class \code{speciesRaster}
##' @param sites locations of sites, see details. 
##'
##' @details If sites are site coordinates, 
##' 	then dataframe or matrix with two columns; 
##' 	if sites are cell indices, then numeric vector;
##' 	if \code{sites = 'all'}, then all cells will be returned as sites. 
##'
##' @return community matrix, with sites as rows and species as columns
##'
##' @author Pascal Title
##'
##' @examples
##' library(raster)
##' tamiasSpRas
##'
##' # from cell indices
##' cells <- c(22440, 41446, 39283, 37697)
##' speciesRasterToPhyloComm(tamiasSpRas, cells)
##' 
##' # from coordinates
##' xy <- xyFromCell(tamiasSpRas[[1]], cells)
##' speciesRasterToPhyloComm(tamiasSpRas, xy)
##' 
##' @export

speciesRasterToPhyloComm <- function(x, sites) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	# depending on input, make sites a vector of cell indices
	if (class(sites) %in% c('matrix', 'data.frame')) {
		sites <- raster::cellFromXY(x[[1]], sites)
	} else if (class(sites) == 'character') {
		sites <- 1:raster::ncell(x[[1]])
	}
	
	# extract relevant cells from speciesRaster
	dat <- sapply(x[['cellCommInd']][sites], function(y) x[['speciesList']][[y]])
	uniqueSp <- sort(unique(unlist(dat)))
	uniqueSp <- uniqueSp[stats::complete.cases(uniqueSp)]
	
	# build phylocomm matrix
	resMat <- matrix(nrow = length(dat), ncol = length(uniqueSp))
	rownames(resMat) <- paste0('site', 1:length(dat))
	colnames(resMat) <- uniqueSp
	
	for (i in 1:nrow(resMat)) {
		resMat[i,] <- as.numeric(uniqueSp %in% dat[[i]])
	}
	
	return(resMat)
}