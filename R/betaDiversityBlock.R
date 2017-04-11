##' @title Map turnover in beta diversity
##'
##' @description Mean community dissimilarity is calculated for each cell within a moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##' @param dist resolution of window size, where the multisite metric will be calculated for all smaller cells found within these coarser cells
##' @param metric choice of metric, see details.
##'
##' @details Should put more information about the specific metric here, 
##' eventually implement a few complimentary ones.
##' Currently: Multi-site metrics implemented for: Sorensen dissimilarity (betaSOR), Simpson
##' dissimlarity / turnover component of Sorensen dissimilarity (betaSIM), nestedness-resultant
##' component of Sorensen dissimilarity (betaSNE), jaccard dissimilarity (betaJAC), turnover component
##' of jaccard dissimilarity (betaJTU), nestedness-resultant component of Jaccard dissimilarity
##' (betaJNE).
##' 
##' @return Returns a new \code{speciesRaster} object, with mean community dissimilarity for each cell.
##' 
##' @author Pascal Title
##' 
##' @examples
##' tamiasSpRas
##' 
##' # using a moving window of roughly 50km radius
##' beta <- betaDiversityBlock(tamiasSpRas, dist = 50000, metric = 'betaSOR')
##' plot(beta)
##' 
##' @export




betaDiversityBlock <- function(x, dist = 50000, metric = 'betaSOR') {
	# radius is distance in map units = meters, will round up to the nearest cell
		
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}

	# check metric validity
	if (length(metric) > 1) {
		stop('Only one metric can be specified.')
	}
	metric <- match.arg(metric, choices = c('betaSOR', 'betaSIM', 'betaSNE', 'betaJAC', 'betaJTU', 'betaJNE'))
	if (!metric %in% c('betaSOR', 'betaSIM', 'betaSNE', 'betaJAC', 'betaJTU', 'betaJNE')) {
		stop('Invalid metric.')
	}

	# define raster template at dist resolution
	template <- raster::raster(ext = raster::extent(x[[1]]), res = c(dist, dist), crs = raster::projection(x[[1]]))

	eList <- lapply(1:raster::ncell(template), function(y) raster::extent(raster::rasterFromCells(template, y, values = FALSE)))
	eList <- lapply(eList, function(y) raster::cellsFromExtent(x[[1]], y))
				
	# calculate cell values
	cellBeta <- calcBetaMultiSiteBlock(x[[2]], eList, metric)
	
	# get new species by cell

	res <- x
	res[[1]] <- template
	raster::values(res[[1]]) <- cellBeta[[1]]
	res[[2]] <- cellBeta[[2]]
	names(res[[1]]) <- metric
	return(res)	
}
