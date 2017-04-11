##' @title Map turnover in beta diversity
##'
##' @description Mean community dissimilarity is calculated for each cell within a moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##' @param radius Radius length in map units to define the moving window. 
##' 	This will be converted by rounding up to a number of neighboring 
##' 	cells in all directions. 
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
##' beta <- betaDiversity(tamiasSpRas, radius = 50000, metric = 'betaSOR')
##' plot(beta)
##' 
##' @export


# create mapping of community turnover with the jaccard index.
# first, apply a moving window, where for each cell, the mean of the jaccard index between that cell and each of its neighbors within the window is calculated. 
# maybe other option: calculate all cell pairwise jaccard indices, to then run nmds on it.

betaDiversity <- function(x, radius = 50000, metric = 'betaSOR') {
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
	
	#Identify cells for moving window
	#create neighbor matrix from window size
	window <- ceiling(radius / raster::res(x[[1]])[[1]])
	side <- window * 2 + 1
	counter <- window
	nb <- vector()
	for (i in 1:side) {
		nb <- rbind(nb, c(rep(NA,counter), rep(1,(side - counter*2)), rep(NA, counter)))
		if (i < side/2) {counter <- counter - 1}
		if (i > side/2) {counter <- counter + 1}
	}
	nb[window + 1, window + 1] <- 0

	#return neighbor cells for each cell
	# the order of the cells matches the order of spGridList
	all_nb <- raster::adjacent(x[[1]], cells = 1:raster::ncell(x[[1]]), directions = nb, include=TRUE)
	nbList <- split(all_nb[,2], all_nb[,1])
			
	cellBeta <- calcBetaMultiSite(x[[2]], nbList, metric)
	
	res <- x
	raster::values(res[[1]]) <- cellBeta
	names(res[[1]]) <- metric
	return(res)	
}
