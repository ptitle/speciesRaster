##' @title Map turnover in beta diversity
##'
##' @description Taxonomic and phylogenetic community dissimilarity is calculated for each cell 
##' 	within a moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##' @param dist resolution of window size, where the multisite metric will be calculated 
##' 	for all smaller cells found within these coarser cells
##' 
##' @param metric choice of metric, see details.
##'
##' @param nthreads number of threads to use for parallelization of the function. 
##' 	The R package \code{parallel} must be loaded for \code{nthreads > 1}.
##'
##' @param verbose Intended primarily for debugging, prints progress to the console
##'
##'
##' @details This function calculates taxonomic and phylogenetic turnover in species communities,
##' 	as defined in the \code{speciesRaster} object. A more coarse-grained raster, with cell size defined
##' 	by the parameter \code{dist}, is overlaid on top of the \code{speciesRaster}, and an indexing is done
##' 	to determine which cells fall within the larger coarse-grained cells. For every coarse-grained 
##'		cell, a multi-site metric is calculated from the sites that fall within that cell.
##'		
##'		Following Baselga (2010) and Leprieur et al. (2012), the overall beta diversity metric can be 
##' 	decomposed into two additive metrics: the turnover component and the nestedness component of
##'		beta diversity. Equivalent metrics are available from two families: Sorensen and Jaccard.
##'		See \code{\link{beta.multi}} and \code{\link{phylo.beta.multi}} from the betapart package for
##'		more information.
##'
##' 	Metric options include:
##'
##'		Sorensen:
##' 	\itemize{
##' 		\item{betaSOR:} {Sorensen dissimilarity}
##' 		\item{betaSIM:} {turnover component of Sorensen dissimilarity}
##' 		\item{betaSNE:} {nestedness-resultant component of Sorensen dissimilarity}
##' 		\item{phyloBetaSOR:} {phylogenetic Sorensen dissimilarity}
##' 		\item{phyloBetaSIM:} {phylogenetic turnover component of Sorensen dissimilarity}
##' 		\item{phyloBetaSNE:} {phylogenetic nestedness-resultant component of Sorensen dissimilarity}
##'		}
##'		Jaccard:
##'		\itemize{
##' 		\item{betaJAC:} {Jaccard dissimilarity}
##' 		\item{betaJTU:} {turnover component of Jaccard dissimilarity}
##' 		\item{betaJNE:} {nestedness-resultant component of Jaccard dissimilarity}
##' 		\item{phyloBetaJAC:} {phylogenetic Jaccard dissimilarity}
##' 		\item{phyloBetaJTU:} {phylogenetic turnover component of Jaccard dissimilarity}
##' 		\item{phyloBetaJNE:} {phylogenetic nestedness-resultant component of Jaccard dissimilarity}
##'		}
##'
##' 		
##' 
##' @return Returns a RasterObject, where the metric can range from 0 to 1.
##' 
##' @author Pascal Title
##'
##' @references
##' Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global 
##' Ecology and Biogeography 19:134-143
##' 
##' Leprieur F, Albouy C, De Bortoli J, Cowman PF, Bellwood DR, et al. (2012) Quantifying Phylogenetic 
##' Beta Diversity: Distinguishing between "True" Turnover of Lineages and Phylogenetic Diversity Gradients. 
##' PLoS ONE 7(8): e42760. doi:10.1371/journal.pone.0042760
##' 
##' @examples
##' library(raster)
##' tamiasSpRas
##' tamiasTree
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' 
##' # using a moving window of roughly 50km radius
##' betaSOR <- betaDiversityBlock(tamiasSpRas, dist = 50000, metric = 'betaSOR')
##' 
##' phyloBetaSOR <- betaDiversityBlock(tamiasSpRas, dist = 50000, metric = 'phyloBetaSOR')
##' 
##' colramp <- colorRampPalette(c('blue','yellow','red'))
##' par(mfrow=c(1,2))
##' plot(betaSOR, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
##' plot(phyloBetaSOR, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
##' addRasterLegend(phyloBetaSOR, location = 'left', minmax = c(0, 1), ramp = c(c('blue','yellow','red')))


##' @export




betaDiversityBlock <- function(x, dist = 50000, metric = 'betaSOR', nthreads = 1, verbose = FALSE) {
	# radius is distance in map units = meters, will round up to the nearest cell
		
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}	

	# check metric validity
	if (length(metric) > 1) {
		stop('Only one metric can be specified.')
	}
	metric <- match.arg(metric, choices = c('betaSOR', 'betaSIM', 'betaSNE', 'betaJAC', 'betaJTU', 'betaJNE', 'phyloBetaSOR', 'phyloBetaSIM', 'phyloBetaSNE', 'phyloBetaJAC', 'phyloBetaJTU', 'phyloBetaJNE'))
	if (!metric %in% c('betaSOR', 'betaSIM', 'betaSNE', 'betaJAC', 'betaJTU', 'betaJNE', 'phyloBetaSOR', 'phyloBetaSIM', 'phyloBetaSNE', 'phyloBetaJAC', 'phyloBetaJTU', 'phyloBetaJNE')) {
		stop('Invalid metric.')
	}
	
	phyloMetric <- ifelse(grepl('phylo', metric), TRUE, FALSE)
	family <- ifelse(grepl('betaS', metric, ignore.case = TRUE), 'sorensen', 'jaccard')
	
	if (phyloMetric & is.null(x[['phylo']])) {
		stop('speciesRaster object does not contain a phylo object!')
	}
	
	# define raster template at dist resolution
	if (verbose) cat('\t...determining which cells belong to each window...\n')
	template <- raster::raster(ext = raster::extent(x[[1]]), res = c(dist, dist), crs = raster::projection(x[[1]]))
	
	cellPoly <- raster::rasterToPolygons(template)
	cellPoly <- lapply(1:length(cellPoly), function(y) raster::extent(cellPoly[y,]))
	ras <- x[[1]]
	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads)
		parallel::clusterExport(cl = cl, varlist = c('cellPoly', 'ras', 'cellsFromExtent'), envir = environment())
		eList <- parallel::parLapply(cl = cl, cellPoly, function(y) raster::cellsFromExtent(ras, y))
		parallel::stopCluster(cl)
	} else {
		eList <- lapply(cellPoly, function(y) raster::cellsFromExtent(ras, y))
	}
	
	nbLengths <- vapply(eList, length, FUN.VALUE = 1L)
	if (verbose) cat('\t...range of cells per window:', range(nbLengths)[1], '-', range(nbLengths)[2], '\n')
					
	
	if (phyloMetric) {
		# phylogenetic beta diversity metric
		metricTranslate <- setNames(c('phyloBetaSOR', 'phyloBetaSIM', 'phyloBetaSNE', 'phyloBetaJAC', 'phyloBetaJTU', 'phyloBetaJNE'), c('phylo.beta.SOR', 'phylo.beta.SIM', 'phylo.beta.SNE', 'phylo.beta.JAC', 'phylo.beta.JTU', 'phylo.beta.JNE'))
		metric <- names(metricTranslate)[which(metricTranslate == metric)]
		tree <- x[['phylo']]
		if (verbose) cat('\t...Phylogenetic metric selected:', metric, '\n')
		
		if (verbose) cat('\t...dropping species that are not in phylo data...\n')
	 	# prune speciesRaster object down to species shared with phylogeny
		x[[2]] <- intersectList(x[[2]], x[['phylo']]$tip.label)
		
		if (nthreads > 1) {
			cl <- parallel::makePSOCKcluster(nthreads)
			parallel::clusterExport(cl = cl, varlist = c('x', 'eList', 'speciesRasterToPhyloComm'), envir = environment())
			
			commList <- parallel::parLapply(cl, eList, function(y) speciesRasterToPhyloComm(x, y))
			
			parallel::clusterEvalQ(cl = cl, library(betapart))
			parallel::clusterExport(cl = cl, varlist = c('commList', 'tree', 'family', 'metric'), envir = environment())
			cellBeta <- rep(NA, length(commList))
			cellBeta[sapply(commList, ncol) == 1] <- 0
			cellBeta[sapply(commList, ncol) > 1] <- parallel::parSapply(cl = cl, commList[sapply(commList, ncol) > 1], function(y) {
				return(betapart::phylo.beta.multi(y, tree, index.family = family)[[metric]])
			})
			parallel::stopCluster(cl)
		} else {
			commList <- lapply(eList, function(y) speciesRasterToPhyloComm(x, y))
			cellBeta <- rep(NA, length(commList))
			cellBeta[sapply(commList, ncol) == 1] <- 0
			cellBeta[sapply(commList, ncol) > 1] <- sapply(commList[sapply(commList, ncol) > 1], function(y) {
				betapart::phylo.beta.multi(y, tree, index.family = family)[[metric]]
			})
		}
		
	} else {
		# taxonomic beta diversity metric
		if (verbose) cat('\t...Taxonomic metric selected:', metric, '\n')
		cellBeta <- calcBetaMultiSiteBlock(x[[2]], eList, metric)
	}
	
	res <- template
	raster::values(res) <- cellBeta
	return(res)	
}