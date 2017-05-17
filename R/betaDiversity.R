##' @title Map turnover in beta diversity
##'
##' @description Taxonomic and phylogenetic community dissimilarity is calculated for each cell 
##' 	within a moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##'
##' @param cellNum resolution of neighborhood size (horizontal and vertical), where the multisite metric 
##' 	will be calculated for all smaller cells found within these coarser cells
##' 
##' @param metric choice of metric, see details.
##'
##' @param minNB the minimum number of communities to consider for calculation of the beta metric 
##' 	metric. See Details. 
##'
##' @param nthreads number of threads to use for parallelization of the function. 
##' 	The R package \code{parallel} must be loaded for \code{nthreads > 1}.
##'
##' @param verbose Intended primarily for debugging, prints progress to the console
##'
##'
##' @details This function calculates taxonomic and phylogenetic turnover in species communities,
##' 	as defined in the \code{speciesRaster} object. A more coarse-grained raster, with cell size defined
##' 	by the parameter \code{cellNum}, is overlaid on top of the \code{speciesRaster}, and an indexing is done
##' 	to determine which cells fall within the larger coarse-grained cells. For every coarse-grained 
##'		cell, a multi-site metric is calculated from the sites that fall within that cell. As differences
##' 	in the number of communities could bias the beta diversity metrics, only coarse-grained cells with 
##' 	the same number of communities are retained.
##'		
##'		Following Baselga (2010) and Leprieur et al. (2012), the overall beta diversity metric can be 
##' 	decomposed into two additive metrics: the turnover component and the nestedness component of
##'		beta diversity. Equivalent metrics are available from two families: Sorensen and Jaccard.
##'		See \code{\link{beta.multi}} and \code{\link{phylo.beta.multi}} from the betapart package for
##'		more information.
##'
##'		\code{minNB:} Depending on the value of \code{cellNum} and the spatial distribution of non-empty
##'		cells in the speciesRaster, the number of communities contributing to the calculation of the beta 
##'		diversity metric might not always be the same. At the edges of the raster, or along coastlines, for
##'		example, fewer cells may contribute to the metric calculation, which might lead to bias. If the user
##'		would like to restrict calculations to a minimum number of communities, this can be specified with 
##'		\code{minNB} (minimum neighborhood size). By default, this is set to half of the maximum number of
##'		communities that are possible. To include all communities, set to 1.
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
##' # using a neighborhood size of 4 cells
##' betaSOR <- betaDiversity(tamiasSpRas, cellNum = 4, metric = 'betaSOR')
##' 
##' phyloBetaSOR <- betaDiversity(tamiasSpRas, cellNum = 4, metric = 'phyloBetaSOR')
##' 
##' colramp <- colorRampPalette(c('blue','yellow','red'))
##' par(mfrow=c(1,2))
##' plot(betaSOR, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
##' plot(phyloBetaSOR, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
##' addRasterLegend(phyloBetaSOR, location = 'left', minmax = c(0, 1), ramp = c(c('blue','yellow','red')))


##' @export




betaDiversity <- function(x, cellNum = 4, metric = 'betaSOR', minNB = (cellNum ^ 2) / 2, nthreads = 1, verbose = FALSE) {
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
	
	# define raster template at cellNum resolution
	if (verbose) cat('\t...resolution of output grid:', cellNum * raster::res(x[[1]])[1], 'by', cellNum * raster::res(x[[1]])[1], '...\n')
	if (verbose) cat('\t...determining which cells belong to each neighborhood...')
	template <- raster::aggregate(x[[1]], cellNum)
	raster::values(template) <- 1:raster::ncell(template)
	
	# determine which large cell each small cell falls in
	cellPts <- raster::xyFromCell(x[[1]], 1:raster::ncell(x[[1]]))
	e <- raster::extract(template, cellPts)
	
	eList <- lapply(1:raster::ncell(template), function(y) which(e == y))
	
	if (verbose) cat('done.\n')
	
	if (minNB > 1) {
		datCells <- which(!is.na(raster::values(x[[1]])))
		eList2 <- lapply(eList, function(y) intersect(y, datCells))
		nbLengths <- vapply(eList2, length, FUN.VALUE = 1L)
	}
		
	# how is neighborhood size spatially distributed?
	# mapping <- template
	# raster::values(mapping) <- nbLengths
	
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
	
	if (minNB > 1) {
		# drop cells that do not have the full number of communities
		cellBeta[which(nbLengths < minNB)] <- NA
	}
	
	res <- template
	raster::values(res) <- cellBeta
	return(res)	
}


