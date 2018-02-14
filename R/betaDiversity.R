##' @title Map turnover in beta diversity
##'
##' @description Mean community dissimilarity is calculated for each cell within a moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##' @param radius Radius in terms of number of cells in horizontal and vertical directions from 
##' 	the focal cell to define the moving window. 
##' @param metric choice of metric, see details.
##' @param verbose Primarily intended for debugging, print progress to the console.
##'
##' @details
##' Sorensen's Beta diversity \code{sorensen} is a purely taxonomic measure of turnover.
##' 
##' Phylogenetic Sorensen's Similarity \code{phylosor} is a measure of phylogenetic turnover, 
##' 	ranging from 0 to 1. This function returns \code{1 - phylosor}, such that little 
##' 	change in phylogenetic community structure results in values closer to 0. 
##'
##' Range-weighted turnover \code{RWTurnover} measures turnover but where taxa are weighted
##' 	according to the inverse of their range size.
##' 
##' Phylogenetic range-weighted turnover \code{phyloRWTurnover} measures turnover in phylogenetic diversity
##' 	where phylogenetic branches are weighted by the inverse of their geographic distribution.
##'
##' 
##' The radius value defines the number of cells counted out from the focal cell in the horizontal and vertical
##' 	direction, creating a square window of \code{side = radius * 2 + 1}. For each cell, all calculations
##' 	between the focal cell and a moving window cell are run, and the results are averaged and assigned 
##' 	to that focal cell. 
##' 
##' @return Returns a new \code{speciesRaster} object, with mean community dissimilarity for each cell.
##' 
##' @author Pascal Title
##'
##' @references
##' Laffan, SW, et al. Range-weighted metrics of species and phylogenetic turnover can better 
##' resolve biogeographic transition zones. Methods in Ecology and Evolution 7 (2016): 580-588.
##'
##' Rosauer, D, Laffan, SW, Crisp, MD, Donnellan, SC, Cook, LG. Phylogenetic endemism: a new approach 
##' for identifying geographical concentrations of evolutionary history. Molecular Ecology 
##' 18 (2009): 4061-4072.
##' 
##' @examples
##' library(raster)
##' tamiasSpRas
##' 
##' beta_taxonomic <- betaDiversity_speciesRaster(tamiasSpRas, radius = 4, metric = 'sorensen')
##' 
##' beta_RW <- betaDiversity_speciesRaster(tamiasSpRas, radius = 4, metric = 'RWTurnover')
##' 
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' beta_phyloRW <- betaDiversity_speciesRaster(tamiasSpRas, radius = 4, metric = 'phyloRWTurnover')
##' 
##' colramp <- colorRampPalette(c('blue','yellow','red'))
##' par(mfrow=c(1,3))
##' plot(beta_taxonomic, col = colramp(100))
##' plot(beta_RW, col = colramp(100))
##' plot(beta_phyloRW, col = colramp(100))
##' 
##' 
##' @export








betaDiversity_speciesRaster <- function(x, radius = 3, metric, verbose = FALSE) {
	# radius is distance in cells
		
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}

	# check metric validity
	if (length(metric) > 1) {
		stop('Only one metric can be specified.')
	}
	
	if (radius < 1 | radius != trunc(radius)) {
		stop('radius must be an integer >= 1.')
	}
	
	metric <- match.arg(metric, choices = c('sorensen', 'RWTurnover', 'phyloRWTurnover', 'phylosor'))
	if (!metric %in% c('sorensen', 'RWTurnover', 'phyloRWTurnover', 'phylosor')) {
		stop('Invalid metric.')
	}
	
	# if phylogenetic, there must be a phylo object
	# check that there is a phylogeny in speciesRaster object
	if (metric %in% c('phylosor', 'phyloRWTurnover')) {
		if (is.null(x[['phylo']])) {
			stop('speciesRaster object does not contain a phylo object!')
		}
		
		if (verbose) cat('\t...dropping species that are not in phylo data...\n')
 		# prune speciesRaster object down to species shared with phylogeny
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']]$tip.label)
		raster::values(x[[1]]) <- lengths(x[['speciesList']])
		raster::values(x[[1]])[which(sapply(x[['speciesList']], anyNA) == TRUE)] <- NA
		
		if (metric == 'phyloRWTurnover') {
			spEdges <- getRootToTipEdges(x[['phylo']])
		
			if (!'edgeArea' %in% names(x)) {
				if (verbose) cat('\t...calculating branch-specific range sizes...\n')
				x[['edgeArea']] <- do.call(cbind, phyloBranchRanges(x[['phylo']], convertNAtoEmpty(x[['speciesList']]), spEdges))
			}
		}		
	} 
		
	# identify non-empty cells
	nonNAcells <- which(!is.na(raster::values(x[[1]])))

	# return neighbor cells for each cell	
	if (verbose) cat('\t...identifying valid cell neighbors...\n')
	all_nb <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = 8)
	
	# remove any cells that have no neighbors (isolated cells)
	nonNAcells <- intersect(nonNAcells, all_nb[,1])

	# subset species cell list to non-empty cells
	spCellList <- x[['speciesList']][nonNAcells]
	
	if (anyNA(unlist(spCellList))) {
		stop('Empty communities in list. Please inform the developer.')
		
		# zz <- which(sapply(x[[2]], anyNA) == TRUE)
		# test <- raster::values(x[[1]])[zz]
		# unique(test)
	}

	# remap the cell numbers to the list with empty cells removed
	if (verbose) cat('\t...remapping cell indexing...\n')
	cellMap <- x[[1]]
	cellMap <- rep(-1, raster::ncell(x[[1]]))
	cellMap[nonNAcells] <- 0:(length(nonNAcells) - 1)
	
	# display progress bar if large dataset
	#showProgress <- length(spCellList) > 100000
	showProgress <- TRUE
	
	if (verbose) cat(paste0('\t...calculating metric ', metric, '\n'))
	if (metric == 'sorensen') {
		cellBeta <- calcRWTurnover_taxonomic(spCellList, radius, rasterNRow = nrow(x[[1]]), rasterNCol = ncol(x[[1]]), cellMap, nonNAcells, showProgress = showProgress)
	} else if (metric == 'RWTurnover') {
		cellBeta <- calcRWTurnover_rangeWeighted(spCellList, radius, rasterNRow = nrow(x[[1]]), rasterNCol = ncol(x[[1]]), cellMap, nonNAcells, 1 / x[['cellCount']], showProgress = showProgress)
	} else if (metric == 'phyloRWTurnover') {
		cellBeta <- calcRWTurnover_phyloRangeWeighted(spCellList, radius, rasterNRow = nrow(x[[1]]), rasterNCol = ncol(x[[1]]), cellMap, nonNAcells, x[['phylo']], spEdges, x[['edgeArea']], showProgress = showProgress)
	} else if (metric == 'phylosor') {
		cellBeta <- 1 - calcPhylosor(spCellList, radius, rasterNRow = nrow(x[[1]]), rasterNCol = ncol(x[[1]]), cellMap, nonNAcells, x[['phylo']], showProgress = showProgress)
	}
	
	cellVal <- rep(NA, raster::ncell(x[[1]]))
	cellVal[nonNAcells] <- cellBeta
	res <- x[[1]]
	raster::values(res) <- cellVal
	names(res) <- metric
	return(res)	
}



# betaDiversity_speciesRaster <- function(x, radius = 3, metric, verbose = FALSE) {
	# # radius is distance in cells
		
	# if (!'speciesRaster' %in% class(x)) {
		# stop('x must be of class speciesRaster.')
	# }

	# # check metric validity
	# if (length(metric) > 1) {
		# stop('Only one metric can be specified.')
	# }
	
	# if (radius < 1 | radius != trunc(radius)) {
		# stop('radius must be an integer >= 1.')
	# }
	
	# metric <- match.arg(metric, choices = c('Sorensen', 'RWTurnover', 'phyloRWTurnover'))
	# if (!metric %in% c('Sorensen', 'RWTurnover', 'phyloRWTurnover')) {
		# stop('Invalid metric.')
	# }
	
	# # if phylogenetic, there must be a phylo object
	# # check that there is a phylogeny in speciesRaster object
	# if (metric == 'phyloRWTurnover') {
		# if (is.null(x[['phylo']])) {
			# stop('speciesRaster object does not contain a phylo object!')
		# }
		
		# if (verbose) cat('\t...dropping species that are not in phylo data...\n')
 		# # prune speciesRaster object down to species shared with phylogeny
		# x[[2]] <- intersectList(x[[2]], x[['phylo']]$tip.label)
		
		# spEdges <- getRootToTipEdges(x[['phylo']])
		
		# if (!'edgeArea' %in% names(x)) {
			# if (verbose) cat('\t...calculating branch-specific range sizes...\n')
			# x[['edgeArea']] <- do.call(cbind, phyloBranchRanges(x[['phylo']], convertNAtoEmpty(x[['speciesList']]), spEdges))
		# }
	# } 
	
	# #Identify cells for moving window
	# #create neighbor matrix from window size
	# if (verbose) cat('\t...configuring the moving window...\n')
	# side <- radius * 2 + 1
	# counter <- radius
	# nb <- matrix(rep(1, side * side), nrow = side, ncol = side)
	# nb[ceiling(side / 2), ceiling(side / 2)] <- 0
	
	# # identify non-empty cells
	# nonNAcells <- which(!is.na(raster::values(x[[1]])))

	# # return neighbor cells for each cell	
	# if (verbose) cat('\t...identifying cell neighbors...\n')
	# all_nb <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = 8)
	
	# # remove any cells that have no neighbors (isolated cells)
	# nonNAcells <- intersect(nonNAcells, all_nb[,1])
	
	# # run neighborhood identification again but now with only valid cells
	# all_nb <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = nb)
	
	# # subset species cell list to non-empty cells
	# spCellList <- x[['speciesList']][nonNAcells]

	# # remap the cell numbers to the list with empty cells removed
	# # cellNumVec is a vector that converts original cell numbers (the names) to the new list indices where empty cells have been removed (scientific notation was creating problems)
	# if (verbose) cat('\t...remapping cell neighbors...')
	# scipenVal <- getOption('scipen')
	# options('scipen' = 999)
	# cellNumVec <- rep(NA, raster::ncell(x[[1]]))
	# names(cellNumVec) <- 1:raster::ncell(x[[1]])
	# cellNumVec[nonNAcells] <- 1:length(nonNAcells)
	# all_nb[,2] <- cellNumVec[as.character(all_nb[,2])]
	# nbList <- split(all_nb[,2], all_nb[,1])	
	# options('scipen' = scipenVal)
	# if (verbose) cat('done...\n')
	
	# if (verbose) cat(paste0('\t...calculating metric ', metric, '\n'))
	# if (metric == 'Sorensen') {
		# cellBeta <- calcRWTurnover_taxonomic(spCellList, nbList)
	# } else if (metric == 'RWTurnover') {
		# cellBeta <- calcRWTurnover_rangeWeighted(spCellList, nbList, 1 / x[['cellCount']])
	# } else if (metric == 'phyloRWTurnover') {
		# cellBeta <- calcRWTurnover_phyloRangeWeighted(spCellList, nbList, x[['phylo']], spEdges, x[['edgeArea']])
	# }
	
	# cellVal <- rep(NA, raster::ncell(x[[1]]))
	# cellVal[nonNAcells] <- cellBeta
	# res <- x[[1]]
	# raster::values(res) <- cellVal
	# names(res) <- metric
	# return(res)	
# }



# betaDiversity_speciesRaster <- function(x, radius = 4, metric, nthreads = 1, verbose = TRUE) {
	# # radius is distance in cells
		
	# if (!'speciesRaster' %in% class(x)) {
		# stop('x must be of class speciesRaster.')
	# }

	# # check metric validity
	# if (length(metric) > 1) {
		# stop('Only one metric can be specified.')
	# }
	
	# if (radius < 1 | radius != trunc(radius)) {
		# stop('radius must be an integer >= 1.')
	# }
	
	# metric <- match.arg(metric, choices = c('Sorensen', 'RWTurnover', 'phyloRWTurnover'))
	# if (!metric %in% c('Sorensen', 'RWTurnover', 'phyloRWTurnover')) {
		# stop('Invalid metric.')
	# }
	
	# # if phylogenetic, there must be a phylo object
	# # check that there is a phylogeny in speciesRaster object
	# if (metric == 'phyloRWTurnover') {
		# if (is.null(x[['phylo']])) {
			# stop('speciesRaster object does not contain a phylo object!')
		# }
		
		# if (verbose) cat('\t...dropping species that are not in phylo data...\n')
 		# # prune speciesRaster object down to species shared with phylogeny
		# x[[2]] <- intersectList(x[[2]], x[['phylo']]$tip.label)
		
		# spEdges <- getRootToTipEdges(x[['phylo']])
		
		# if (!'edgeArea' %in% names(x)) {
			# if (verbose) cat('\t...calculating branch-specific range sizes...\n')
			# x[['edgeArea']] <- do.call(cbind, phyloBranchRanges(x[['phylo']], convertNAtoEmpty(x[['speciesList']]), spEdges))
		# }
	# } 
	
	# #Identify cells for moving window
	# #create neighbor matrix from window size
	# if (verbose) cat('\t...configuring the moving window...\n')
	# side <- radius * 2 + 1
	# counter <- radius
	# nb <- matrix(rep(1, side * side), nrow = side, ncol = side)
	# nb[ceiling(side / 2), ceiling(side / 2)] <- 0
	
	# # identify non-empty cells
	# nonNAcells <- which(!is.na(raster::values(x[[1]])))
	
	# # remove cells that have no neighbors
	# nb_firstpass <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = 8)
	# nonNAcells <- intersect(nonNAcells, nb_firstpass[,1])

	# # subset species cell list to non-empty cells
	# spCellList <- x[['speciesList']][nonNAcells]

	# # remap the cell numbers to the list with empty cells removed
	# # cellNumVec is a vector that converts original cell numbers (the names) to the new list indices where empty cells have been removed (scientific notation was creating problems)
	# if (verbose) cat('\t...remapping cell neighbors...')
	# scipenVal <- getOption('scipen')
	# options('scipen' = 999)
	# cellNumVec <- rep(NA, raster::ncell(x[[1]]))
	# names(cellNumVec) <- 1:raster::ncell(x[[1]])
	# cellNumVec[nonNAcells] <- 1:length(nonNAcells)
		
	# if (metric == 'Sorensen') {
		
		# if (nthreads > 1) {
			# cl <- parallel::makePSOCKcluster(nthreads)
			# parallel::clusterExport(cl = cl, varlist = c('x', 'nonNAcells', 'nb', 'spCellList', 'cellNumVec', 'calcRWTurnover_taxonomic_singleCell'), envir = environment())	
		# } else {
			# cl <- 1
		# }
		
		# cellBeta <- pbapply::pbsapply(nonNAcells, function(y) {
			# nbVec <- raster::adjacent(x[[1]], cells = y, target = nonNAcells, directions = nb, pairs = FALSE)
			# nbVec <- cellNumVec[as.character(nbVec)]
			# calcRWTurnover_taxonomic_singleCell(spCellList[[which(nonNAcells == y)]], spCellList[nbVec])
		# }, cl = cl)
		
	# } else if (metric == 'RWTurnover') {
		
		# if (nthreads > 1) {
			# cl <- parallel::makePSOCKcluster(nthreads)
			# parallel::clusterExport(cl = cl, varlist = c('x', 'nonNAcells', 'nb', 'spCellList', 'cellNumVec', 'calcRWTurnover_rangeWeighted_singleCell'), envir = environment())	
		# } else {
			# cl <- 1
		# }
		
		# cellBeta <- pbapply::pbsapply(nonNAcells, function(y) {
			# nbVec <- raster::adjacent(x[[1]], cells = y, target = nonNAcells, directions = nb, pairs = FALSE)
			# nbVec <- cellNumVec[as.character(nbVec)]
			# cellBeta <- calcRWTurnover_rangeWeighted_singleCell(spCellList[[which(nonNAcells == y)]], spCellList[nbVec], 1 / x[['cellCount']])
		# }, cl = cl)
		
	# } else if (metric == 'phyloRWTurnover') {

		# if (nthreads > 1) {
			# cl <- parallel::makePSOCKcluster(nthreads)
			# parallel::clusterExport(cl = cl, varlist = c('x', 'nonNAcells', 'nb', 'spCellList', , 'cellNumVec', 'spEdges', 'calcRWTurnover_phyloRangeWeighted_singleCell'), envir = environment())	
		# } else {
			# cl <- 1
		# }
		
		# cellBeta <- pbapply::pbsapply(nonNAcells, function(y) {
			# nbVec <- raster::adjacent(x[[1]], cells = y, target = nonNAcells, directions = nb, pairs = FALSE)
			# nbVec <- cellNumVec[as.character(nbVec)]
			# cellBeta <- calcRWTurnover_phyloRangeWeighted_singleCell(spCellList[[which(nonNAcells == y)]], spCellList[nbVec], x[['phylo']]$tip.label, spEdges, x[['edgeArea']])
		# }, cl = cl)
	# }

	# if (nthreads > 1) {
		# parallel::stopCluster(cl)
	# }	
	
	# options('scipen' = scipenVal)
	
	# cellVal <- rep(NA, raster::ncell(x[[1]]))
	# cellVal[nonNAcells] <- cellBeta
	# res <- x[[1]]
	# raster::values(res) <- cellVal
	# names(res) <- metric
	# return(res)	
# }

