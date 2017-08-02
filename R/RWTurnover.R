##' @title Map turnover in beta diversity
##'
##' @description Mean community dissimilarity is calculated for each cell within a moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##' @param radius Radius in terms of number of cells in all directions from the focal cell 
##' 	to define the moving window. 
##' @param metric choice of metric, see details.
##' @param verbose Primarily intended for debugging, print progress to the console.
##'
##' @details
##' Simpson's Beta diversity \code{SimpsonBeta} is a purely taxonomic measure of turnover.
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
##' beta_taxonomic <- betaDiversity_speciesRaster(tamiasSpRas, radius = 4, metric = 'SimpsonBeta')
##' 
##' beta_RW <- betaDiversity_speciesRaster(tamiasSpRas, radius = 4, metric = 'RWTurnover')
##' 
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' beta_phyloRW <- betaDiversity_speciesRaster(tamiasSpRas, radius = 4, metric = 'phyloRWTurnover')
##' 
##' colramp <- colorRampPalette(c('blue','yellow','red'))
##' par(mfrow=c(1,3))
##' plot(beta_taxonomic, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
##' plot(beta_RW, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
##' plot(beta_phyloRW, col = colramp(100), breaks=seq(0, 1,length.out=99), legend = FALSE)
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
	
	metric <- match.arg(metric, choices = c('SimpsonBeta', 'RWTurnover', 'phyloRWTurnover'))
	if (!metric %in% c('SimpsonBeta', 'RWTurnover', 'phyloRWTurnover')) {
		stop('Invalid metric.')
	}
	
	# if phylogenetic, there must be a phylo object
	# check that there is a phylogeny in speciesRaster object
	if (metric == 'phyloRWTurnover') {
		if (is.null(x[['phylo']])) {
			stop('speciesRaster object does not contain a phylo object!')
		}
		
		if (verbose) cat('\t...dropping species that are not in phylo data...\n')
 		# prune speciesRaster object down to species shared with phylogeny
		x[[2]] <- intersectList(x[[2]], x[['phylo']]$tip.label)
		
		spEdges <- getRootToTipEdges(x[['phylo']])
		
		if (!'edgeArea' %in% names(x)) {
			if (verbose) cat('\t...calculating branch-specific range sizes...\n')
			x[['edgeArea']] <- do.call(cbind, phyloBranchRanges(x[['phylo']], convertNAtoEmpty(x[['speciesList']]), spEdges))
		}
	} 
	
	#Identify cells for moving window
	#create neighbor matrix from window size
	if (verbose) cat('\t...configuring the moving window...\n')
	side <- radius * 2 + 1
	counter <- radius
	nb <- matrix(rep(1, side * side), nrow = side, ncol = side)
	nb[ceiling(side / 2), ceiling(side / 2)] <- 0
	
	# identify non-empty cells
	nonNAcells <- which(!is.na(raster::values(x[[1]])))

	#return neighbor cells for each cell	
	if (verbose) cat('\t...identifying cell neighbors...\n')
	all_nb <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = nb)
	
	# remove any cells that have no neighbors (isolated cells)
	nonNAcells <- intersect(nonNAcells, all_nb[,1])
	
	# run neighborhood identification again but now with only valid cells
	all_nb <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = nb)
	
	# subset species cell list to non-empty cells
	spCellList <- x[['speciesList']][nonNAcells]

	# remap the cell numbers to the list with empty cells removed
	# cellNumVec is a vector that converts original cell numbers (the names) to the new list indices where empty cells have been removed (scientific notation was creating problems)
	if (verbose) cat('\t...remapping cell neighbors...')
	scipenVal <- getOption('scipen')
	options('scipen' = 999)
	cellNumVec <- rep(NA, raster::ncell(x[[1]]))
	names(cellNumVec) <- 1:raster::ncell(x[[1]])
	cellNumVec[nonNAcells] <- 1:length(nonNAcells)
	all_nb[,2] <- cellNumVec[as.character(all_nb[,2])]
	nbList <- split(all_nb[,2], all_nb[,1])	
	options('scipen' = scipenVal)
	if (verbose) cat('done...\n')
	
	if (verbose) cat(paste0('\t...calculating metric ', metric, '\n'))
	if (metric == 'SimpsonBeta') {
		cellBeta <- calcRWTurnover_taxonomic(spCellList, nbList)
	} else if (metric == 'RWTurnover') {
		cellBeta <- calcRWTurnover_rangeWeighted(spCellList, nbList, 1 / x[['cellCount']])
	} else if (metric == 'phyloRWTurnover') {
		cellBeta <- calcRWTurnover_phyloRangeWeighted(spCellList, nbList, x[['phylo']], spEdges, x[['edgeArea']])
	}
	
	cellVal <- rep(NA, raster::ncell(x[[1]]))
	cellVal[nonNAcells] <- cellBeta
	res <- x[[1]]
	raster::values(res) <- cellVal
	names(res) <- metric
	return(res)	
}




