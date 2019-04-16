##' @title Map turnover in species communities v2
##'
##' @description Mean community dissimilarity in terms of species composition, phylogeny or traits is 
##' 	calculated for each cell within a circular moving window of neighboring cells. 
##' 
##' @param x object of class \code{speciesRaster}.
##' @param radius Radius of the moving window in map units.
##' @param metric choice of metric, see details.
##' @param component which component of beta diversity to use, can be \code{"turnover"}, 
##' 	\code{"nestedness"} or \code{"full"}
##'
##' @details
##' 	For each cell, mean dissimilarity is calculated from the focal cell to each of its neighbors.
##' 
##' 	The following metrics are available. All metrics are based on Sorensen dissimilarity and range 
##' 	from 0 to 1:
##' 	\itemize{
##' 		\item{taxonomic}
##'			\item{phylogenetic}
##'			\item{trait}
##'		}
##'		For each metric, the following components can be specified. These components are additive,
##' 	such that the full metric = turnover + nestedness. 
##' 	\itemize{
##'			\item{turnover}: species turnover without the influence of richness differences
##'			\item{nestedness}: species turnover due to differences in richness
##'			\item{full}: the combined turnover due to both differences in richness and pure turnover
##'		}
##'
##' Trait turnover \code{trait} measures the turnover in community trait data. This metric is 
##' identical to phylogenetic turnover, but where the phylogeny is replaced with a neighbor-joining tree 
##' that is generated from the trait distance matrix.
##'
##'	The following are currently not implemented, but might be in the future:
##' Range-weighted turnover \code{RWTurnover} measures turnover but where taxa are weighted
##' 	according to the inverse of their range size.
##' 
##' Phylogenetic range-weighted turnover \code{phyloRWTurnover} measures turnover in phylogenetic diversity
##' 	where phylogenetic branches are weighted by the inverse of their geographic distribution.
##'
##' 
##' 
##' @return Returns a raster with mean community dissimilarity for each cell.
##' 
##' @author Pascal Title
##'
##' @references
##' 
##' Baselga, A. The relationship between species replacement, dissimilarity derived from nestedness, 
##' and nestedness. Global Ecology and Biogeography 21 (2012): 1223–1232.
##' 
##' Laffan, SW, et al. Range-weighted metrics of species and phylogenetic turnover can better 
##' resolve biogeographic transition zones. Methods in Ecology and Evolution 7 (2016): 580-588.
##'
##' Leprieur, F, Albouy, C, De Bortoli, J, Cowman, PF, Bellwood, DR & Mouillot, D. Quantifying 
##' Phylogenetic Beta Diversity: Distinguishing between "True" Turnover of Lineages and Phylogenetic 
##' Diversity Gradients. PLoS ONE 7 (2012): e42760–12.
##'
##' Rosauer, D, Laffan, SW, Crisp, MD, Donnellan, SC, Cook, LG. Phylogenetic endemism: a new approach 
##' for identifying geographical concentrations of evolutionary history. Molecular Ecology 
##' 18 (2009): 4061-4072.
##' 
##' @examples
##' \donttest{
##' library(raster)
##' tamiasSpRas
##'
##' tamiasSpRas <- addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' tamiasSpRas <- addTraits_speciesRaster(tamiasSpRas, tamiasTraits)
##' 
##' # taxonomic turnover
##' beta_taxonomic_turnover <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'taxonomic', component = 'turnover')
##' beta_taxonomic_nestedness <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##'		metric = 'taxonomic', component = 'nestedness')
##' beta_taxonomic_full <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##'		metric = 'taxonomic', component = 'full')
##'
##' # phylogenetic turnover
##' beta_phylo_turnover <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'phylogenetic', component = 'turnover')
##' beta_phylo_nestedness <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'phylogenetic', component = 'nestedness')
##' beta_phylo_full <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'phylogenetic', component = 'full')
##'
##' # trait turnover
##' beta_trait_turnover <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'trait', component = 'turnover')
##' beta_trait_nestedness <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'trait', component = 'nestedness')
##' beta_trait_full <- betaDiversity2_speciesRaster(tamiasSpRas, radius = 70000, 
##' 	metric = 'trait', component = 'full')
##' 
##' 
##' colramp <- colorRampPalette(c('blue','yellow','red'))
##' par(mfrow=c(1,3))
##' plot(beta_taxonomic_full, col = colramp(100), zlim = c(0,1))
##' plot(beta_phylo_full, col = colramp(100), zlim = c(0,1))
##' plot(beta_trait_full, col = colramp(100), zlim = c(0,1))
##' 
##' }
##' @export



betaDiversity2_speciesRaster <- function(x, radius, metric, component = 'full') {
	# radius is in map units
		
	if (!inherits(x, 'speciesRaster')) {
		stop('x must be of class speciesRaster.')
	}

	# check metric validity
	if (length(metric) > 1) {
		stop('Only one metric can be specified.')
	}
	
	if (radius < max(raster::res(x[[1]]))) {
		stop(paste0('The radius must at minimum be ', max(raster::res(x[[1]])), '.'))
	}
	
	metric <- match.arg(metric, choices = c('taxonomic', 'phylogenetic', 'trait'))
	if (!metric %in% c('taxonomic', 'phylogenetic', 'trait')) {
		stop('Invalid metric.')
	}
	
	component <- match.arg(component, choices = c('turnover', 'nestedness', 'full'))
	if (!component %in% c('turnover', 'nestedness', 'full')) {
		stop('component must be turnover, nestedness or full.')
	}
	
	if (metric == 'phylogenetic') {
		if (is.null(x[['phylo']])) {
			stop('speciesRaster object does not contain a phylo object!')
		}
		
		# prune speciesRaster object down to species shared with phylogeny
		if (!identical(sort(x[['geogSpecies']]), sort(x[['phylo']]$tip.label))) {
			x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']]$tip.label)
		}
	}
	
	if (metric == 'trait') {
		if (is.null(x[['data']])) {
			stop('speciesRaster object does not contain any associated trait data!')
		}
		
		# prune speciesRaster object down to species shared with data
		if (is.vector(x[['data']])) {
			if (!identical(sort(x[['geogSpecies']]), sort(names(x[['data']])))) {
				x[['speciesList']] <- intersectList(x[['speciesList']], names(x[['data']]))
			}
		} else {
			if (!identical(sort(x[['geogSpecies']]), sort(rownames(x[['data']])))) {
				x[['speciesList']] <- intersectList(x[['speciesList']], rownames(x[['data']]))
			}
		}
		
		
		# calculate pairwise distances
		d <- dist(x[['data']])
		
		# generate tree structure from distance matrix
		traitTree <- ape::bionj(d)	
	}
	
	# convert NA communities to "empty" (easier to handle in Rcpp)
	isNAind <- which(sapply(x[['speciesList']], anyNA) == TRUE)
	if (length(isNAind) > 0) {
		for (i in 1:length(isNAind)) {
			x[['speciesList']][[isNAind[i]]] <- 'empty'
		}
	}	
		
	if (metric == 'taxonomic') {
		pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = component)
	} else if (metric == 'phylogenetic') {
		pairwiseD <- calcPairwisePhylosor2(x[['speciesList']], x[['phylo']], component = component)
	} else if (metric == 'trait') {
		pairwiseD <- calcPairwisePhylosor2(x[['speciesList']], traitTree, component = component)
	}
	
	pairwiseD[pairwiseD < 0] <- NA

	wMat <- raster::focalWeight(x[[1]], d = radius, type = 'circle')
	wMat <- wMat > 0
	mode(wMat) <- 'numeric'
	wMat[wMat == 0] <- NA
	wMat[mean(1:nrow(wMat)), mean(1:ncol(wMat))] <- 0
	
	# identify non-empty cells
	nonNAcells <- which(!is.na(raster::values(x[[1]])))

	# return neighbor cells for each cell	
	all_nb <- raster::adjacent(x[[1]], cells = nonNAcells, target = nonNAcells, directions = wMat)
	
	# remove any cells that have no neighbors (isolated cells)
	nonNAcells <- intersect(nonNAcells, all_nb[,1])
	
	all_nb <- split(all_nb[,2], all_nb[,1])
	
	cellVals <- rep(NA, raster::ncell(x[[1]]))
	for (i in 1:length(nonNAcells)) {
		# cat(i, '\n')
		focalCell <- nonNAcells[i]
		nbCells <- all_nb[[i]]
		
		# convert to cell indices
		focalCell <- x[['cellCommInd']][focalCell]
		nbCells <- x[['cellCommInd']][nbCells]
		
		cellVals[nonNAcells[i]] <- mean(pairwiseD[focalCell, nbCells], na.rm=TRUE)		
	}

	ret <- raster::raster(x[[1]])
	raster::values(ret) <- cellVals
	names(ret) <- metric
	return(ret)
}

## TESTING ##
# verify distances with betapart package to ensure calculations are correct

# # sourceCpp('~/Dropbox/speciesRaster/src/betaDiversity.cpp')
# testWithBetaPart <- function(x, metric, component) {
	# if (!metric %in% c('taxonomic', 'phylogenetic') | !component %in% c('turnover', 'nestedness', 'full')) {
		# stop()
	# }
	
	# commMat <- matrix(0, nrow = length(x[['speciesList']]), ncol = length(na.omit(unique(unlist(x[['speciesList']])))))
	# colnames(commMat) <- na.omit(unique(unlist(x[['speciesList']])))
	# for (i in 1:length(x[['speciesList']])) {
		# if (!anyNA(x[['speciesList']][[i]])) {
			# commMat[i,] <- as.numeric(colnames(commMat) %in% x[['speciesList']][[i]])
		# }
	# }
	# if ('empty' %in% colnames(commMat)) {
		# commMat[, 'empty'] <- 0
	# }
	
	# emptyInd <- which(rowSums(commMat) == 0)
	
	# if (metric == 'phylogenetic') {
		# beta.mat <- betapart::phylo.beta.pair(commMat, x[['phylo']])
		# if (component == 'turnover') {
			# beta.mat <- as.matrix(beta.mat$phylo.beta.sim)
		# }
		# if (component == 'nestedness') {
			# beta.mat <- as.matrix(beta.mat$phylo.beta.sne)
		# }
		# if (component == 'full') {
			# beta.mat <- as.matrix(beta.mat$phylo.beta.sor)
		# }
	# } else if (metric == 'taxonomic') {
		# beta.mat <- betapart::beta.pair(commMat)
		# if (component == 'turnover') {
			# beta.mat <- as.matrix(beta.mat$beta.sim)
		# }
		# if (component == 'nestedness') {
			# beta.mat <- as.matrix(beta.mat$beta.sne)
		# }
		# if (component == 'full') {
			# beta.mat <- as.matrix(beta.mat$beta.sor)
		# }			
	# }
	
	# dtest <- matrix(nrow = length(x[['speciesList']]), ncol = length(x[['speciesList']]))
	# for (i in 1:length(x[['speciesList']])) {
		# for (j in 1:length(x[['speciesList']])) {
			# if (i >= j) {
				# dtest[i,j] <- beta.mat[i,j]
				# dtest[j,i] <- dtest[i,j]
			# }
		# }
	# }
	# dtest[is.nan(dtest)] <- NA
	# dtest[emptyInd,] <- NA
	# dtest[, emptyInd] <- NA
	# return(dtest)
# }



# # taxonomic

# pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = 'turnover'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'taxonomic', component = 'turnover')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)

# pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = 'nestedness'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'taxonomic', component = 'nestedness')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)

# pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = 'full'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'taxonomic', component = 'full')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)



# # phylogenetic
# pairwiseD <- calcPairwisePhylosor2(x[['speciesList']], x[['phylo']], component ='turnover'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'phylogenetic', component = 'turnover')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)
# identical(round(pairwiseD, 6), round(betapartD, 6))

# pairwiseD <- calcPairwisePhylosor(x[['speciesList']], x[['phylo']], component ='nestedness'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'phylogenetic', component = 'nestedness')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)
# identical(round(pairwiseD, 6), round(betapartD, 6))

# pairwiseD <- calcPairwisePhylosor(x[['speciesList']], x[['phylo']], component ='full'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'phylogenetic', component = 'full')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)
# identical(round(pairwiseD, 6), round(betapartD, 6))




# test with morphdist

# # morphDist <- function(x) {
	# dtest <- matrix(nrow = length(x[['speciesList']]), ncol = length(x[['speciesList']]))

	# d <- as.matrix(dist(x[['data']]))
	# diag(d) <- NA

	# for (i in 1:length(x[['speciesList']])) {
		# for (j in 1:length(x[['speciesList']])) {
			# comm1 <- x[['speciesList']][[i]]
			# comm2 <- x[['speciesList']][[j]]
			# if (i >= j) {
				# if (!anyNA(c(comm1, comm2)) & !('empty' %in% c(comm1, comm2))) {
					# if (identical(comm1, comm2)) {
						# dtest[i,j] <- 0
						# dtest[j,i] <- 0
					# } else {
						# dtest[i,j] <- mean(d[comm1, comm2], na.rm=TRUE)			
						# dtest[j,i] <- dtest[i,j]
					# }
				# }
			# }
		# }
	# }
	# return(dtest)
# }


