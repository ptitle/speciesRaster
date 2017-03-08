##' @title addPhylo_speciesRaster
##'
##' @description Add a phylogeny to speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param tree a phylogeny of class \code{phylo}
##' @param replace boolean; if a tree is already a part of \code{x},
##' should it be replaced?
##'
##' @details If any species in the phylogeny are not found in the speciesRaster
##' geographical data, then those species will be dropped from the phylogeny, and
##' a warning will be issued. 
##'
##' @return object of class \code{speciesRaster}, with a \code{phylo}
##' object as the list element named \code{phylo}. 
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasSpRas
##' tamiasTree
##'
##' addPhylo_speciesRaster(tamiasSpRas, tamiasTree)
##' 
##' @export

addPhylo_speciesRaster <- function(x, tree, replace = FALSE) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	if (!'phylo' %in% class(tree)) {
		stop('tree must be a phylo object.')
	}
	
	if (class(x[[5]]) %in% 'phylo' & !replace) {
		stop('Phylogeny already present. If phylogeny is to be replaced, set replace = TRUE')
	}
	
	# if needed, prune species in phylogeny down to species with geog data
	inPhyloNotGeog <- setdiff(tree$tip.label, x[[3]])
	inGeogNotPhylo <- setdiff(x[[3]], tree$tip.label)
	tree <- ape::drop.tip(tree, inPhyloNotGeog)
	x[['phylo']] <- tree
	
	if (length(inGeogNotPhylo) > 0) {
		cat('Warning: The following species were were pruned from the phylogeny because they lack geographic data:\n')
		for (i in 1:length(inGeogNotPhylo)) {
			cat('\t', inGeogNotPhylo[i], '\n')
		}
	}
	
	return(x)
}

