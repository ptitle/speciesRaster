##' @title speciesRaster summary
##'
##' @description Generates a summary of a speciesRaster object.
##'
##' @param object object of class \code{speciesRaster}
##' @param ... further arguments passed to \code{\link{summary}}
##'
##' @details
##' Summary information includes
##' 
##' @return A list containing the summary information is returned invisibly.
##'
##' @author Pascal Title
##' 
##' @examples
##' summary(tamiasSpRas)
##' attr <- summary(tamiasSpRas)
##' attr
##' 
##' @rdname summary
##' @export

summary.speciesRaster <- function(object, ...) {
	
	if (!'speciesRaster' %in% class(object)) {
		stop('Object must be of class speciesRaster.')
	}
	
	# if data present in object, then report info
	if (class(object[['data']]) %in% c('numeric', 'matrix', 'data.frame')) {
		if (is.vector(object[['data']])) {
			data <- length(intersect(object[['geogSpecies']], names(object[['data']])))
		} else {
			data <- length(intersect(object[['geogSpecies']], rownames(object[['data']])))
		}
	} else {
		data <- NA
	}
	
	# if phylogeny present in object, then report info
	if (class(object[['phylo']]) %in% 'phylo') {
		phylo <- length(intersect(object[['geogSpecies']], object[['phylo']]$tip.label))
	} else {
		phylo <- NA
	}
	
	metric <- names(object[[1]])
	ncells <- raster::ncell(object[[1]])
	rasterExtent <- raster::extent(object[[1]])
	resolution <- raster::res(object[[1]])
	proj <- sf::st_crs(object[[1]])
	lengthUniqueSp <- length(object[['geogSpecies']])
	minSp <- min(sapply(object[[2]], length))
	maxSp <- max(sapply(object[[2]], length))
	
	cat('\n\tSummary of speciesRaster object:\n\n')
	cat('\tMetric:', metric, '\n')
	cat('\tnumber of raster cells:', ncells, '\n')
	cat('\traster resolution:', resolution[1], 'by', resolution[2], '\n')
	cat('\traster projection:', proj$proj4string, '\n\n')
	cat(paste0('\tnumber of unique species: ', lengthUniqueSp, ' (richness range: ', minSp, ' - ', maxSp, ')'), '\n')
	cat('\tdata present:', ifelse(is.na(data), 'No', 'Yes'), '\n')
	if (!is.na(data)) {
		cat('\tnumber of species shared between data and raster:', data, '\n')
	}

	cat('\tphylogeny present:', ifelse(is.na(phylo), 'No', 'Yes'), '\n')
	if (!is.na(phylo)) {
		cat('\tnumber of species shared between phylogeny and raster:', phylo)
	}
		
	obj <- list(
				ncells = ncells, 
				extent = rasterExtent, 
				resolution = resolution, 
				crs = proj, 
				numberUniqueSpecies = lengthUniqueSp, 
				minSp = minSp,
				maxSp = maxSp,
				overlapWithPhylogeny = phylo)
				
	invisible(obj)
}
