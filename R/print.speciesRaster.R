##' @export

print.speciesRaster <- function(x, ...) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('Object must be of class speciesRaster.')
	}
	
	# if data present in object, then report info
	if (class(x[['data']]) %in% c('numeric', 'matrix', 'data.frame')) {
		if (is.vector(x[['data']])) {
			data <- length(intersect(x[['geogSpecies']], names(x[['data']])))
		} else {
			data <- length(intersect(x[['geogSpecies']], rownames(x[['data']])))
		}
	} else {
		data <- NA
	}
	
	# if phylogeny present in object, then report info
	if (class(x[['phylo']]) %in% 'phylo') {
		phylo <- length(intersect(x[['geogSpecies']], x[['phylo']]$tip.label))
	} else {
		phylo <- NA
	}
	
	metric <- names(x[[1]])
	ncells <- raster::ncell(x[[1]])
	rasterExtent <- raster::extent(x[[1]])
	resolution <- raster::res(x[[1]])
	proj <- sf::st_crs(x[[1]])
	lengthUniqueSp <- length(x[['geogSpecies']])
	minSp <- min(sapply(x[['speciesList']], length))
	maxSp <- max(sapply(x[['speciesList']], length))
	
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
	
	cat('\n')
	
}