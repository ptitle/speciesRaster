##' @export

print.speciesRaster <- function(x, ...) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('Object must be of class speciesRaster.')
	}
	
	# if data present in object, then report info
	if (class(x[[4]]) %in% c('numeric', 'matrix', 'data.frame')) {
		if (is.vector(x[[4]])) {
			data <- length(intersect(x[[3]], names(x[[4]])))
		} else {
			data <- length(intersect(x[[3]], rownames(x[[4]])))
		}
	} else {
		data <- NA
	}
	
	# if phylogeny present in object, then report info
	if (class(x[[5]]) %in% 'phylo') {
		phylo <- length(intersect(x[[3]], x[[5]]$tip.label))
	} else {
		phylo <- NA
	}
	
	metric <- names(x[[1]])
	ncells <- raster::ncell(x[[1]])
	rasterExtent <- raster::extent(x[[1]])
	resolution <- raster::res(x[[1]])
	proj <- raster::projection(x[[1]])
	lengthUniqueSp <- length(x[[3]])
	minSp <- min(sapply(x[[2]], length))
	maxSp <- max(sapply(x[[2]], length))
	
	cat('\n\tSummary of speciesRaster object:\n\n')
	cat('\tMetric:', metric, '\n')
	cat('\tnumber of raster cells:', ncells, '\n')
	cat('\traster resolution:', resolution[1], 'by', resolution[2], '\n')
	cat('\traster projection:', proj, '\n\n')
	cat(paste0('\tnumber of unique species: ', lengthUniqueSp, ' (richness range: ', minSp, ' - ', maxSp, ')'), '\n')
	cat('\tdata present:', ifelse(is.na(data), 'No', 'Yes'), '\n')
	if (!is.na(data)) {
		cat('\tnumber of species shared between data and raster:', data)
	}

	cat('\tphylogeny present:', ifelse(is.na(phylo), 'No', 'Yes'), '\n')
	if (!is.na(phylo)) {
		cat('\tnumber of species shared between phylogeny and raster:', phylo)
	}
	
}