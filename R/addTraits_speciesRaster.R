##' @title addTraits_speciesRaster
##'
##' @description Add univariate or multivariate trait data to a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param data named numeric vector, matrix or dataframe with rownames corresponding to species in \code{x}. 
##' @param replace boolean; if data is already a part of \code{x},
##' should it be replaced?
##'
##' @details If any species in \code{data} are not found in the speciesRaster
##' geographical data, then those species will be dropped from \code{data}, and
##' a warning will be issued. 
##'
##' @return object of class \code{speciesRaster}, with trait data
##' as the list element named \code{data}. 
##'
##' @author Pascal Title
##' 
##' @examples
##' tamiasSpRas
##' tamiasTraits
##'
##' addTraits_speciesRaster(tamiasSpRas, tamiasTraits)
##' 
##' @export


addTraits_speciesRaster <- function(x, data, replace = FALSE) {
	
	if (!'speciesRaster' %in% class(x)) {
		stop('x must be of class speciesRaster.')
	}
	
	if (!class(data) %in% c('numeric', 'matrix', 'data.frame')) {
		stop('data must be either a numeric vector, matrix or dataframe.')
	}
	
	if (class(x[['data']]) %in% c('numeric', 'matrix', 'data.frame') & !replace) {
		stop('Data already present. If data are to be replaced, set replace = TRUE')
	}
	
	# drop species from trait vector if missing from raster
	if (is.vector(data)) {
		traitSpecies <- intersect(x$geogSpecies, names(data))
		inGeogNotData <- setdiff(x$geogSpecies, names(data))
		inDataNotGeog <- setdiff(names(data), x$geogSpecies)
		x[['data']] <- data[traitSpecies]
	} else if (class(data) %in% c('matrix','data.frame')) {
		traitSpecies <- intersect(x$geogSpecies, rownames(data))
		inGeogNotData <- setdiff(x$geogSpecies, rownames(data))
		inDataNotGeog <- setdiff(rownames(data), x$geogSpecies)
		x[['data']] <- as.matrix(data[traitSpecies,])
	}
	
	if (length(inDataNotGeog) > 0) {
		cat('Warning: The following species were dropped from the trait data because they lacked geographic data:\n')
		for (i in 1:length(inDataNotGeog)) {
			cat('\t', inDataNotGeog[i], '\n')
		}
	}
	
	return(x)
}

