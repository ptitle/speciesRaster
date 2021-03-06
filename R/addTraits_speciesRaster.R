##' @title addTraits_speciesRaster
##'
##' @description Add univariate or multivariate trait data to a speciesRaster object.
##'
##' @param x object of class \code{speciesRaster}
##' @param data named numeric vector, matrix or dataframe with rownames corresponding to species in \code{x}
##' 	or pairwise matrix with row and column names corresponding to species in \code{x}. If pairwise matrix,
##' 	the upper triangle of the matrix will be used for calculations.
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
	
	if (!inherits(x, 'speciesRaster')) {
		stop('x must be of class speciesRaster.')
	}
	
	if (!inherits(data, c('numeric', 'matrix', 'data.frame'))) {
		stop('data must be either a numeric vector, matrix or dataframe.')
	}
	
	if (inherits(x[['data']], c('numeric', 'matrix', 'data.frame')) & !replace) {
		stop('Data already present. If data are to be replaced, set replace = TRUE')
	}
	
	# drop species from trait vector if missing from raster
	if (is.vector(data)) {
		if (is.null(names(data))) {
			stop('Data must have names.')
		}
		traitSpecies <- intersect(x$geogSpecies, names(data))
		inGeogNotData <- setdiff(x$geogSpecies, names(data))
		inDataNotGeog <- setdiff(names(data), x$geogSpecies)
		if (length(traitSpecies) == 0) {
			stop('There are no common species in geographic and trait data.')
		}
		x[['data']] <- data[traitSpecies]
	} else if (class(data) %in% c('matrix','data.frame')) {
		if (is.null(rownames(data))) {
			stop('Data must have rownames.')
		}
		traitSpecies <- intersect(x$geogSpecies, rownames(data))
		inGeogNotData <- setdiff(x$geogSpecies, rownames(data))
		inDataNotGeog <- setdiff(rownames(data), x$geogSpecies)
		if (length(traitSpecies) == 0) {
			stop('There are no common species in geographic and trait data.')
		}
		
		if (identical(rownames(data), colnames(data))) {
			# pairwise matrix
			x[['data']] <- as.data.frame(data[traitSpecies, traitSpecies], stringsAsFactors = FALSE)
		} else {
			x[['data']] <- as.data.frame(data[traitSpecies,], stringsAsFactors = FALSE)
		}
	}
	
	if (length(inDataNotGeog) > 0) {
		msg <- paste0('The following species were dropped from the trait data because they lacked geographic data:\n\t', paste(inDataNotGeog, collapse='\n\t'))
		warning(msg)
	}
	
	return(x)
}

