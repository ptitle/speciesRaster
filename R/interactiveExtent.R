##' @title Interactively choose extent
##'
##' @description Given a list of polygons, sets up an interactive plot
##' 	to allow the user to draw the desired extent.
##'
##' @param polyList a list of Simple Feature polygons.
##'
##' @return A list with a polygon, and its WKT string 
##'
##' @author Pascal Title
##' 
##' 




interactiveExtent <- function(polyList) {
	
	#get overall extent
	masterExtent <- getExtentOfList(polyList)
	masterExtent <- list(minLong = masterExtent@xmin, maxLong = masterExtent@xmax, minLat = masterExtent@ymin, maxLat = masterExtent@ymax)
	
	# coarse template
	# if projected, use 100km, if not, use 20 degrees
	quickRes <- ifelse(sf::st_is_longlat(polyList[[1]]), 20, 100000)
	proj <- sf::st_crs(polyList[[1]])
	quickTemplate <- raster::raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(quickRes, 2), crs = proj$proj4string)
	quick <- lapply(polyList, function(x) fasterize::fasterize(x, quickTemplate))
	rich <- raster::calc(raster::stack(quick), fun=sum, na.rm = TRUE)
	
	# add map for context
	if (!sf::st_is_longlat(polyList[[1]])) {
		wrld <- sf::st_transform(worldmap, crs = sf::st_crs(polyList[[1]]))
	} else {
		wrld <- worldmap
	}
	
	raster::plot(rich, legend = FALSE)
	message('\n\tAn interactive coarse-grain map has been displayed.\n')
	message('\n\tPlease wait until plot is completed......')
	
	graphics::plot(wrld, add = TRUE, lwd = 0.5)
	message('done!\n')
	graphics::title(main = 'Define your extent polygon.')

	message('\tClick on the map to create a polygon that will define the extent of the rasterStack.')
	message('\tRight-clicking will close the polygon and terminate the interactive plot.\n\n')
	
	userPoly <- raster::drawPoly(sp = TRUE, col = 'red', xpd = NA)
	userPoly <- sf::st_as_sf(userPoly)
	sf::st_crs(userPoly) <- sf::st_crs(polyList[[1]])
	
	# display call so user can use this extent in the future
	wkt <- sf::st_as_text(sf::st_geometry(userPoly))
	
	grDevices::dev.off()
	
	return(list(poly = userPoly, wkt = wkt))
}