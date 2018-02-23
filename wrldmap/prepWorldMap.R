# prep world quick map
# from: http://www.naturalearthdata.com/downloads/110m-cultural-vectors/110m-admin-0-countries/

require(cleangeo)
require(maptools)
require(rgeos)

worldmap <- readShapeSpatial('~/Dropbox/speciesRaster/wrldmap/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp', force_ring=T, delete_null_obj=T, repair=T, proj4string=CRS('+proj=longlat +datum=WGS84'))


# worldmap <- clgeo_Clean(worldmap, verbose=TRUE)

clgeo_IsValid(worldmap)
clgeo_SummaryReport(clgeo_CollectionReport(worldmap))

# simplify
worldmap <- gSimplify(worldmap, tol=0.1)

plot(worldmap, lwd=0.5)

# save as sysdata.rda

save(worldmap, file='~/Dropbox/speciesRaster/R/sysdata.rda')

tools::checkRdaFiles('~/Dropbox/speciesRaster/R/sysdata.rda')
tools::resaveRdaFiles('~/Dropbox/speciesRaster/R/sysdata.rda', compress='auto')
tools::checkRdaFiles('~/Dropbox/speciesRaster/R/sysdata.rda')


