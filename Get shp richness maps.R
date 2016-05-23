################################################################
#### Erosion of phylogenetic and functional diversity in Amphibians
#### Brunno Oliveira, 2015
#### Universidade Federal do Rio Grande do Norte - Brasil
################################################################


library(sp)
library(rgdal)
library(raster)
library(letsR)
library(maptools)

setwd('F:/Arboreality Fossoriality')
#setwd("~/Arboreality Fossoriality")


### Get presence/absence grid for 2dg cells
# Load the species shapefile: XX.shp.
# Between the "" write the name of the shp file.
data <- readOGR(dsn= "AMPHIBIANS.shp",layer="AMPHIBIANS")

occr <- lets.presab(data, xmn = -180, xmx = 180, ymn = -90, ymx = 90, resol = 2, crs = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

rm(data)

save.image("Get shp richness maps.RData")