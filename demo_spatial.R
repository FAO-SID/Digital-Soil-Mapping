# Isabel Luotto 14/06/22
# Demo spatial data

# Empty the global env. and clean the cache
rm(list=ls())
gc()

# setting the work environment ----
setwd('C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping/')

# load packages
library(raster)

# Load in tif file
landcover <- raster("01-Data/land cover/LandCover.tif")
#Let's check the CRS and extent of 
#our landcover map
landcover 
plot(landcover)#Let's plot our lc
summary(landcover)

DEM <- raster("demonstrations/QGIS Introduction/DEMENV5.tif")
#You can customize your map similarly 
#like we did for the graphs
plot(DEM, col = bpy.colors(255), 
     main= "Digital elevation model", 
     xlab = "Longitude", ylab = "Latitude")

#let's try stacking our landcover and DEM map
lc_and_dem <-stack(landcover, DEM)

#Match the extent and resolution to DEM
landcover <- projectRaster(landcover, DEM, method = 'ngb')

lc_and_dem <-stack(landcover, DEM)
plot(lc_and_dem)
#Now you can plot two 
#rasters at the same time

writeRaster(landcover, 
            file= "demonstrations/QGIS Introduction/example.tif", overwrite=TRUE)

