#######################################################
#
#  Mask out waterbodies
#  and visualize maps
#
# GSP-Secretariat
# Contact: Isabel.Luotto@fao.org
#
#######################################################

#Empty environment and cache ----
rm(list = ls());
gc()

#######################################################
#
#  User defined variables:

# Working directory
wd <- 'C:/Users/luottoi/Documents/GitHub/Digital-Soil-Mapping'
#wd <- 'C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping'



# Area of interest
AOI <- '01-Data/MKD.shp'
# GEE Resolution 
res = 250

soilatt <- 'ocs_0_30'

#
#
#######################################################


# Load libraries ----
library(tidyverse)
library(terra)
library(rgee)
library(sf)
library(tmap)


# Set working directory ----
setwd(wd)

# Upload own AOI shapefile ----
AOI <- read_sf(AOI)
# convert AOI to a box polygon
#AOI <- st_as_sfc(st_bbox(AOI))
#AOI <- st_as_sf(AOI)


#Initialize GEE ----
ee_Initialize()

#Initial setup either convert shp to gee geometry ----

#Convert shp to gee geometry
region <- sf_as_ee(AOI)
region <- region$geometry()


# Water bodies mask  ----
# Go to https://code.earthengine.google.com/
# browse for the dataset you're interested in
# find and copy/paste the path  

image1 <- ee$Image("JRC/GSW1_0/GlobalSurfaceWater") %>%
  ee$Image$select("occurrence")%>%
  ee$Image$clip(region)%>%
  ee$Image$toDouble()


proj <- image1$projection()$getInfo()

crs <- proj$crs

image1 = image1$resample('bilinear')$reproject(
  crs= crs,
  scale= res)


# CTRL + shift + C to comment
wtbs <- ee_as_raster(
  image = image1,
  scale= res,
  region = region,
  via = "drive"
)

wtbs <-rast(wtbs)
plot(wtbs)
wtbs[wtbs$occurrence <=90] <-NaN
wtbs[wtbs$occurrence >=90] <-1


#Upload final maps and mask out waterbodies
pred_mean <- rast(paste0('02-Outputs/maps/',soilatt,'_QRF.tif'))
pred_sd <- rast(paste0('02-Outputs/maps/',soilatt,'_QRF_SD.tif'))


pred_mean <- mask(pred_mean, wtbs, inverse=T)
pred_sd <- mask(pred_sd, wtbs, inverse=T)

plot(pred_mean)
plot(pred_sd)

#Export masked layers
writeRaster(pred_mean, paste0('02-Outputs/maps/',soilatt,'_QRF_mask.tif'), overwrite=T)
writeRaster(pred_sd, paste0('02-Outputs/maps/',soilatt,'_QRF_SD_mask.tif'), overwrite=T)

# 3 - Zonal statistics =========================================================

zones <- vect('01-Data/soil map/SoilTypes.shp')
plot(zones, "Symbol")

x <- terra::extract(pred_mean, zones)
zones$ocs_0_30 <- x$lyr1

plot(zones, "ocs_0_30")

writeVector(zones, '02-Outputs/maps/zones.shp', overwrite = TRUE)

# 3.1 - Plot zonal stat

(zones <- as_tibble(zones))

unique(select(zones, FAO, Symbol))

ggplot(zones, aes(x = Symbol, y = ocs_0_30, fill = Symbol)) + 
  geom_boxplot() + theme(legend.position = "") +
  ylab("Organic Carbon Stock 0-30cm tn/ha") + 
  xlab("Cartographic Symbol") + 
  coord_flip()

# Plot maps 
pred <- c(pred_mean, pred_sd)
names(pred) <- c("predicted mean", "predicted sd")
tm_shape(pred) +
  tm_raster(title = "OCS 0-30cm \ntn/ha", style = "cont") +
  tm_layout(legend.outside=TRUE, legend.just = c("left", "top"), 
            legend.width=-0.17, legend.height=-0.3,
            inner.margins = c(0,0,0,0)) 
  
tm_shape(pred_mean, style = "natural") +
  tm_raster(title = "OCS 0-30cm \ntn/ha", style = "cont") +
  tm_compass(type = "rose", size = 1.5, color.dark = "#656565", 
             position = c("RIGHT", "TOP")) +
  tm_scale_bar(breaks = c(0,50), color.dark = "#656565",
               position = c("RIGHT", "BOTTOM")) +
  tm_layout(legend.outside=TRUE, legend.just = c("left", "top"), 
            legend.width=-0.17, legend.height=-0.3,
            inner.margins = c(0,0,0,0)) 

