#######################################################
#
#  Mask our waterbodies
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
# GEE Resolution (CRS defined based on the first TerraClimate layer WGS84 )
res = 1000

soilatt <- 'OCS'

#
#
#######################################################


# Load libraries ----
library(data.table)
library(terra)
library(rgee)
library(sf)

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
region = region$geometry()


# Water bodies mask  ----
# Go to https://code.earthengine.google.com/
# browse for the dataset you're interested in
# find and copy/paste the path  

image1 <- ee$Image("JRC/GSW1_0/GlobalSurfaceWater") %>%
  ee$Image$select("occurrence")%>%
  ee$Image$clip(region)%>%
  ee$Image$toDouble()


proj = image1$projection()$getInfo()

crs =proj$crs

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
wtbs[wtbs$occurrence <=90] <-0
wtbs[wtbs$occurrence >90] <-1
wtbs[wtbs$occurrence ==0] <-NA


#Upload final maps and mask out waterbodies
mean <- rast(paste0(wd,'/02-Outputs/maps/',soilatt,'_QRF.tif'))
sd <- rast(paste0(wd,'/02-Outputs/maps/',soilatt,'_QRF_SD.tif'))


mean <- mask(mean, wtbs, inverse=T)
sd <- mask(sd, wtbs, inverse=T)


#Export masked layers
writeRaster(mean, paste0(wd,'/02-Outputs/maps/',soilatt,'_QRF_mask.tif'), overwrite=T)
writeRaster(sd, paste0(wd,'/02-Outputs/maps/',soilatt,'_QRF_SD_mask.tif'), overwrite=T)