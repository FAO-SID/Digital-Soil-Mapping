#######################################################
#
#  Process and download 22 covariates
#  from GEE to R (10 minutes to execute)
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
# wd <- 'C:/Users/luottoi/Documents/GitHub/Digital-Soil-Mapping'
 wd <- 'C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping'
 
# Area of interest
 AOI <- '01-Data/MKD.shp'
#Start and End time 
 start_T <- "2017-01-01"
 end_T <- "2017-12-31"
 
# Resolution (CRS defined based on the first TerraClimate layer WGS84 )
 res = 1000
#
#
 #######################################################
 
  
# Load libraries ----
library(data.table)
library(terra)
library(sf)
library(rgee)
# Set working directory Module I ----


setwd(wd)

# Country shapefile
AOI <- read_sf(AOI)
# convert AOI to a box polygon
AOI <- st_as_sfc(st_bbox(AOI))
AOI <- st_as_sf(AOI)

#List of covariates to prepare
# Mean annual temperature
# Total annual Precipitation
# Precipitation of wettest month
# Precipitation of driest month

# TAGEE 13 soil attributes (list below)

# MODIS EVI & NDVI

# Daytime temperature SD

# Landsat 8 RED and NIR standard deviation

# OpenLandMap 
# soil water content 0-10-30
# Potential FAPAR Monthly


# Multi-Scale Topographic Position Index




# Mean annual temperature (daytime) ----
ee_Initialize()

region <- sf_as_ee(AOI)
region = region$geometry()


image1 <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate(start_T, end_T) %>%
  ee$ImageCollection$select("tmmx")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

  # from imagecollection to image
image2 <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate(start_T, end_T) %>%
  ee$ImageCollection$select("tmmn")%>%
  ee$ImageCollection$toBands()


image1 <- image1$multiply(0.1)
image2 <- image2$multiply(0.1)

diff <- image1$add(image2)
avT = diff$divide(2)
avT = avT$reduce(ee$Reducer$mean())

proj = avT$projection()$getInfo()

crs = proj$wkt

avT = avT$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

avtr <- ee_as_raster(
  image = avT,
  scale= res,
  region = region,
  via = "drive"
)

writeRaster(avtr, '01-Data/covs/avtr.tif', overwrite=T)

# Total annual Precipitation ----

image <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate(start_T, end_T) %>%
  ee$ImageCollection$select("pr")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$sum()%>%
  ee$Image$toDouble()

Pr = image$resample('bilinear')$reproject(
  crs= crs,
  scale= res)


Prr <- ee_as_raster(
  image = Pr,
  scale= res,
  region = region,
  via = "drive"
)



writeRaster(Prr, '01-Data/covs/Prr.tif', overwrite=T)

# Precipitation of wettest month ----
image <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate(start_T, end_T) %>%
  ee$ImageCollection$select("pr")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

image = image$resample('bilinear')$reproject(
  crs= crs,
  scale= res)


Prr_all <- ee_as_raster(
  image = image,
  scale= res,
  region = region,
  via = "drive"
)

sums <- cellStats(Prr_all, 'sum')
wettest <- names(sums[sums ==max(sums)])
Prr_wet <- Prr_all[[wettest]]

writeRaster(Prr_wet, '01-Data/covs/Prr_wet.tif', overwrite=T)

# Precipitation of driest month ----

dryest <- names(sums[sums ==min(sums)])
Prr_dry <- Prr_all[[dryest]]

writeRaster(Prr_wet, '01-Data/covs/Prr_dry.tif', overwrite=T)


# Terrain attributes using TAGEE - 13 covariates
# Attribute	| Unit	| Description
# Elevation  |	meter	|Height of terrain above sea level
# Slope	degree	| Slope |  gradient
# Aspect |	degree |	Compass direction
# Hillshade	| dimensionless	|Brightness of the illuminated terrain
# Northness |	dimensionless |	Degree of orientation to North
# Eastness	| dimensionless |	Degree of orientation to East
# Horizontal curvature |	meter |	Curvature tangent to the contour line
# Vertical curvature |	meter	| Curvature tangent to the slope line
# Mean curvature |	meter	| Half-sum of the two orthogonal curvatures
# Minimal curvature	meter |	Lowest value of curvature
# Maximal curvature	meter |	Highest value of curvature
# Gaussian curvature	 | meter |	Product of maximal and minimal curvatures
# Shape Index	| dimensionless |	Continuous form of the Gaussian landform classification

# install TAGEE
system("pip install tagee")
# Import
TAGEE <- import("tagee")

image <- ee$Image("MERIT/DEM/v1_0_3") %>%
  ee$Image$clip(region)%>%
  ee$Image$toDouble()

image = image$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

DEMAttributes = TAGEE$terrainAnalysis( image, region)


tagee_test <- ee_as_raster(
  image = DEMAttributes,
  scale= res,
  region = region,
  via = "drive"
)
plot(tagee_test)
writeRaster(tagee_test, '01-Data/covs/Terrain_attributes.tif', overwrite=T)

# EVI & NDVI ----
EVI <- ee$ImageCollection("MODIS/061/MOD13Q1") %>%
  ee$ImageCollection$filterDate(start_T, end_T) %>%
  ee$ImageCollection$select("EVI")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

EVI = EVI$reduce(ee$Reducer$mean())

EVI = EVI$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

NDVI <- ee$ImageCollection("MODIS/061/MOD13Q1") %>%
  ee$ImageCollection$filterDate(start_T, end_T) %>%
  ee$ImageCollection$select("NDVI")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

NDVI = NDVI$reduce(ee$Reducer$mean())

NDVI = NDVI$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

ndvir <- ee_as_raster(
  image = NDVI,
  scale= res,
  region = region,
  via = "drive"
)

evir <- ee_as_raster(
  image = EVI,
  scale= res,
  region = region,
  via = "drive"
)


writeRaster(ndvir, '01-Data/covs/NDVI.tif', overwrite=T)
writeRaster(evir, '01-Data/covs/EVI.tif', overwrite=T)

# Land daytime surface temperature (SD) ----

image <- ee$ImageCollection("MODIS/061/MOD11A1") %>%
  ee$ImageCollection$filterDate("2018-01-01", "2018-12-31") %>%
  ee$ImageCollection$select("LST_Day_1km")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

d_T = image$subtract(273.15)
sd_d_T = d_T$reduce(ee$Reducer$stdDev())

sd_d_T = sd_d_T$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

sd_d_Tr <- ee_as_raster(
  image = sd_d_T,
  scale= res,
  region = region,
  via = "drive"
)


writeRaster(sd_d_Tr, '01-Data/covs/sd_d_Tr.tif', overwrite=T)

# Landsat bands mean and sd ----
image <- ee$ImageCollection("LANDSAT/LC08/C02/T1_RT") %>%
  ee$ImageCollection$filterDate("2018-01-01", "2018-12-31") %>%
  ee$ImageCollection$select(c("B4"))%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()


land_red = image$reduce(ee$Reducer$stdDev())

land_red = land_red$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

land_sd_red <- ee_as_raster(
  image = land_red,
  scale= res,
  region = region,
  via = "drive"
)

writeRaster(land_sd_red, '01-Data/covs/land_sd_red.tif', overwrite=T)

image <- ee$ImageCollection("LANDSAT/LC08/C02/T1_RT") %>%
  ee$ImageCollection$filterDate("2018-01-01", "2018-12-31") %>%
  ee$ImageCollection$select(c("B5"))%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()


land_nir = image$reduce(ee$Reducer$stdDev())
land_nir = land_nir$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

land_nirr <- ee_as_raster(
  image = land_nir,
  scale= res,
  region = region,
  via = "drive"
)

writeRaster(land_nirr, '01-Data/covs/land_sd_nir.tif', overwrite=T)

# OpenLandMap soil water content 0-10-30 ----

image <- ee$Image("OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01") %>%
  ee$Image$select(c('b0','b10','b30'))%>%
  ee$Image$clip(region)

soil_wt = image$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

soil_wtr <- ee_as_raster(
  image = soil_wt,
  scale= res,
  region = region,
  via = "drive"
)

#Harmonize to 0-30 depth with a weighted average
WeightedAverage<-function(r){return(r[[1]]*(1/30)+r[[2]]*(9/30)+r[[3]]*(20/30))}

soil_wtr<-overlay(soil_wtr,fun=WeightedAverage)

writeRaster(soil_wtr, '01-Data/covs/soil_wtr.tif', overwrite=T)

# OpenLandMap Potential FAPAR Monthly ----

image <- ee$Image("OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01") %>%
  ee$Image$select(c('jan','feb','mar','apr','may','jun','jul',
                   'aug','sep','oct','nov','dec' ))%>%
  ee$Image$clip(region)

fapar = image$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

fapar=fapar$reduce(ee$Reducer$mean())

faparr <- ee_as_raster(
  image = fapar,
  scale= res,
  region = region,
  via = "drive"
)

writeRaster(faparr, '01-Data/covs/faparr_mean.tif', overwrite=T)



# Multi-Scale Topographic Position Index ----
image <- ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI") %>%
  ee$Image$select('AVE')%>%
  ee$Image$clip(region)

top_pos = image$resample('bilinear')$reproject(
  crs= crs,
  scale= res)

top_posr <- ee_as_raster(
  image = top_pos,
  scale= res,
  region = region,
  via = "drive"
)

writeRaster(top_posr, '01-Data/covs/top_posr.tif', overwrite=T)
