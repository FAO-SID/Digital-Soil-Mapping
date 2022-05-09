#########################################
#
# GSOCmap v2 & GSOCseq V2 | Module I
# SOC map and Clay layer preparation 
# Covariates
#
# GSP-Secretariat
# Contact: Isabel.Luotto@fao.org
#
#########################################
#Empty environment and cache ----
rm(list = ls());
gc()

# Load libraries ----
library(data.table)
library(terra)
library(sf)
library(rgee)
# Set working directory Module I ----
wd <- 'C:/Users/luottoi/Documents/GitHub/Digital-Soil-Mapping'

setwd(wd)

# Country shapefile
AOI <- read_sf('01-Data/MKD.shp')


#List of covariates to prepare
# Mean annual temperature
# Total annual Precipitation
# Precipitation of wettest month
# Precipitation of driest month

# TAGEE 13 soil attributes (list below)

# MODIS EVI & NDVI

# Daytime temperature SD

# Landsat 8 RED and NIR standard deviation

# Mean annual temperature (daytime) ----
ee_Initialize()

region <- sf_as_ee(AOI)
region = region$geometry()


image1 <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate("2017-01-01", "2017-12-31") %>%
  ee$ImageCollection$select("tmmx")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

  # from imagecollection to image
image2 <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate("2017-01-01", "2017-12-31") %>%
  ee$ImageCollection$select("tmmx")%>%
  ee$ImageCollection$toBands()


image1 <- image1$multiply(0.1)
image2 <- image2$multiply(0.1)

diff <- image1$add(image2)
avT = diff$divide(2)
avT = avT$reduce(ee$Reducer$mean())


avtr <- ee_as_raster(
  image = avT,
  scale= 4000,
  region = region,
  via = "drive"
)

writeRaster(avtr, '01-Data/covs/avtr.tif', overwrite=T)

# Total annual Precipitation ----
ee_Initialize()



image <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate("2017-01-01", "2017-12-31") %>%
  ee$ImageCollection$select("pr")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$sum()%>%
  ee$Image$toDouble()




Prr <- ee_as_raster(
  image = image,
  scale= 4000,
  region = region,
  via = "drive"
)

writeRaster(Prr, '01-Data/covs/Prr.tif', overwrite=T)

# Precipitation of wettest month ----
image <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate("2017-01-01", "2017-12-31") %>%
  ee$ImageCollection$select("pr")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()


Prr_all <- ee_as_raster(
  image = image,
  scale= 4000,
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


# Terrain attributes using TAGEE 13 covariates
# Attribute	| Unit	| Description
# Elevation  |	meter	|Height of terrain above sea level
# Slope	degree	| Slope |  gradient
# Aspect |	degree	Compass direction
# Hillshade	| dimensionless	Brightness of the illuminated terrain
# Northness |	dimensionless	Degree of orientation to North
# Eastness	| dimensionless	Degree of orientation to East
# Horizontal curvature |	meter	Curvature tangent to the contour line
# Vertical curvature |	meter	Curvature tangent to the slope line
# Mean curvature |	meter	Half-sum of the two orthogonal curvatures
# Minimal curvature	meter |	Lowest value of curvature
# Maximal curvature	meter |	Highest value of curvature
# Gaussian curvature	 | meter	Product of maximal and minimal curvatures
# Shape Index	| dimensionless	Continuous form of the Gaussian landform classification


TAGEE <- import("tagee")

image <- ee$Image("MERIT/DEM/v1_0_3") %>%
  ee$Image$clip(region)%>%
  ee$Image$toDouble()


DEMAttributes = TAGEE$terrainAnalysis( image, region)

tagee_test <- ee_as_raster(
  image = DEMAttributes,
  scale= 1000,
  region = region,
  via = "drive"
)

writeRaster(tagee_test, '01-Data/covs/Terrain_attributes.tif', overwrite=T)

# EVI & NDVI ----
EVI <- ee$ImageCollection("MODIS/061/MOD13Q1") %>%
  ee$ImageCollection$filterDate("2017-01-01", "2017-12-31") %>%
  ee$ImageCollection$select("EVI")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

NDVI <- ee$ImageCollection("MODIS/061/MOD13Q1") %>%
  ee$ImageCollection$filterDate("2017-01-01", "2017-12-31") %>%
  ee$ImageCollection$select("NDVI")%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()

ndvir <- ee_as_raster(
  image = NDVI,
  scale= 1000,
  region = region,
  via = "drive"
)

evir <- ee_as_raster(
  image = EVI,
  scale= 1000,
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


sd_d_Tr <- ee_as_raster(
  image = sd_d_T,
  scale= 1000,
  region = region,
  via = "drive"
)


writeRaster(sd_d_Tr, '01-Data/covs/sd_d_Tr.tif', overwrite=T)

# Landsat bands mean and sd 
image <- ee$ImageCollection("LANDSAT/LC08/C02/T1_RT") %>%
  ee$ImageCollection$filterDate("2018-01-01", "2018-12-31") %>%
  ee$ImageCollection$select(c("B4"))%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()


land_red = image$reduce(ee$Reducer$stdDev())

land_sd_red <- ee_as_raster(
  image = land_red,
  scale= 1000,
  region = region,
  via = "drive"
)

image <- ee$ImageCollection("LANDSAT/LC08/C02/T1_RT") %>%
  ee$ImageCollection$filterDate("2018-01-01", "2018-12-31") %>%
  ee$ImageCollection$select(c("B5"))%>%
  ee$ImageCollection$filterBounds(region)%>%
  ee$ImageCollection$toBands()


land_nir = image$reduce(ee$Reducer$stdDev())

land_nirr <- ee_as_raster(
  image = land_nir,
  scale= 1000,
  region = region,
  via = "drive"
)

writeRaster(land_nirr, '01-Data/covs/land_sd_nir.tif', overwrite=T)

