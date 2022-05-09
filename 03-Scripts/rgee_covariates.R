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
setwd("C:/Users/hp/Documents/FAO/GSOCseq/GSOCseq v2/Module I")

AOI <- read_sf('01-Data/MKD.shp')


#List of covariates to prepare
# Mean annual temperature
# Total annual Precipitation
# Precipitation of wettest month
# Precipitation of driest month

# Terrain slope
# Aspect
# convergence index
# distance to channel network
# channel network base level
# relative slope position
# LS-factor
# Multiresolution index of valley bottom flatness (MRVBF)
# Valley depth
# Wetness index

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

writeRaster(avtr, '01-Data/covs/avtr.tif')

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

writeRaster(Prr, '01-Data/covs/Prr.tif')

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

writeRaster(Prr_wet, '01-Data/covs/Prr_wet.tif')

# Precipitation of driest month ----

dryest <- names(sums[sums ==min(sums)])
Prr_dry <- Prr_all[[dryest]]

writeRaster(Prr_wet, '01-Data/covs/Prr_dry.tif')

# Terrain slope ----
image <- ee$Image("MERIT/DEM/v1_0_3") %>%
ee$Image$clip(region)

slope = ee$Terrain$slope(image)
  
ter_slop <- ee_as_raster(
  image = slope,
  scale= 1000,
  region = region,
  via = "drive"
)

writeRaster(ter_slop, '01-Data/covs/ter_slop.tif')

# Profile aspect ----
asp = ee$Terrain$aspect(image)

aspect <- ee_as_raster(
  image = asp,
  scale= 1000,
  region = region,
  via = "drive"
)


writeRaster(ter_slop, '01-Data/covs/aspect.tif')

