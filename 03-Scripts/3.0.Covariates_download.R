#######################################################
#
#  Process and download covariates
#  from GEE R 
#
#  Export both raw covariates and PCAs
#
#  GSP-Secretariat
#  Contacts:  Marcos.Angelini@fao.org
#             Isabel.Luotto@fao.org
#
#######################################################

#Empty environment and cache ----
rm(list = ls());
gc()

#######################################################
#
#  User defined variables:

# Working directory
#wd <- 'C:/Users/luottoi/Documents/GitHub/Digital-Soil-Mapping'
wd <- 'C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping'

# Output covariate folder
#output_dir <-''
output_dir <-'01-Data/covs/'

# Area of interest: either own shapefile or 3-digit ISO code to extract from UN 2020 boundaries
#AOI <- '01-Data/MKD.shp'
AOI <- 'MKD'
#Start and End time 
start_T <- "2017-01-01"
end_T <- "2017-12-31"

# Resolution and projectio
res = 250
crs = "EPSG:4326"

#
#
#######################################################


# Load libraries ----
library(raster)
library(sf)
library(rgee)

# Set working directory ----
setwd(wd)

# Upload own AOI shapefile ----
#AOI <- read_sf(AOI)
# convert AOI to a box polygon
#AOI <- st_as_sfc(st_bbox(AOI))
#AOI <- st_as_sf(AOI)


#Overview of covariates ----
# CLIMATIC VARIABLES from CHELSA
# VEGETATION INDICES, FPAR and LAND SURFACE TEMPERATURE from MODIS
# LAND COVER LAYERS from Dynamic World 10m near-real-time (NRT) 
# TERRAINE attributes from OpenLandMap

# for more information about the single covariates: open covariates.xslx in the 
# training material folder

#Initialize GEE ----
ee_Initialize()

#Initial setup either convert shp to gee geometry or extract from UN 2020 map----

#Convert shp to gee geometry
#region <- sf_as_ee(AOI)
#region = region$geometry()

#Extract from UN 2020 map using ISO code ----
region <-ee$FeatureCollection("projects/digital-soil-mapping-gsp-fao/assets/UN_BORDERS/BNDA_CTY")%>%
  ee$FeatureCollection$filterMetadata('ISO3CD', 'equals', AOI)

AOI_shp <-ee_as_sf(region)
write_sf(AOI_shp, paste0('01-Data/',AOI,'.shp'))

region = region$geometry()

# CLIMATIC VARIABLES ----
#Obtain list of climatic variables
assetname_clim <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/CHELSA")

# Loop over the names of assets to clip and dowload the covariates
for (i in unique(assetname_clim$ID)){
  
  #Extract filename 
  filename <- sub('.*\\/', '', i)
  
  #Clip image to the extent of the AOI
  image <- ee$Image(i) %>%
    ee$Image$clip(region)
  
  # Resample to target resolution
    image = image$resample('bilinear')$reproject(
    crs= crs,
    scale= res)
  
  #Export clipped covariate as raster
  raster <- ee_as_raster(
    image = image,
    scale= res,
    region = region,
    via = "drive"
  )
  
  writeRaster(raster, paste0(output_dir,filename, '.tif'), overwrite=T)
  print(paste(filename, 'exported successfully'))
  
}

# VEGETATION INDICES ----
#Obtain list of modis variables
assetname_veg <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/MODIS")

# Loop over the names of assets to clip and dowload the covariates
for (i in unique(assetname_veg$ID)){
  
  #Extract filename 
  filename <- sub('.*\\/', '', i)
  
  #Clip image to the extent of the AOI
  image <- ee$Image(i) %>%
    ee$Image$clip(region)
  
  # Resample to target resolution
    image = image$resample('bilinear')$reproject(
    crs= crs,
    scale= res)
  
  #Export clipped covariate as raster
  raster <- ee_as_raster(
    image = image,
    scale= res,
    region = region,
    via = "drive"
  )
  
  writeRaster(raster, paste0(output_dir,filename, '.tif'), overwrite=T)
  print(paste(filename, 'exported successfully'))
  
}


#  LAND COVER LAYERS ----
#Obtain list of landcover variables
assetname_lc <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/LANDCOVER")

# Loop over the names of assets to clip and dowload the covariates
for (i in unique(assetname_lc$ID)){
  
  #Extract filename 
  filename <- sub('.*\\/', '', i)
  
  #Clip image to the extent of the AOI
  image <- ee$Image(i) %>%
    ee$Image$clip(region)
  
  # # Resample to target resolution
  #   image = image$resample('bilinear')$reproject(
  #   crs= crs,
  #   scale= res)
  
  #Export clipped covariate as raster
  raster <- ee_as_raster(
    image = image,
    scale= res,
    region = region,
    via = "drive"
  )
  
  writeRaster(raster, paste0(output_dir,filename, '.tif'), overwrite=T)
  print(paste(filename, 'exported successfully'))
  
}

# TERRAIN ATTRIBUTES ----
#Obtain list of terrain attributes
assetname_ta <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/OPENLANDMAP")

# Loop over the names of assets to clip and download the covariates
for (i in unique(assetname_ta$ID)){
  
  #Extract filename 
  filename <- sub('.*\\/', '', i)
  
  #Clip image to the extent of the AOI
  image <- ee$Image(i) %>%
    ee$Image$clip(region)
  
  # # Resample to target resolution
  #   image = image$resample('bilinear')$reproject(
  #   crs= crs,
  #   scale= res)
  
  #Export clipped covariate as raster
  raster <- ee_as_raster(
    image = image,
    scale= res,
    region = region,
    via = "drive"
  )
  
  writeRaster(raster, paste0(output_dir,filename, '.tif'), overwrite=T)
  print(paste(filename, 'exported successfully'))
  
}



