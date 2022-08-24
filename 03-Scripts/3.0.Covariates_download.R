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

# Loop over the names of assets to clip and dowload the covariates
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

###########################################################################
# Export PCs based on Principle Component Analysis ----
assetname_clim <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/CHELSA")
assetname_lc <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/LANDCOVER")
assetname_ta <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/OPENLANDMAP")
assetname_veg <-  ee_manage_assetlist(path_asset = "projects/digital-soil-mapping-gsp-fao/assets/MODIS")


assetname <- rbind(assetname_clim, assetname_lc, assetname_ta,assetname_veg)

# clip images
for (i in unique(assetname$ID)){
  #Extract filename 
  filename <- sub('.*\\/', '', i)
  
  #Clip image to the extent of the AOI
  image <- ee$Image(i) %>%
    ee$Image$clip(region)
  
  assign(filename, image)
}


# Stack bands 
assetname$names <- sub('.*\\/', '', assetname$ID)

SG <- bio1

assetname <- assetname[assetname$names!= 'bio1',]

for (i in unique(assetname$names)){
  
  image <- get(i)
  SG <-SG$addBands(image)
}

inputBandNames <- SG$bandNames()$getInfo()
print(inputBandNames)

dimOne <- length(inputBandNames)
cat("Number of input bands:", dimOne)

# Calculate scale standardize each band to mean=0, s.d.=1. ----
scale <- SG$select("b1")$projection()$nominalScale()
cat("Nominal scale: ", scale$getInfo())

SGmean <- SG$reduceRegion(ee$Reducer$mean(),
                          geometry = region, scale = scale, bestEffort = TRUE)
head(SGmean$getInfo(),3) # an example of the means

SGsd <- SG$reduceRegion(ee$Reducer$stdDev(),
                        geometry = region, scale = scale, bestEffort = TRUE)
head(SGsd$getInfo(),3)

# saveRDS(as.vector(SGmean$getInfo()), file = "./inputBandMeans.rds")
# saveRDS(as.vector(SGsd$getInfo()), file = "./inputBandSDs.rds")

SGmean.img <- ee$Image$constant(SGmean$values(inputBandNames))
SGsd.img <- ee$Image$constant(SGsd$values(inputBandNames))

#Standardize
SGstd <- SG$subtract(SGmean.img)
SGstd <- SGstd$divide(SGsd.img)


SGstd.minMax <- SGstd$reduceRegion(ee$Reducer$minMax(),
                                   geometry = region, scale = scale,
                                   maxPixels = 1e9, bestEffort = TRUE)
minMaxNames <- names(SGstd.minMax$getInfo())
minMaxVals <- SGstd.minMax$values()$getInfo()
# per property min/max
head(data.frame(property_depth_q=minMaxNames, value=minMaxVals), 12)

# overall min/max
print(c(min(minMaxVals), max(minMaxVals)) )

# convert image to an array of pixels, for matrix calculations
# dimensions are N x P; i.e., pixels x bands
arrays <- SGstd$toArray()


# Compute the covariance of the bands within the region.
covar <- arrays$reduceRegion(
  ee$Reducer$centeredCovariance(),
  geometry = region, scale = scale,
  maxPixels = 1e6,
  bestEffort = TRUE
)


# Get the covariance result and cast to an array.
# This represents the band-to-band covariance within the region.
covarArray <- ee$Array(covar$get('array'))
# note we know the dimensions from the inputs
# Perform an eigen analysis and slice apart the values and vectors.
eigens <- covarArray$eigen()
# the first item of each slice (PC) is the corresponding eigenvalue
# the remaining items are the eigenvalues (rotations) for that PC

# by removing the first axis (eigenvalues) and converting to a vector
# array$slice(axis=0, start=0, end=null, step=1)
# here we only have one axis, indexed by the PC
eigenValues <- eigens$slice(1L, 0L, 1L)
cat('Eigenvalues:')
## Eigenvalues:
print(eigenValues.vect <- unlist(eigenValues$getInfo()))

# compute and show proportional eigenvalues
eigenValuesProportion <-
  round(100*(eigenValues.vect/sum(eigenValues.vect)),2)
cat('PCs percent of variance explained:')
## PCs percent of variance explained:
print(eigenValuesProportion)

plot(eigenValuesProportion, type="h", xlab="PC",
     ylab="% of variance explained",
     main="Standardized PCA")

evsum <- cumsum(eigenValuesProportion)
cat('PCs percent of variance explained, cumulative sum:')
## PCs percent of variance explained, cumulative sum:
print(round(evsum,1))

cat(npc95 <- which(evsum > 95)[1],
    "PCs are needed to explain 95% of the variance")

#saveRDS(eigenValues.vect, file = "./eigenValuesVector.rds")


# The eigenvectors (rotations); this is a PxP matrix with eigenvectors in rows.
eigenVectors <- eigens$slice(1L, 1L)
# show the first few rotations
cat('Eigenvectors 1--3:')
## Eigenvectors 1--3:
print(data.frame(band=inputBandNames,
                 rotation1 = eigenVectors$getInfo()[[1]],
                 rotation2 = eigenVectors$getInfo()[[2]],
                 rotation3 = eigenVectors$getInfo()[[3]]
))

#export the eigenvectors ----
eVm <- matrix(unlist(eigenVectors$getInfo()), byrow = TRUE, nrow = dimOne)
#saveRDS(eVm, file = "./eigenvectorMatrix.rds")

arrayImage <- arrays$toArray(1L)
PCsMatrix <- ee$Image(eigenVectors)$matrixMultiply(arrayImage)

# arrayProject: "Projects the array in each pixel to a lower dimensional #
#space by specifying the axes to retain"
# arrayFlatten: "Converts a single band image of equal-shape
# multidimensional pixels
# to an image of scalar pixels,
# with one band for each element of the array."
PCbandNames <- as.list(paste0('PC', seq(1L:dimOne)))
PCs <- PCsMatrix$arrayProject(list(0L))$arrayFlatten(
  list(PCbandNames)
)

PCs95 <- PCs$select(0L:(npc95-1L)) # indexing starts from 0
PCs95$bandNames()$getInfo()

# Export PCs as raster stack ----
PCAs_covs <- ee_as_raster(
  image = PCs95,
  scale= res,
  region = region,
  via = "drive"
)



writeRaster(PCAs_covs, paste0(output_dir, 'PCAs_covs.tif'), overwrite=T)

