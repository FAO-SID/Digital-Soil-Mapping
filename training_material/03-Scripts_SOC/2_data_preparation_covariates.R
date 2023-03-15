# This script is a part of the FAO/GSP training course on Digital Soil Mapping.
# It is used to prepare and select covariates for digital soil mapping of SOC stocks.
# It can be modified to be used for mapping other soil properties.
#################################################################################################################

# Remove existing objects in memory and clear memory cache
rm(list = ls());
gc()

# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

#####################################  Import and stack the covariates ########################################
library (sp)
library(raster)

# list all .tif files the 'covs' folder and load them in a raster stack.
# all rasters should have same extent, resoultion and coordinate system
files <- list.files(path = "01-Data/covs", pattern = "tif*$", full.names = TRUE)
covs <- stack(files)

# Explore the data
names(covs)
plot(covs$DEMENV5)

#####################################  Select the covariates using correlation analysis #############################

# Load the processed data for digital soil mapping. 
# This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv("02-Outputs/dat_train.csv")

# Upgrade points data frame to SpatialPointsDataFrame and define their coordinate system
coordinates(dat) <- ~ X + Y
proj4string(dat) = CRS("+init=epsg:4326") # WGS84

# Check that the points overlay with the rasters
plot(covs$DEMENV5)
points(dat)

# Extract values from covariates to the soil points
dat <- extract(x = covs, y = dat, sp = TRUE)
summary(dat)

# Remove NA values
dat<-as.data.frame(dat)
dat <- dat[complete.cases(dat),]

# Test correlation between each covariate and the 'OCSKGMlog'
names(dat)
test_covs <- cor(x = as.matrix(dat[,8]),
                     y = as.matrix(dat[,c(9:33)]))
test_covs 

# Select only the covariates that have correlation higher than 0.3
library(reshape)

x <- subset(melt(test_covs), value != 1 | value != NA)
x <- x[with(x, order(-abs(x$value))),]
(x <- subset(x,abs(x$value)>0.3))
selection <- as.character(x$X2)

# Leave only selected covariates in the 'covs' raster stack
covs <- covs[[selection]]
names(covs)

############################################### Add categorical covariates #########################################

# Import land cover layer
landcover <- raster ('01-Data/land cover/LandCover.tif')
plot(landcover)

# Try to stack LandCover raster with the rest
covs <- stack(covs, landcover)

# Error indicates different extent. Check coordinate system (crs) of covs and LandCover
(covs@crs) 
(landcover@crs)  

# 'covs' has crs WGS 84, while landcover has crs UTM Zone 34
# We need to reproject the 'Lancover' raster using one of the 'covs' as a template 
landcover <- projectRaster(from = landcover, to = covs$DEMENV5, method = "ngb")

# Save the reprojected raster
writeRaster(landcover, '02-Outputs/Landcover_WGS84.tif', overwrite=TRUE)

# Now we can stack it with the other rasters
covs <- stack(covs, landcover)

# Import and explore the soil map (vector polygon data)
soilmap<-shapefile("01-Data/Soil map/SoilTypes.shp")
plot(soilmap)
summary(soilmap)

# We will use 'Symbol' attribute as an indication of soil type. It has to be 'factor' type.
soilmap@data$Symbol <- as.factor(soilmap@data$Symbol)

# Rasterize the soil map using one of the 'covs' as a template and 'Symbol' as value field
soilmap.r <- rasterize(x = soilmap, fun='first', y = covs$DEMENV5, field = "Symbol")

# Explore the result
plot(soilmap.r, col=rainbow(20))
legend("bottomright",legend = levels(soilmap$Symbol), fill=rainbow(20),
       cex=0.5)

# Save the rasterized map
writeRaster(soilmap.r, '02-Outputs/Soil_map.tif', overwrite=TRUE)

# Now we can stack it with the other rasters
covs <- stack(covs, soilmap.r)
names(covs)
# correct the name 
names(covs)[names(covs)=='layer'] <- "soilmap"

########################################### Prepare the final stack of covariates for modelling #######################

# Mask the covariates with the country mask
plot(covs$DEMENV5)
mask <- raster("01-Data/mask.tif")
plot(mask)
covs <- mask(x = covs, mask = mask)
plot(covs$DEMENV5)

# export all the covariates in the R-object
save(covs, file = "02-Outputs/covariates.RData")
