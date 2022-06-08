#######################################################
#
#  Select covariates
#  based on PCA and Correlation
#
# GSP-Secretariat
# Contact: Isabel.Luotto@fao.org
#
#######################################################

#Empty environment and cache 
rm(list = ls());
gc()

#######################################################
#
#  User defined variables: ----

# Working directory
#wd <- 'C:/Users/luottoi/Documents/GitHub/Digital-Soil-Mapping'
wd <- 'C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping'
#List of soil attributes prepared in script #2
soilatt <- c('OCS','clay', 'pH')
#
#
#######################################################
# Set working directory
setwd(wd)
# load packages

library(sf)
library(terra)
library(data.table)
library(reshape)

# Selection based on correlation ----


# list all .tif files the 'covs' folder and load them in a raster stack.
# all rasters should have same extent, resoultion and coordinate system
files <- list.files(path = "01-Data/covs", pattern = "tif*$", full.names = TRUE)


covs <- rast(files)


# Explore the data
names(covs)

for (i in unique(soilatt)){

# Load the processed data for digital soil mapping. 
# This table was prepared in the 'data_preparation_profiles' script
dat <- fread(paste0("02-Outputs/",i,"_dat.csv"))

#Use coordinates to turn point data into a vector with coordinates
# WGS84
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

dat <- st_as_sf(dat,                         
               coords = c("X", "Y"),
               crs = projcrs)

# Extract values from covariates to the soil points
pv <- extract(x = covs, y = vect(dat),xy=F)
dat <- cbind(dat,pv)

summary(dat)

# Remove NA values
dat <-as.data.frame(dat)
dat[, 'geometry'] <-NULL

dat <- dat[complete.cases(dat),]

# Test correlation between each covariate and the soil attribute
names(dat)


test_covs <- cor(x = as.matrix(dat[,i]),
                 y = as.matrix(dat[,names(covs)]))
test_covs 

# Select only the covariates that have correlation higher than 0.3

x <- subset(melt(test_covs), value != 1 | complete.cases(test_covs))
x <- x[with(x, order(-abs(x$value))),]
(x <- subset(x,abs(x$value)>0.20))
selection <- as.character(x$X2)

# Leave only selected covariates in the 'covs' raster stack
covs <- covs[[selection]]
names(covs)


writeRaster(covs, file = paste0("02-Outputs/", i,"_covariates.tif"),overwrite=TRUE)
print(i)
}


#Add categorical layers

# # Import land cover layer
# landcover <- raster('01-Data/land cover/LandCover.tif')
# landcover <- projectRaster(landcover, covs$EVI,)
# 
# # Import and explore the soil map (vector polygon data)
# soilmap<-shapefile("01-Data/Soil map/SoilTypes.shp")
# plot(soilmap)
# summary(soilmap)
# 
# # We will use 'Symbol' attribute as an indication of soil type. It has to be 'factor' type.
# soilmap@data$Symbol <- as.factor(soilmap@data$Symbol)
# 
# # Rasterize the soil map using one of the 'covs' as a template and 'Symbol' as value field
# soilmap.r <- rasterize(x = soilmap, fun='first', y = covs$Terrain_attributes.1, field = "Symbol")
# 
# # Save the rasterized map
# writeRaster(soilmap.r, '02-Outputs/Soil_map.tif', overwrite=TRUE)
# 
# 
# for (i in unique(soilatt)){
#   
#  load(file = paste0("02-Outputs/", i,"_covariates.RData"))
# 
#   # stack it with the other rasters
#   covs <- stack(covs, landcover)
#   
#   # Now we can stack it with the other rasters
#   covs <- stack(covs, soilmap.r)
#   names(covs)
#   # correct the name 
#   names(covs)[names(covs)=='layer'] <- "soilmap"
#   
#   # Mask the covariates with the country mask
#   mask <- shapefile('01-Data/MKD.shp')
#   covs <- mask(x = covs, mask = mask)
#   
#   save(covs, file = paste0("02-Outputs/", i,"_covariates.RData"))
# 
# } 





