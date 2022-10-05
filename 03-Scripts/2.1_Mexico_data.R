#Empty environment and cache 
rm(list = ls())
gc()

# 0 - User-defined variables ===================================================

setwd('C:/GIT/Digital-Soil-Mapping') # change the path accordingly

library(tidyverse) # for data management and reshaping
library(readxl) # for importing excel files
library(mapview) # for seeing the profiles in a map
library(sf) # to manage spatial data (shp vectors) 
library(terra)


dxy <- read_csv("01-Data/Mexico/L0_FAO_SAGARPA_DATASET - suelos puntos de muestreo.csv")
names(dxy)
dxy <- dxy %>% dplyr::select(LabID=ID_SAMPLE, x=Y_COORD, y=X_COORD, STATE, MUNICIPALITY, p_bray=P, k=K, ca=Ca)
dxy <- dxy %>% na.omit()

AOI <- st_read("01-Data/Mexico/AOI.shp")
plot(AOI)
eco <- st_read("01-Data/Mexico/ecort08cw.shp")
croplands <- st_read("01-Data/Mexico/Frontera_Agricola_Sin_pastos.shp")

eco <- st_transform(eco,crs = st_crs(AOI))
croplands <- st_transform(croplands,crs = st_crs(AOI))

eco <- st_crop(eco, AOI)

r <- rast("01-Data/covs/bio12.tif")
croplands <- vect(croplands)
crop <- rasterize(croplands, r)
plot(crop)
s <- vect(dxy, geom = c("x","y"), crs = 4326)

eco <- vect(eco)

d <- terra::intersect(s,eco)
d <- d[,1:6]
eco <- eco[is.related(eco, d, "intersects") |> which()]
eco$DESECON4 <- as.factor(eco$DESECON4)
eco$stratum <- as.numeric(eco$DESECON4)
# plot(eco)

ecor <- rasterize(eco, crop, field = "stratum")
plot(ecor)
ecor_crop <- terra::mask(ecor, crop)
plot(ecor_crop)
plot(d, add=TRUE)


d <- cbind(d,terra::extract(ecor_crop, d, ID=FALSE, xy=TRUE))

d <- d %>% as_tibble() %>% na.omit()

write_csv(d, "01-Data/Mexico/p_bray_data_coord.csv")
writeRaster(ecor_crop, "01-Data/Mexico/strata.tif")
