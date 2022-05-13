#######################################################
#
#  Process and download topographic 
#  OpenLandMap
#  from zenodo to R 
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

# Folder to store global layers
output_dir <-'C:/Users/hp/Documents/FAO/data/OpenLandMap/'
# Area of interest
AOI <- '01-Data/MKD.shp'

#Resolution 2km 1km, 250 m or 500 m
res <- '1km'

#
#######################################################
#set working directory 
setwd(wd)

#install.packages('zen4R')
library(zen4R)
library(terra)

# Get links for the OpenLandMap topographic attributes
# to download from Zenodo 

#instantiate Zen4R client
zenodo <- ZenodoManager$new()

# Get file record
my_rec <- zenodo$getRecordByDOI("10.5281/zenodo.1447210")

sel.tif = my_rec$files

# Extract list of links
links <- data.frame()
for (i in 1:length(sel.tif)){

  x<-as.data.frame(sel.tif[[i]][["links"]][["download"]])
    links <- rbind(links,x)


}

names(links) <-'Links'

# Selection of topographic attributes
seltop <- c('slope', 'twi', 'vbf','_curvature','downlslope.curvature','dvm2','dvm','mrn','tpi')
# Download global layers (to be done once) ---

for (i in unique(seltop)){
link <- links[grep(i, links$Links),]
link <- link[grep(res, link)]

download.file(link, paste0(output_dir, 'olm_',i,'_',res,'.tif'),mode = "wb")
}

# Clip and store covariates in working directory
AOI <- vect(AOI)

files <- list.files(path = output_dir, pattern = res, full.names = T)
  covs <- rast(files)

  covs <- crop(covs, AOI)
  covs <- mask(covs, AOI)
 
  
  writeRaster(covs, '01-Data/covs/olm_covs.tif', overwrite=T)  

