#######################################################
#
# Quantile Regression Forest
# Soil Property Mapping
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

AOI <- '01-Data/MKD.shp'
soilatt <- 'OCS'


#
#
#######################################################
# Set working directory
setwd(wd)
#load packages
library(tidyverse)
library(data.table)
library(caret)
library(quantregForest)
library(terra)
library(sf)
library(doParallel)


# Make a selection of covariates, calibrate a QRF model, make a prediction ----

#Load data and covariates 
files <- list.files(path= '01-Data/covs/', pattern = '.tif*', full.names = T)
files <- files[!(grepl('aux', files))]
covs <-rast(files)

# Load the processed data for digital soil mapping (Script 2) as a spatial object
dat <- read.csv(paste0("02-Outputs/",soilatt,"_dat.csv"))
dat <- vect(dat, geom=c("X", "Y"))
names(dat)

# Extract values from covariates to the soil points
pv <- terra::extract(x = covs, y = dat,xy=F)
dat <- cbind(dat,pv)

summary(dat)

# Remove NA values
dat <-as.data.frame(dat)
dat <- dat[complete.cases(dat),]

# Define a formila  and QRF parameters (set 3 times repeated 10-fold cross-validation)

fitControl <- rfeControl(functions = rfFuncs,
                         method = "repeatedcv",
                         number = 10,         ## 10 -fold CV
                         repeats = 3,        ## repeated 3 times
                         verbose = TRUE,
                         saveDetails = TRUE)


tuneGrid <-  expand.grid(mtry = c(100, 200, 500))

#Calibrate the model using multiple cores
fm = as.formula(paste(soilatt," ~", paste0(names(covs),
                                             collapse = "+")))


#Calibrate the model using multiple cores
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)


# Select covariates based on recursive feature elimination ----
covsel <- rfe(fm,
             data = dat,  
             sizes = 10:(length(dat)-2),
             #method = "rf",
             rfeControl = fitControl,
             verbose = TRUE,
             tuneGrid = tuneGrid,
             keep.inbag = T)
#Extract selection of covariates and subset covs
covsel <- covsel$optVariables 
covs <- covs[[covsel]]

#update formula
#Calibrate the model using multiple cores
fm = as.formula(paste(soilatt," ~", paste0(names(covs),
                                           collapse = "+")))


# Run quantile regression forest
model <- caret::train(fm,
                      data = dat[, c(soilatt,names(covs))],
                      method = "qrf",
                      trControl = fitControl,
                      verbose = TRUE,
                      tuneGrid = tuneGrid,
                      keep.inbag = T)
print(model)


#Extract predictor importance as relative values (%)
model$importance <-
  data.frame(var=rownames(importance(model$finalModel)),
             importance(model$finalModel)/sum(importance(model$finalModel))*100) %>%
  arrange(desc(IncNodePurity))

# print(model)
saveRDS(model, file = paste0("02-Outputs/models/model_",soilatt,".rds"))

stopCluster(cl)
gc()



# Generation of maps (prediction of soil attributes) ----
#Produce tiles
r <-covs[[1]]
t <- rast(nrows = 5, ncols = 5, extent = ext(r), crs = crs(r))
tile <- makeTiles(r, t,overwrite=TRUE,filename="02-Outputs/tiles/tiles.tif")



#Predict soil attributes per tiles 
# loop to predict on each tile


for (j in seq_along(tile)) {
  gc()
  t <- rast(tile[j])
  covst <- crop(covs, t)
  
  
  # plot(r)# 
  pred_mean <- terra::predict(covst, model = final_mod, na.rm=TRUE,  
                              cpkgs="randomForest", what=mean,
                              # filename = paste0("02-Outputs/Final Maps/tiles/OCS_tile_", j, ".tif"),
                              # overwrite = TRUE
  )
  pred_sd <- terra::predict(covst, model = final_mod, na.rm=TRUE,  
                            cpkgs="randomForest", what=sd)  
  
  
  
  # ###### Raster package solution (in case terra results in many NA pixels)
  # library(raster)
  # covst <- stack(covst)
  # class(final_mod$finalModel) <-"quantregForest"
  # # Estimate model uncertainty
  # pred_sd <- predict(covst,model=final_mod$finalModel,type=sd)
  # # OCSKGMlog prediction based in all available data
  # pred_mean <- predict(covst,model=final_mod)
  # 
  # 
  # ##################################  
  
  
  
  writeRaster(pred_mean, filename = paste0("02-Outputs/tiles/soilatt_tiles/",soilatt,"_tile_", j, ".tif"), 
              overwrite = TRUE)
  writeRaster(pred_sd, filename = paste0("02-Outputs/tiles/soilatt_tiles/",soilatt,"_tileSD_", j, ".tif"), 
              overwrite = TRUE)
  
  rm(pred_mean)
  rm(pred_sd)
  
  
  print(paste("tile",tile[j]))
}

#Merge tiles both prediction and st.Dev
f_mean <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", pattern = paste0(soilatt,"_tile_"), full.names = TRUE)
f_sd <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", pattern =  paste0(soilatt,"_tileSD_"), full.names = TRUE)
r_mean_l <- list()
r_sd_l <- list()
for (g in 1:length(f_mean)){
  
  r <- rast(f_mean[g])
  r_mean_l[g] <-r
  rm(r)
}

for (g in 1:length(f_sd)){
  
  r <- rast(f_sd[g])
  r_sd_l[g] <-r
  rm(r)
}

r_mean <-sprc(r_mean_l)
r_sd <-sprc(r_sd_l)
pred_mean <- mosaic(r_mean)
pred_sd <- mosaic(r_sd)

AOI <- vect(AOI)
pred_mean <- mask(pred_mean,AOI)
pred_sd <- mask(pred_sd,AOI)


plot(pred_mean)
plot(pred_sd)


# Export final rasters ---- 
writeRaster(pred_mean, paste0("02-Outputs/maps/",soilatt,"_QRF.tif"),overwrite=TRUE)
writeRaster(pred_sd, paste0("02-Outputs/maps/",soilatt,"_QRF_SD.tif"),overwrite=TRUE)






