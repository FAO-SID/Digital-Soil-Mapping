#######################################################
#
# Quantile Regression Forest
# Bulk soil attributes
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

AOI <- '01-Data/MKD.shp'
#
#
#######################################################
# Set working directory
setwd(wd)


#Load packages
library(tidyverse)
library(data.table)
library(caret)
library(quantregForest)
library(terra)
library(sf)
library(doParallel)




# Calibrate QRF models, perform QA and select the best model ----
for(i in unique(soilatt)){
gc()
#Load data and covariates stacks based on correlation
covs <-rast( paste0("02-Outputs/", i,"_covariates.tif"))

# Load PC covariates (the same for all soil attributes)
covs_pc <- rast('02-Outputs/PCA_covariates.tif')
covs_all <-  c(covs, covs_pc)  
# Load the processed data for digital soil mapping (Script 2)
dat <- read.csv(paste0("02-Outputs/",i,"_dat.csv"))
names(dat)

# Promote dataframe to spatial
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
dat <- st_as_sf(dat,                         
                coords = c("X", "Y"),
                crs = projcrs)


# Extract values from covariates to the soil points
pv <- terra::extract(x = covs_all, y = vect(dat),xy=F)
dat <- cbind(dat,pv)

summary(dat)

# Remove NA values
dat <-as.data.frame(dat)
dat[, 'geometry'] <-NULL

dat <- dat[complete.cases(dat),]
rm(covs_all)
# Define a formila  and QRF parameters (set 3 times repeated 10-fold cross-validation)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,         ## 10 -fold CV
                           repeats = 3,        ## repeated 3 times
                           verboseIter = TRUE,
                           returnData = TRUE)

fm = as.formula(paste(i," ~", paste0(names(covs),
                                     collapse = "+")))

fm_pc = as.formula(paste(i," ~", paste0(names(covs_pc),
                                        collapse = "+")))
tuneGrid <-  expand.grid(mtry = c(100, 200, 500))


#Calibrate the model using multiple cores
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)



#Correlation covariates
model <- caret::train(fm,
                      data = dat[, c(i,names(covs))],
                      method = "qrf",
                      trControl = fitControl,
                      verbose = TRUE,
                      tuneGrid = tuneGrid,
                      keep.inbag = T)
print(model)
print(paste("qrf model = ", i))

#Extract predictor importance as relative values (%)
model$importance <-
  data.frame(var=rownames(importance(model$finalModel)),
             importance(model$finalModel)/sum(importance(model$finalModel))*100) %>%
  arrange(desc(IncNodePurity))

# print(model)
saveRDS(model, file = paste0("02-Outputs/models/model_",i, ".rds"))

#PC covariates
model_pc <- caret::train(fm_pc,
                      data = dat[, c(i, names(covs_pc))],
                      method = "qrf",
                      trControl = fitControl,
                      verbose = TRUE,
                      tuneGrid = tuneGrid,
                      keep.inbag = T)
print(model_pc)
print(paste("qrf model = ", soilatt[i]))

#Extract predictor importances as relative values (%)
model_pc$importance <-
  data.frame(var=rownames(importance(model_pc$finalModel)),
             importance(model_pc$finalModel)/sum(importance(model_pc$finalModel))*100) %>%
  arrange(desc(IncNodePurity))

# print(model)
saveRDS(model_pc, file = paste0("02-Outputs/models/model_PC_", i, ".rds"))

stopCluster(cl)
gc()


#Quality Assuarance (RMSE,Rsquared,MAE)
model_pc <-readRDS(file = paste0("02-Outputs/models/model_PC_", i, ".rds"))
model <-readRDS(file = paste0("02-Outputs/models/model_", i, ".rds"))


#Select final model (PC vs correlation) based on RMSE,Rsquared,MAE
mod_sel <- data.frame(model=0,model_pc=0)

if(max(model$results$Rsquared)>max(model_pc$results$Rsquared)){
  mod_sel$model <- 1
  mod_sel$model_pc <- 0
}else{
  mod_sel$model <- 0
  mod_sel$model_pc <- 1 }
  if(max(model$results$RMSE)<max(model_pc$results$RMSE)){
    mod_sel$model <-mod_sel$model +1
 }else{
    mod_sel$model_pc <-mod_sel$model_pc +1
}
if(max(model$results$MAE)<max(model_pc$results$MAE)){
    mod_sel$model <-mod_sel$model +1
 }else{
    mod_sel$model_pc <-mod_sel$model_pc +1
}



if(mod_sel$model > mod_sel$model_pc){
  mod_sel <- 'Model based on correlated covaraites'
  final_mod <-model$finalModel
  write.csv(model$results, paste0('02-Outputs/',i,'model_stats.csv'))
  stats <-model$results
  VarImportance <- model$importance
  
  
  rm(covs_pc)
  
  rm(model_pc)
  
  
}else{
  mod_sel <- 'Model based on PCs'
  final_mod <-model_pc$finalModel
  write.csv(model_pc$results, paste0('02-Outputs/',i,'model_pc_stats.csv'))
  stats <-model_pc$results
  VarImportance <- model_pc$importance
  covs <- covs_pc
  rm(covs_pc)
  rm(model)
}

print(paste('Final Model to select:',mod_sel ))
stats
VarImportance

# Generation of maps (prediction of soil attributes)
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
                       # filename = paste0("02-Outputs/Final Maps/tiles/", i,"_tile_", j, ".tif"),
                       # overwrite = TRUE
  )
  pred_sd <- terra::predict(covst, model = final_mod, na.rm=TRUE,  
                     cpkgs="randomForest", what=sd)  
  
  writeRaster(pred_mean, filename = paste0("02-Outputs/tiles/soilatt_tiles/", i,"_tile_", j, ".tif"), 
              overwrite = TRUE)
  writeRaster(pred_sd, filename = paste0("02-Outputs/tiles/soilatt_tiles/", i,"_tileSD_", j, ".tif"), 
              overwrite = TRUE)
  
  rm(pred_mean)
  rm(pred_sd)
  
  
  print(paste("tile",tile[j]))
}

#Merge tiles both prediction and st.Dev
f_mean <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", pattern = paste0(i,"_tile_"), full.names = TRUE)
f_sd <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", pattern =  paste0(i,"_tileSD_"), full.names = TRUE)
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

plot(pred_mean)
plot(pred_sd)

writeRaster(pred_mean, paste0("02-Outputs/maps/", i,"_QRF.tif"),overwrite=TRUE)
writeRaster(pred_sd, paste0("02-Outputs/maps/", i,"_QRF_SD.tif"),overwrite=TRUE)


}


  
