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
wd <- 'C:/Users/luottoi/Documents/GitHub/Digital-Soil-Mapping'
#wd <- 'C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping'
#List of soil attributes prepared in script #2
soilatt <- c('OCS','clay', 'pH')

AOI <- '01-Data/MKD.shp'
#
#
#######################################################
# Set working directory
setwd(wd)
# load packages

# install packages (if not installed before)
#install.packages(c("quantregForest", "snow", "Metrics"))

############################### Prapare the final table for modelling (regression matrix) 
library(tidyverse)
library(data.table)
library(caret)
library(quantregForest)
library(terra)
library(sf)
library(doParallel)

# Load PC covariates (the same for all soil attributes)
covs_pc <- rast('02-Outputs/PCA_covariates.tif')

# Calibrate QRF models, perform QA and select the best model ----
for(i in unique(soilatt)){

#Load data and covariates stacks based on correlation
covs <-rast( paste0("02-Outputs/", i,"_covariates.tif"))
  
covs_all <-  c(covs, covs_pc)  
# Load the processed data for digital soil mapping. This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv(paste0("02-Outputs/",i,"_dat_train.csv"))
names(dat)

# extract values from covariates to the soil points
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
dat <- st_as_sf(dat,                         
                coords = c("X", "Y"),
                crs = projcrs)


# Extract values from covariates to the soil points
pv <- extract(x = covs_all, y = vect(dat),xy=F)
dat <- cbind(dat,pv)

summary(dat)

# Remove NA values
dat <-as.data.frame(dat)
dat[, 'geometry'] <-NULL

dat <- dat[complete.cases(dat),]
rm(covs_all)
# Defined formulas and QRF parameters
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
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
# set 3 times repeated 10-fold cross-validation


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

#Extract predictor importances as relative values (%)
model$importance <-
  data.frame(var=rownames(importance(model$finalModel)),
             importance(model$finalModel)/sum(importance(model$finalModel))*100) %>%
  arrange(desc(IncNodePurity))

# print(model)
saveRDS(model, file = paste0("02-Outputs/model_",i, ".rds"))

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
saveRDS(model_pc, file = paste0("02-Outputs/model_PC_", i, ".rds"))

stopCluster(cl)
gc()


#Quality Assuarance (RMSE,Rsquared,MAE)
model_pc <-readRDS(file = paste0("02-Outputs/model_PC_", i, ".rds"))
model <-readRDS(file = paste0("02-Outputs/model_", i, ".rds"))


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
  rm(model_pc)
}else{
  mod_sel <- 'Model based on PCs'
  final_mod <-model_pc$finalModel
  rm(model)
}

print(paste('Final Model to select:',mod_sel ))
final_mod$results




}


  
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
  pred_mean <- predict(covst, model = final_mod, na.rm=TRUE,  
                       cpkgs="randomForest", what=mean,
                      # filename = paste0("02-Outputs/Final Maps/tiles/", i,"_tile_", j, ".tif"),
                      # overwrite = TRUE
                      )
  pred_sd <- predict(covst, model = final_mod, na.rm=TRUE,  
                     cpkgs="randomForest", what=sd)  
  
  writeRaster(pred_mean, filename = paste0("02-Outputs/Final Maps/tiles/", i,"_tile_", j, ".tif"), 
              overwrite = TRUE)
  writeRaster(pred_sd, filename = paste0("02-Outputs/Final Maps/tiles/", i,"_tileSD", j, ".tif"), 
              overwrite = TRUE)
            
  rm(pred_mean)
  rm(pred_sd)

  
  print(paste("tile",tile[j]))
}

#Merge tiles both prediction and st.Dev
f_mean <- list.files(path = "02-Outputs/Final Maps/tiles/", pattern = "tile_", full.names = TRUE)
f_sd <- list.files(path = "02-Outputs/Final Maps/tiles/", pattern = "tileSD", full.names = TRUE)
r_mean_l <- list()
r_sd_l <- list()
for (i in 1:length(f_mean)){
  
  r <- rast(f_mean[i])
  r_mean_l[i] <-r
rm(r)
}

for (i in 1:length(f_sd)){
  
  r <- rast(f_sd[i])
  r_sd_l[i] <-r
  rm(r)
}

r_mean <-sprc(r_mean_l)
r_sd <-sprc(r_sd_l)
pred_mean <- mosaic(r_mean)
pred_sd <- mosaic(r_sd)

plot(pred_mean)
plot(pred_sd)
