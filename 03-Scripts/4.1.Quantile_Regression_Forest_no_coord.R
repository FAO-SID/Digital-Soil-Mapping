#_______________________________________________________________________________
#
# Quantile Regression Forest
# Soil Property Mapping
#
# GSP-Secretariat
# Contact: Isabel.Luotto@fao.org
#          Marcos.Angelini@fao.org
#_______________________________________________________________________________

#Empty environment and cache 
rm(list = ls())
gc()

# Content of this script =======================================================
# 0 - Set working directory, soil attribute, and packages
# 1 - Merge soil data with environmental covariates 
# 2 - Covariate selection
# 3 - Model calibration
# 4 - Uncertainty assessment
# 5 - Prediction
# 6 - Export final maps
#_______________________________________________________________________________


# 0 - Set working directory, soil attribute, and packages ======================

# Working directory
setwd('C:/GIT/Digital-Soil-Mapping')

# Load Area of interest (shp)
AOI <- '01-Data/buenos_aires.shp'

# Terget soil attribute
soilatt <- 'p_bray'

# Function for Uncertainty Assessment
load(file = "03-Scripts/eval.RData")

#load packages
library(tidyverse)
library(data.table)
library(caret)
library(quantregForest)
library(terra)
library(sf)
library(doParallel)
library(mapview)
# install.packages("terra", repos = "https://rspatial.r-universe.dev")

# 1 - Merge soil data with environmental covariates ============================

val_data <- read_csv("01-Data/data_with_coord.csv") %>% 
  vect(geom=c("x", "y"), crs = "epsg:4326")

covs <- rast("01-Data/covs/covs.tif")
ncovs <- names(covs)

ex <- terra::extract(x = covs, y = val_data, xy=F)
val_data <- as.data.frame(val_data)
val_data <- cbind(val_data, ex)
val_data <- na.omit(val_data)

## 1.1 - Load covariates -------------------------------------------------------
# files <- list.files(path= '01-Data/covs/', pattern = '.tif$', full.names = T)
# covs <- rast(files)
# ncovs <- str_remove(files, "01-Data/covs/")
# ncovs <- str_remove(ncovs, ".tif")
# ncovs <- str_replace(ncovs, "-", "_")
# names(covs) <- ncovs 

# aoi <- vect("01-Data/AOI_Arg.shp")
# mask <- rasterize(aoi, covs[[1]])
# mask <- trim(mask)
# covs <- crop(covs, mask) %>% mask(mask)
# plot(covs)
# writeRaster(covs, "01-Data/covs/covs.tif")



## 1.2 - Load the soil data (Script 2) -----------------------------------------
# dat <- read_csv("02-Outputs/soil_data.csv")

# Convert soil data into a spatial object (check https://epsg.io/6204)

df <- NULL
dat <- read_csv("01-Data/sim_data_no_coord.csv")
N <- dat %>% group_by(stratum) %>% summarise(N=n())

for (i in 1:50) {
  stratum <- rast("01-Data/land cover/SE_districts_croplands.tif")
  names(stratum) <- "stratum"
  
  y <- freq(stratum)
  
  x <- spatSample(x = stratum, max(N$N), method = "stratified", as.df=T, xy=T, replace=F)
  table(x$stratum)
  
  dat <- read_csv("01-Data/sim_data_no_coord.csv")
  
  strata <- unique(x$stratum)
  d <- NULL
  for (j in seq_along(strata)) {
    z <- x[x$stratum==strata[j],]
    n <- N[N$stratum==strata[j],"N"][[1]]
    z <- z[sample(1:max(N$N), size = n),]
    if(strata[j] %in% dat$stratum){
      z <- cbind(z,dat[dat$stratum==strata[j],soilatt])
      d <- rbind(d,z)}
  }
  
  # d %>% 
  #   st_as_sf(coords = c("x", "y"), crs = 4326) %>% # convert to spatial object
  #   mapview(zcol = "p", cex = 2, lwd = 0.1) 
  dat <- vect(d, geom=c("x", "y"), crs = "epsg:4326")
  
  # Reproject point coordinates to match coordinate system of covariates
  # dat <- terra::project(dat, covs)
  # names(dat)
  
  ## 1.3 - Extract values from covariates to the soil points ---------------------
  pv <- terra::extract(x = covs, y = dat, xy=F)
  dat <- cbind(dat,pv)
  dat <- as.data.frame(dat)
  
  
  
  ## 1.4 - Target soil attribute + covariates ------------------------------------
  d <- select(dat, soilatt, names(covs))
  d <- na.omit(d)
  
  
  # 2 - Covariate selection with RFE =============================================
  ## 2.1 - Setting parameters ----------------------------------------------------
  # Repeatedcv = 3-times repeated 10-fold cross-validation
  # fitControl <- rfeControl(functions = rfFuncs,
  #                          method = "cv",
  #                          number = 2,         ## 10 -fold CV
  #                          # repeats = 2,        ## repeated 3 times
  #                          verbose = TRUE,
  #                          saveDetails = TRUE, 
  #                          returnResamp = "all")
  
  # Set the regression function
  fm = as.formula(paste(soilatt," ~", paste0(ncovs,
                                             collapse = "+")))
  
  # Calibrate the model using multiple cores
  # cl <- makeCluster(detectCores()-1)
  # registerDoParallel(cl)
  
  
  ## 2.2 - Calibrate a RFE model to select covariates ----------------------------
  # covsel <- rfe(fm,
  #               data = d,  
  #               sizes = seq(from=10, to=length(ncovs)-1, by = 10),
  #               rfeControl = fitControl,
  #               verbose = TRUE,
  #               keep.inbag = T)
  # stopCluster(cl)
  # 
  # ## 2.3 - Plot selection of covariates ------------------------------------------
  # trellis.par.set(caretTheme())
  # plot(covsel, type = c("g", "o"))
  # 
  # # Extract selection of covariates and subset covs
  # opt_covs <- predictors(covsel)[1:10]
  
  # 3 - QRF Model calibration ====================================================
  ## 3.1 - Update formula with the selected covariates ---------------------------
  # fm <- as.formula(paste(soilatt," ~", paste0(opt_covs, collapse = "+")))
  
  # parallel processing
  print(paste("start fitting model", i))
  timestamp()
  # cl <- makeCluster(detectCores()-1)
  # registerDoParallel(cl)
  
  ## 3.2 - Set training parameters -----------------------------------------------
  fitControl <- trainControl(method = "none",
                             # number = 10,         ## 10 -fold CV
                             # repeats = 3,        ## repeated 3 times
                             savePredictions = TRUE)
  
  # Tune mtry hyperparameters
  mtry <- round(length(ncovs)/3)
  tuneGrid <-  expand.grid(mtry = c(mtry))
  
  ## 3.3 - Calibrate the QRF model -----------------------------------------------
  model <- caret::train(fm,
                        data = d,
                        method = "qrf",
                        trControl = fitControl,
                        verbose = TRUE,
                        tuneGrid = tuneGrid,
                        keep.inbag = T,
                        importance = TRUE, )
  # stopCluster(cl)
  # gc()
  print(paste("end fitting", i))
  timestamp()
  
  ## 3.4 - Extract predictor importance as relative values (%)
  x <- randomForest::importance(model$finalModel)
  model$importance <- x
  ## 3.5 - Print and save model --------------------------------------------------
  print(model)
  saveRDS(model, file = paste0("02-Outputs/models/model_",soilatt,"_",i,".rds"))
  
  # 4 - Uncertainty assessment ===================================================
  # extract observed and predicted values
  o <- val_data$p
  p <- predict(model$finalModel, val_data, what = mean)
  df <- rbind(df, data.frame(o,p, iteration=i))
  
  ## 4.1 - Plot and save scatterplot --------------------------------------------- 
  (g1 <- ggplot(df, aes(x = o, y = p)) + 
     geom_point(alpha = 0.5) + 
     geom_abline(slope = 1, intercept = 0, color = "red")+
     ylim(c(min(o), max(o))) + theme(aspect.ratio=1)+ 
     labs(title = soilatt) + 
     xlab("Observed") + ylab("Predicted"))
  ggsave(g1, filename = paste0("02-Outputs/residuals_",soilatt,"_",i,".png"), 
         scale = 1,
         units = "cm", width = 12, height = 12)
  
  ## 4.2 - Print accuracy coeficients --------------------------------------------
  # https://github.com/AlexandreWadoux/MapQualityEvaluation
  print(paste("model", i))
  print(eval(p,o))
  
  ## 4.3 - Plot Covariate importance ---------------------------------------------
  (g2 <- varImpPlot(model$finalModel, main = soilatt, type = 1))
  
  png(filename = paste0("02-Outputs/importance_",soilatt,"_",i,".png"), 
      width = 15, height = 15, units = "cm", res = 600)
  g2
  dev.off()
}
write_csv(df, paste("residuals_",soilatt,".csv"))
# 5 - Prediction ===============================================================
# Generation of maps (prediction of soil attributes) 
## 5.1 - Produce tiles ---------------------------------------------------------
r <-covs[[1]]
t <- rast(nrows = 5, ncols = 5, extent = ext(r), crs = crs(r))
tile <- makeTiles(r, t,overwrite=TRUE,filename="02-Outputs/tiles/tiles.tif")

## 5.2 - Predict soil attributes per tiles -------------------------------------
# loop to predict on each tile

for (j in seq_along(tile)) {
  for (i in 1:50) {
    
    gc()
    t <- rast(tile[j])
    covst <- crop(covs, t)
    
    model <- readRDS(paste0("02-Outputs/models/model_",soilatt,"_", i,".rds"))
    # plot(r)# 
    meanFunc = function(x) mean(x, na.rm=TRUE)
    pred_mean <- terra::predict(covst, model = model$finalModel, na.rm = TRUE,   
                                cpkgs="quantregForest", what=meanFunc)
    sdFunc = function(x) sd(x, na.rm=TRUE)
    pred_sd <- terra::predict(covst, model = model$finalModel, na.rm=TRUE,  
                              cpkgs="quantregForest", what=sdFunc)  
    
    writeRaster(pred_mean, 
                filename = paste0("02-Outputs/tiles/soilatt_tiles/",
                                  soilatt,"_tile_", j,"_model_",i,".tif"), 
                overwrite = TRUE)
    writeRaster(pred_sd, 
                filename = paste0("02-Outputs/tiles/soilatt_tiles/",
                                  soilatt,"_tileSD_", j,"_model_",i,".tif"), 
                overwrite = TRUE)
    
    rm(pred_mean)
    rm(pred_sd)
    
    
    print(paste("tile",tile[j],"model", i))
  }
}

## 5.3 - Merge tiles both prediction and st.Dev --------------------------------
for (j in 1:25) {
  f_mean <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", 
                       pattern = paste0(soilatt,"_tile_",j,"_model_"), full.names = TRUE)
  f_sd <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", 
                     pattern =  paste0(soilatt,"_tileSD_",j,"_model_"), full.names = TRUE)
  pred_mean <- rast(f_mean) %>% mean(na.rm=TRUE)
  # q05 <- rast(f_mean) %>% quantile(.05)
  # q95 <- rast(f_mean) %>% quantile(.95)
  r_sd_l <- rast(f_sd) %>% mean(na.rm=TRUE)
  names(r_sd_l) <- "sd"
  
  writeRaster(r_sd_l, 
              filename = paste0("02-Outputs/tiles/aggregated_tiles/",
                                soilatt,"_tilesd_", j,".tif"), 
              overwrite = TRUE)
  writeRaster(q05, 
              filename = paste0("02-Outputs/tiles/aggregated_tiles/",
                                soilatt,"_tileq05_", j,".tif"), 
              overwrite = TRUE)
  writeRaster(q50, 
              filename = paste0("02-Outputs/tiles/aggregated_tiles/",
                                soilatt,"_tileq50_", j,".tif"), 
              overwrite = TRUE)
  writeRaster(q95, 
              filename = paste0("02-Outputs/tiles/aggregated_tiles/",
                                soilatt,"_tileq95_", j,".tif"), 
              overwrite = TRUE)
}

name_tiles <- c("sd_", "q05_","q50_", "q95_")

for (i in seq_along(name_tiles)) {
  f <- list.files(path = "02-Outputs/tiles/aggregated_tiles/", 
                  pattern = name_tiles[i], full.names = TRUE)
  r <- list()
  for (g in seq_along(f)){
    r[[g]] <- rast(f[g])
    print(g)
  }
  r <- sprc(r)
  r <- mosaic(r)
  plot(r)
  writeRaster(r, paste0("02-Outputs/maps/",soilatt,name_tiles[i],".tif"),
              overwrite=TRUE)
  rm(r)
  gc()
}

val_data <- read_csv("01-Data/data_with_coord.csv") %>% 
  vect(geom=c("x", "y"), crs = "epsg:4326")

q50 <- rast("02-Outputs/maps/p_tileq50.tif")

ex <- terra::extract(x = q50, y = val_data, xy=F, ID=FALSE)
val_data <- as.data.frame(val_data)
val_data <- cbind(val_data, ex)
val_data <- na.omit(val_data)

eval(val_data$q0.50,val_data$p)

##############################

v <- st_read("01-Data/AOI_Arg.shp")
names(v)[5] <- "stratum"
darea <- readxl::read_excel("01-Data/data_without_coord.xls")
darea <-  darea %>% 
  select(stratum=code, median)
darea$stratum <- as.numeric(darea$stratum)
v <- left_join(v, darea)
v <- vect(v)
names(v)[8] <- "p"
plot(v, "p")
r <- rast("01-Data/land cover/SE_croplands.tif")
r <- rasterize(v,r,field="p")
plot(r)
val_data <- read_csv("01-Data/data_with_coord.csv") %>% 
  vect(geom=c("x", "y"), crs = "epsg:4326")

val_data$pred <- terra::extract(r, val_data, ID=FALSE)[[1]]

val <- val_data %>% as.data.frame() %>% na.omit()

eval(val$pred, val$p)

##################


val_data <- read_csv("01-Data/data_with_coord.csv") %>% 
  vect(geom=c("x", "y"), crs = "epsg:4326")

covs <- rast("01-Data/covs/covs.tif")
psd <- rast("02-Outputs/maps/p_tileSD.tif")
p05 <- rast("02-Outputs/maps/p_tileq05.tif")
p50 <- rast("02-Outputs/maps/p_tileq50.tif")
p95 <- rast("02-Outputs/maps/p_tileq95.tif")
covs <- c(covs, psd, p05, p50, p95)

ncovs <- names(covs)

ex <- terra::extract(x = covs, y = val_data, xy=F)
val_data <- as.data.frame(val_data)
val_data <- cbind(val_data, ex)
val_data <- na.omit(val_data)

plot(covs[[71]])

## 1.4 - Target soil attribute + covariates ------------------------------------
d <- select(val_data, soilatt, ncovs)
d <- na.omit(d)


# 2 - Covariate selection with RFE =============================================
## 2.1 - Setting parameters ----------------------------------------------------
# Repeatedcv = 3-times repeated 10-fold cross-validation
fitControl <- rfeControl(functions = rfFuncs,
                         method = "repeatedcv",
                         number = 10,         ## 10 -fold CV
                         repeats = 3,        ## repeated 3 times
                         verbose = TRUE,
                         saveDetails = TRUE, 
                         returnResamp = "all")

# Set the regression function
fm = as.formula(paste(soilatt," ~", paste0(ncovs,
                                           collapse = "+")))

# Calibrate the model using multiple cores
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)


## 2.2 - Calibrate a RFE model to select covariates ----------------------------
covsel <- rfe(fm,
              data = d,  
              sizes = seq(from=10, to=length(ncovs)-1, by = 5),
              rfeControl = fitControl,
              verbose = TRUE,
              keep.inbag = T)
stopCluster(cl)

## 2.3 - Plot selection of covariates ------------------------------------------
trellis.par.set(caretTheme())
plot(covsel, type = c("g", "o"))

# Extract selection of covariates and subset covs
opt_covs <- predictors(covsel)[1:65]

# 3 - QRF Model calibration ====================================================
## 3.1 - Update formula with the selected covariates ---------------------------
fm <- as.formula(paste(soilatt," ~", paste0(opt_covs, collapse = "+")))

# parallel processing
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

## 3.2 - Set training parameters -----------------------------------------------
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,         ## 10 -fold CV
                           repeats = 3,        ## repeated 3 times
                           savePredictions = TRUE)

# Tune mtry hyperparameters
mtry <- round(length(opt_covs)/3)
tuneGrid <-  expand.grid(mtry = c(mtry-5, mtry, mtry+5))

## 3.3 - Calibrate the QRF model -----------------------------------------------
model <- caret::train(fm,
                      data = d,
                      method = "qrf",
                      trControl = fitControl,
                      verbose = TRUE,
                      tuneGrid = tuneGrid,
                      keep.inbag = T,
                      importance = TRUE)
stopCluster(cl)
gc()


## 3.4 - Extract predictor importance as relative values (%)
x <- randomForest::importance(model$finalModel)
model$importance <- x
## 3.5 - Print and save model --------------------------------------------------
print(model)
saveRDS(model, file = paste0("02-Outputs/models/model_",soilatt,".rds"))

# 4 - Uncertainty assessment ===================================================
# extract observed and predicted values
o <- model$pred$obs
p <- model$pred$pred
df <- data.frame(o,p)

## 4.1 - Plot and save scatterplot --------------------------------------------- 
(g1 <- ggplot(df, aes(x = o, y = p)) + 
   geom_point(alpha = 0.3) + 
   geom_abline(slope = 1, intercept = 0, color = "red")+
   ylim(c(min(o), max(o))) + theme(aspect.ratio=1)+ 
   labs(title = soilatt) + 
   xlab("Observed") + ylab("Predicted"))
ggsave(g1, filename = paste0("02-Outputs/residuals_",soilatt,".png"), scale = 1, 
       units = "cm", width = 12, height = 12)

## 4.2 - Print accuracy coeficients --------------------------------------------
# https://github.com/AlexandreWadoux/MapQualityEvaluation
eval(p,o)

## 4.3 - Plot Covariate importance ---------------------------------------------
(g2 <- varImpPlot(model$finalModel, main = soilatt, type = 1))

png(filename = paste0("02-Outputs/importance_",soilatt,".png"), 
    width = 15, height = 15, units = "cm", res = 600)
g2
dev.off()

# 5 - Prediction ===============================================================
# Generation of maps (prediction of soil attributes) 
## 5.1 - Produce tiles ---------------------------------------------------------
r <-covs[[1]]
t <- rast(nrows = 5, ncols = 5, extent = ext(r), crs = crs(r))
tile <- makeTiles(r, t,overwrite=TRUE,filename="02-Outputs/tiles/tiles.tif")

## 5.2 - Predict soil attributes per tiles -------------------------------------
# loop to predict on each tile

for (j in seq_along(tile)) {
  gc()
  t <- rast(tile[j])
  covst <- crop(covs, t)
  
  
  # plot(r)# 
  pred_mean <- terra::predict(covst, model = model$finalModel, na.rm=TRUE,  
                              cpkgs="quantregForest", what=mean)
  pred_sd <- terra::predict(covst, model = model$finalModel, na.rm=TRUE,  
                            cpkgs="quantregForest", what=sd)  
  
  
  
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
  
  writeRaster(pred_mean, 
              filename = paste0("02-Outputs/tiles/soilatt_tiles/",
                                soilatt,"_tile_", j, ".tif"), 
              overwrite = TRUE)
  writeRaster(pred_sd, 
              filename = paste0("02-Outputs/tiles/soilatt_tiles/",
                                soilatt,"_tileSD_", j, ".tif"), 
              overwrite = TRUE)
  
  rm(pred_mean)
  rm(pred_sd)
  
  
  print(paste("tile",tile[j]))
}

## 5.3 - Merge tiles both prediction and st.Dev --------------------------------
f_mean <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", 
                     pattern = paste0(soilatt,"_tile_"), full.names = TRUE)
f_mean <- f_mean[str_detect(f_mean, "model", negate = TRUE)]
f_sd <- list.files(path = "02-Outputs/tiles/soilatt_tiles/", 
                   pattern =  paste0(soilatt,"_tileSD_"), full.names = TRUE)
f_sd <- f_sd[str_detect(f_sd, "model", negate = TRUE)]
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


# 6 - Export final maps ========================================================
writeRaster(pred_mean, 
            paste0("02-Outputs/maps/",soilatt,"_xy_QRF.tif"),
            overwrite=TRUE)
writeRaster(pred_sd, 
            paste0("02-Outputs/maps/",soilatt,"_xy_QRF_SD.tif"),
            overwrite=TRUE)

