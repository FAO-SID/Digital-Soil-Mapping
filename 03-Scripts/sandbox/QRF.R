# target soil properties
library(tidyverse)
library(data.table)
library(caret)
library(quantregForest)
library(terra)
library(doParallel)

setwd("/home/nomade/Documents/gbsmap/Global/")

# calibration data with soil properties and covariates
d <- read_csv("file.csv")

target <- c("soc", "ph", "cec")
# names of covariates
cov_names <- names(d)[4:30] # check col numbers


for (i in seq_along(target)) {
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  # set 3 times repeated 10-fold cross-validation
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,         ## 10 -fold CV
                             repeats = 3,        ## repeated 3 times
                             verboseIter = TRUE,
                             returnData = TRUE)
  # tune mtry
  # tuneGrid <-  expand.grid(mtry = c(100, 200, 500))
  fm <- as.formula(paste(target[i], '~', paste0(cov_names, collapse = "+")))
  model <- caret::train(fm,
                        data = d,
                        method = "qrf",
                        trControl = fitControl,
                        verbose = TRUE,
                        tuneGrid = tuneGrid,
                        keep.inbag = T)
  print(model)
  print(paste("qrf model = ", target[i]))

  #Extract predictor importances as relative values (%)
  model$importance <-
    data.frame(var=rownames(importance(model$finalModel)),
               importance(model$finalModel)/sum(importance(model$finalModel))*100) %>%
    arrange(desc(IncNodePurity))

  # print(model)
  saveRDS(model, file = paste0("models/model_", target[i], ".rds"))
  stopCluster(cl = cl)
  gc()
}

# predict at tiles #-----------------------------------------------------------#
library(terra)
terraOptions(tempdir = "/home/nomade/tmp", todisk = TRUE)
library(tidyverse)
library(data.table)
library(caret)
# library(ranger)
# library(tuneRanger)
library(quantregForest)
setwd("/home/nomade/Documents/gbsmap/Global/")
rfmodel <- readRDS(paste0("models/model_", "cec", ".rds"))

# Load shape file with tiles for prediction
tile <- vect("tiles/tiles_small.shp")
tile$lyr.1
# load tif covariates
f <- list.files(path = "covs", pattern = ".tif$", full.names = TRUE)

# loop to predict on each tile
for (j in seq_along(tile)) {
  t <- tile[j]
  r <- list()
  gc()
  # build stack for each tile
  for (i in seq_along(f)) {
    x <- rast(f[i])
    names(x) <- str_remove(str_remove(f[i], ".tif"), "covs/")
    x <- crop(x, t)
    if (i==1) {
      g <- x
    }
    r[[i]] <- project(x, g) %>% mask(g)
    print(i)
    gc()
  }
  r <- rast(r)
  gc()
  # plot(r)# 
  pred_mean <- predict(r, model = model$finalModel, na.rm=TRUE,  
                  cpkgs="randomForest", what="mean",
                  filename = paste0("pred/", k, "_tile_", tile$lyr.1[j], ".tif"),
                  overwrite = TRUE)
  pred_sd <- predict(r, model = model$finalModel, na.rm=TRUE,  
                  cpkgs="randomForest", what="sd",
                  filename = paste0("pred/", "cec", "_tile_", tile$lyr.1[j], ".tif"),
                  overwrite = TRUE)                
  writeRaster(pred_mean, filename = paste0("pred/", "cec_mean", "_tile_", tile$lyr.1[j], ".tif"), 
              overwrite = TRUE)
  writeRaster(pred_sd, filename = paste0("pred/", "cec_sd", "_tile_", tile$lyr.1[j], ".tif"), 
              overwrite = TRUE)
  rm(pred_mean)
  rm(pred_sd)
  gc()
  print(paste("tile",tile$lyr.1[j]))
}
