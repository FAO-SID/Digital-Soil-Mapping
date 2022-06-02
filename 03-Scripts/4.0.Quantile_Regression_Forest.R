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
soilatt <- c('SOC','clay', 'pH')
#
#
#######################################################
# Set working directory
setwd(wd)
# load packages

# install packages (if not installed before)
#install.packages(c("quantregForest", "snow", "Metrics"))

############################### Prapare the final table for modelling (regression matrix) 
library(raster)
library(sp)
library(caret)
library(Metrics)
library(quantregForest)
library(snow)
library(foreach)    
library(doParallel)


# Generate RF maps and uncertainty with the selected covariates ----
for(i in unique(soilatt)){

# Load the covariates stack. It was was prepared in the 'data_preparation_covariates' script

load(file = paste0("02-Outputs/", i,"_covariates.RData"))
names(covs)

# Load the processed data for digital soil mapping. This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv(paste0("02-Outputs/",i,"_dat_train.csv"))
names(dat)

# extract values from covariates to the soil points
coordinates(dat) <- ~ X + Y
dat <- extract(x = covs, y = dat, sp = TRUE)
summary(dat)

# Remove NA values
dat<-as.data.frame(dat)
dat <- dat[complete.cases(dat),]

str(dat)

# LandCover and soilmap are categorical variables, they need to be 'factor' type
#dat$LandCover <- as.factor(dat$LandCover)
#dat$soilmap <- as.factor(dat$soilmap)
str(dat)

#coordinates(dat) <- ~ X + Y

# Generate an empty dataframe
#validation <- data.frame(rmse=numeric(), r2=numeric())

fm = as.formula(paste(i," ~", paste0(names(covs),
                                      collapse = "+")))
fm


ctrl <- trainControl(method = "cv", savePred=T)


# Sensitivity to the dataset
# Start a loop with 10 model realizations
registerDoParallel(makeCluster(4)) # Use 4 cores for parallel CV

# Assuming this is your dataset 

cv <- caret::createFolds(1:nrow(dat), k=10, list=T) # Create 10 folds

pred <- stack()
# 'dopar' here would run this on multiple threads (change to just 'do' for synchronous runs)
results <- foreach(i= 1:length(cv), .packages = c("raster", "sp", "caret"),.export = ls(globalenv())) %dopar% {
  
  # Get the fold data where 'fold' is nothing more than a list of indexes for test observations in the data
  train <- dat[cv[[i]],] # Get the opposite of the test observations to train on
  
  
  # Fit the model and make predictions
  modn <- train(fm, data=train, method = "rf",
                trControl = ctrl)
  pred_cv <-predict(covs, modn)
  stack(pred, pred_cv)
  
}
sensitivity <- calc(stack(results), sd)

# The sensitivity map is the dispersion of all individual models
plot(sensitivity, col=rev(topo.colors(10)),
     main='Sensitivity based on 10 realizations using 25% samples')



#generate mtry variable using caret package 
rfmodel<- train(fm, data=dat, method = "rf", trControl = ctrl,
                importance=TRUE)


# run the quantile regression forest algorithm
model <- quantregForest(y=dat[[i]], x=dat[,names(covs)],
                        ntree=500, keep.inbag=TRUE,
                        mtry = as.numeric(rfmodel$bestTune))


# Define number of cores to use
beginCluster()

# Estimate model uncertainty
unc <- clusterR(covs, predict, args=list(model=model,what=sd))

# OCS prediction based in all available data
mean <- clusterR(covs, predict, args=list(model=model,what=mean))
# The total uncertainty is the sum of sensitivity and model
# uncertainty
unc <- unc + sensitivity
# Express the uncertainty in percent % (divide by the mean)
Total_unc_Percent <- unc/mean*100
endCluster()

# Plot both maps (the predicted OCS and associated uncertainty)
plot(mean, main='based in all data')
plot(Total_unc_Percent, main='Total uncertainty %')

writeRaster(mean, paste0('02-Outputs/Final Maps/',i,'_RF.tif'), overwrite=TRUE)
writeRaster(unc, paste0('02-Outputs/Final Maps/',i,'_RF_sd.tif'), overwrite=TRUE)


# Run analysis with PCs
#Load PCs
covs_pc <- stack('02-Outputs/PCA_covariates.tif')
# Load the processed data for digital soil mapping. This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv(paste0("02-Outputs/",i,"_dat_train.csv"))
names(dat)

# extract values from covariates to the soil points
coordinates(dat) <- ~ X + Y
dat <- extract(x = covs_pc, y = dat, sp = TRUE)
summary(dat)

# Remove NA values
dat<-as.data.frame(dat)
dat <- dat[complete.cases(dat),]

str(dat)





fm = as.formula(paste(i," ~", paste0(names(covs_pc),
                                     collapse = "+")))
fm


# Sensitivity to the dataset
# Start a loop with 10 model realizations
registerDoParallel(makeCluster(4)) # Use 4 cores for parallel CV

# Assuming this is your dataset 

cv <- caret::createFolds(1:nrow(dat), k=10, list=T) # Create 10 folds

pred <- stack()
# 'dopar' here would run this on multiple threads (change to just 'do' for synchronous runs)
results <- foreach(i= 1:length(cv), .packages = c("raster", "sp", "caret"),.export = ls(globalenv())) %dopar% {
  
  # Get the fold data where 'fold' is nothing more than a list of indexes for test observations in the data
  train <- dat[cv[[i]],] # Get the opposite of the test observations to train on
  
  
  # Fit the model and make predictions
  modn <- train(fm, data=train, method = "rf",
                trControl = ctrl)
  pred_cv <-predict(covs_pc, modn)
  stack(pred, pred_cv)
  
}
sensitivity <- calc(stack(results), sd)
plot(sensitivity, col=rev(topo.colors(10)),
     main='Sensitivity based on 10 realizations using 25% samples')



#generate mtry variable using caret package 
rfmodel<- train(fm, data=dat, method = "rf", trControl = ctrl,
                importance=TRUE)


# run the quantile regression forest algorithm
model <- quantregForest(y=dat[[i]], x=dat[,names(covs_pc)],
                        ntree=500, keep.inbag=TRUE,
                        mtry = as.numeric(rfmodel$bestTune))


# Define number of cores to use
beginCluster()
## 8 cores detected, using 7
# Estimate model uncertainty
unc <- clusterR(covs_pc, predict, args=list(model=model,what=sd))

# OCS prediction based in all available data
mean <- clusterR(covs_pc, predict, args=list(model=model,what=mean))
# The total uncertainty is the sum of sensitivity and model
# uncertainty
unc <- unc + sensitivity
# Express the uncertainty in percent % (divide by the mean)
Total_unc_Percent <- unc/mean*100
endCluster()

# Plot both maps (the predicted OCS and associated uncertainty)
plot(mean, main='based in all data')
plot(Total_unc_Percent, main='Total uncertainty %')

writeRaster(mean, paste0('02-Outputs/Final Maps/',i,'_RF_PCA.tif'), overwrite=TRUE)
writeRaster(unc, paste0('02-Outputs/Final Maps/',i,'_RF_sd_PCA.tif'), overwrite=TRUE)




}
