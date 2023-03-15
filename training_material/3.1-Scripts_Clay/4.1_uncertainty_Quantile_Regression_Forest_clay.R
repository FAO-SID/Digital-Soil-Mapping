# This script is a part of the FAO/GSP training course on Digital Soil Mapping.
# It is used for uncertainty assessment of the Random Forest map.
# Please be aware that the algorithm is resource demanding, 
# and may require a lot of time and memory to run.

# Remove existing objects in memory and clear memory cache
rm(list = ls());
gc()

# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

# install packages (if not installed before)
#install.packages(c("quantregForest", "snow", "Metrics"))

############################### Prapare the final table for modelling (regression matrix) ########################
library (raster)
library (sp)

# Load the covariates stack. It was was prepared in the 'data_preparation_covariates' script
load(file = "02-Outputs/covariates_clay.RData")
names(covs)


# Load the processed data for digital soil mapping. This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv("02-Outputs/dat_train_clay.csv")
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
dat$soilmap <- as.factor(dat$soilmap)
dat$LandCover <- as.factor(dat$LandCover)

str(dat)

coordinates(dat) <- ~ X + Y


####################################### Calculate validation statistics #############################################

library(Metrics)
library(ModelMetrics)
library (caret)

set.seed(12345)

# Generate an empty dataframe
validation <- data.frame(rmse=numeric(), r2=numeric())

fm = as.formula(paste("clay ~", paste0(names(covs),
                                         collapse = "+")))
fm

ctrl <- trainControl(method = "cv", savePred=T)


# Search for the best mtry parameter
rfmodel <- train(fm, data=dat@data, method = "rf", trControl = ctrl,
                 importance=TRUE)
pred <- predict(covs, rfmodel)

# Sensitivity to the dataset
# Start a loop with 10 model realizations
for (i in 1:10){
  # We will build 10 models using random samples of 25%
  smp_size <- floor(0.25 * nrow(dat))
  train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
  train <- dat[train_ind, ]
  test <- dat[-train_ind, ]
  modn <- train(fm, data=train@data, method = "rf",
                trControl = ctrl)
  pred <- stack(pred, predict(covs, modn))
  test$pred <- extract(pred[[i+1]], test)
  # Store the results in a dataframe
  validation[i, 1] <- rmse(test$clay, test$pred)
  validation[i, 2] <- cor(test$clay, test$pred)^2
}

# The sensitivity map is the dispersion of all individual models
sensitivity <- calc(pred[[-1]], sd)
plot(sensitivity, col=rev(topo.colors(10)),
     main='Sensitivity based on 10 realizations using 25% samples')

summary(validation)

# Plot of the map based on 75% of data and the sensitivity to data
# variations
#prediction75 <- exp(pred[[1]])
plot(pred[[1]], main='clay prediction based on 75% of data')

# Use quantile regression forest to estimate the full conditional
# distribution of clay, note that we are using the mtry
# parameter that was selected by the train function of the caret
# package, assuming that the 75% of data previously used well
# resembles the statistical distribution of the entire data
# population. Otherwise, repeat the train function with all
# available data (using the object dat that instead of train)
# to select mtry.
library(quantregForest)
names(dat@data)
model <- quantregForest(y=dat@data$clay, x=dat@data[,7:16],
                        ntree=500, keep.inbag=TRUE,
                        mtry = as.numeric(rfmodel$bestTune))

library(snow)
# Estimate model uncertainty at the pixel level using parallel
# computing
# Define number of cores to use
beginCluster()
## 8 cores detected, using 7
# Estimate model uncertainty
unc <- clusterR(covs, predict, args=list(model=model,what=sd))

mean <- clusterR(covs, predict, args=list(model=model,what=mean))
# The total uncertainty is the sum of sensitivity and model
# uncertainty
unc <- unc + sensitivity
# Express the uncertainty in percent % (divide by the mean)
Total_unc_Percent <- unc/mean
endCluster()

# Plot both maps (the predicted clay and associated uncertainty)
plot(mean, main='clay based in all data')
plot(unc, main='Total uncertainty')
plot(Total_unc_Percent, main='Total uncertainty in %')

writeRaster(unc, '02-Outputs/Final Maps/MKD_clay_RF_sd.tif', overwrite=TRUE)
