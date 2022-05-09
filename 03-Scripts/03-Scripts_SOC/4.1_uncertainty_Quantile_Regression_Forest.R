# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

# install packages (if not installed before)
#install.packages(c("quantregForest", "snow", "Metrics"))

############################### Prapare the final table for modelling (regression matrix) 
library (raster)
library (sp)

# Load the covariates stack. It was was prepared in the 'data_preparation_covariates' script
load(file = "02-Outputs/covariates.RData")
names(covs)

# Load the processed data for digital soil mapping. This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv("02-Outputs/dat_train.csv")
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
str(dat)

coordinates(dat) <- ~ X + Y

# Generate an empty dataframe
validation <- data.frame(rmse=numeric(), r2=numeric())

fm = as.formula(paste("OCS ~", paste0(names(covs),
                                      collapse = "+")))
fm

library(caret)
library(Metrics)
ctrl <- trainControl(method = "cv", savePred=T)


# Sensitivity to the dataset
# Start a loop with 10 model realizations
pred <- stack()
for (i in 1:10){
  # We will build 10 models using random samples of 25%
  smp_size <- floor(0.25 * nrow(dat))
  train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
  train <- dat[train_ind, ]
  test <- dat[-train_ind, ]
  modn <- train(fm, data=train@data, method = "rf",
                trControl = ctrl)
  pred_cv <-predict(covs, modn)
  pred <- stack(pred, pred_cv)
  test$pred <- extract(pred[[i]], test)
  # Store the results in a dataframe
  validation[i, 1] <- rmse(test$OCS, test$pred)
  validation[i, 2] <- cor(test$OCS, test$pred)^2
}

# The sensitivity map is the dispersion of all individual models
sensitivity <- calc(pred, sd)
plot(sensitivity, col=rev(topo.colors(10)),
     main='Sensitivity based on 10 realizations using 25% samples')

summary(validation)

#generate mtry variable using caret package 
rfmodel<- train(fm, data=dat@data, method = "rf", trControl = ctrl,
                importance=TRUE)

library(quantregForest)
# run the quantile regression forest algorithm
model <- quantregForest(y=dat@data$OCS, x=dat@data[,8:19],
                        ntree=500, keep.inbag=TRUE,
                        mtry = as.numeric(rfmodel$bestTune))

library(snow)
# Define number of cores to use
beginCluster()
## 8 cores detected, using 7
# Estimate model uncertainty
unc <- clusterR(covs, predict, args=list(model=model,what=sd))

# OCS prediction based in all available data
mean <- clusterR(covs, predict, args=list(model=model,what=mean))
# The total uncertainty is the sum of sensitivity and model
# uncertainty
unc <- unc + sensitivity
# Express the uncertainty in percent % (divide by the mean)
Total_unc_Percent <- unc/mean
endCluster()

# Plot both maps (the predicted OCS and associated uncertainty)
plot(mean, main='OCS based in all data')
plot(unc, main='Total uncertainty')
plot(Total_unc_Percent, main='Total uncertainty in %')

writeRaster(unc, '02-Outputs/Final Maps/MKD_OCS_RF_sd.tif', overwrite=TRUE)
