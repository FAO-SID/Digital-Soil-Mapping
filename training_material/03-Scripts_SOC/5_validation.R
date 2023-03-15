# This script is a part of the FAO/GSP training course on Digital Soil Mapping.
# It is used for the validation of the maps.

# Remove existing objects in memory and clear memory cache
rm(list = ls());
gc()

# Set working directory
setwd("E:/GSP/Trainings/Indonesia/Macedonia")

####################################### Calculate validation statistics #############################################
library(raster)

# Load and stack the maps from the results folder
RKmap<-raster("02-Outputs/Final Maps/MKD_OCS_RK.tif")
RFmap<-raster("02-Outputs/Final Maps/MKD_OCS_RF.tif")
maps <- stack(RKmap, RFmap)

# Explore the maps
names(maps)
summary(maps)
plot(maps)

# Load the validation dataset. It was was prepared in the 'data_preparation_profiles' script
test <- read.csv("02-Outputs/dat_test.csv")

# Promote to spatialPointsDataFrame and set crs
coordinates(test) <- ~ X + Y
test@proj4string <- CRS(projargs = "+init=epsg:4326")

# Extract the predicted values from the maps to the validation dataset
test <- extract(x = maps, y = test, sp = TRUE)
summary(test)

# Remove NA values
test<-as.data.frame(test)
test <- test[complete.cases(test),]

# Calculate prediction errors
test$PE_RK <- test$MKD_OCS_RK - test$OCS
test$PE_RF <- test$MKD_OCS_RF - test$OCS

# Explore prediction errors
summary(test$PE_RK)
summary(test$PE_RF)

hist(test$PE_RK, breaks=50)
hist(test$PE_RF, breaks=50)

# Regression Kriging 

# Mean Error
ME_RK <- mean(test$PE_RK)
# Mean Absolute Error (MAE)
MAE_RK <- mean(abs(test$PE_RK))
# Root Mean Squared Error (RMSE)
RMSE_RK <- sqrt( sum(test$PE_RK^2) / length(test$PE_RK) )
# Amount of Variance Explained (AVE)
AVE_RK <- 1 - sum(test$PE_RK^2) / sum( (test$OCS - mean(test$OCS))^2 )

# Random Forest

# Mean Error
ME_RF <- mean(test$PE_RF)
# Mean Absolute Error (MAE)
MAE_RF <- mean(abs(test$PE_RF))
# Root Mean Squared Error (RMSE)
RMSE_RF <- sqrt(sum(test$PE_RF^2) / length(test$PE_RF))
# Amount of Variance Explained (AVE)
AVE_RF <- 1 - sum(test$PE_RF^2) / sum( (test$OCS - mean(test$OCS))^2 )

# scatter plot
plot(test$MKD_OCS_RK, test$OCS, main="Regression Kriging", xlab="predicted",
     ylab='observed')
# 1:1 line in black
abline(0,1, lty=2, col='black')
# regression line between predicted and observed in blue
abline(lm(test$OCS ~ test$MKD_OCS_RK), col = 'blue', lty=2)

# scatter plot
plot(test$MKD_OCS_RF, test$OCS, main="Random Forest", xlab="predicted",
     ylab='observed')
# 1:1 line in black
abline(0,1, lty=2, col='black')
# regression line between predicted and observed in blue
abline(lm(test$OCS ~ test$MKD_OCS_RF), col = 'blue', lty=2)


