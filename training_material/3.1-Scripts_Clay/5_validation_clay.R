# Set working directory
rm(list = ls())
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")
library(raster)


# Load and stack the maps from the results folder
RKmap<-raster("02-Outputs/Final Maps/MKD_clay_RK.tif")
RFmap<-raster("02-Outputs/Final Maps/MKD_clay_RF.tif")
maps <- stack(RKmap, RFmap)

# Explore the maps
names(maps)
summary(maps)
plot(maps)

# Load the validation dataset. 
# It was was prepared in the 'data_preparation_profiles' script
test <- read.csv("02-Outputs/dat_test_clay.csv")

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
test$PE_RK <- test$MKD_clay_RK - test$clay
test$PE_RF <- test$MKD_clay_RF - test$clay


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
AVE_RK <- 1 - sum(test$PE_RK^2) / sum( (test$clay - mean(test$clay))^2 )


# Random Forest

# Mean Error
ME_RF <- mean(test$PE_RF)
# Mean Absolute Error (MAE)
MAE_RF <- mean(abs(test$PE_RF))
# Root Mean Squared Error (RMSE)
RMSE_RF <- sqrt(sum(test$PE_RF^2) / length(test$PE_RF))
# Amount of Variance Explained (AVE)
AVE_RF <- 1 - sum(test$PE_RF^2) / sum( (test$clay - mean(test$clay))^2 )

table <-  data.frame(Quality_measure = c("AVE", "MAE", "ME","RMSE"))
table$RKmap <- c(AVE_RK, MAE_RK, ME_RK,  RMSE_RK)
table$RFmap <- c(AVE_RF, MAE_RF, ME_RF,  RMSE_RF)
table


