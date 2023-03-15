# This script is a part of the FAO/GSP training course on Digital Soil Mapping.
# It is used to map SOC stocks with Random Forest algorithm.
# It can be modified to be used for mapping other soil properties.
#################################################################################################################

# Remove existing objects in memory and clear memory cache
rm(list = ls());
gc()

# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

############################### Prapare the final table for modelling (regression matrix) ########################

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
dat$LandCover <- as.factor(dat$LandCover)
dat$soilmap <- as.factor(dat$soilmap)
str(dat)

# Cave the final table and all the covariates
write.csv(dat, "02-Outputs/SOC_RegMatrix_clay.csv", row.names = FALSE)

############################### Spatial modelling with Random Forest ##############################################

# Define the random numbers table (to get reproducable result)
set.seed(12042019)

library(sp)

# Promote to spatialPointsDataFrame and set the coordinate system
coordinates(dat) <- ~ X + Y
proj4string(dat) = CRS("+init=epsg:4326") # WGS84

names(dat)

# We need to define a formula for the model
fm = as.formula(paste("clay ~", paste0(names(covs),
                                       collapse = "+")))
fm

# Run the Random Forest model and explore the results
library(randomForest)

rfmodel <- randomForest(fm, data=dat, ntree=500, importance=TRUE) 
rfmodel

# Explore the importance of covariates in the model
varImpPlot(rfmodel)


# Make a prediction across all Macedonia
pred <- predict(covs, rfmodel)


# Explore and save the result as a tiff file
plot(pred)
writeRaster(pred, filename = "02-Outputs/Final Maps/MKD_clay_RF.tif",
            overwrite=TRUE)

