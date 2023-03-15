# This script is a part of the FAO/GSP training course on Digital Soil Mapping.
# It is used to map SOC stocks with Regression Kriging algorithm.
# It can be modified to be used for mapping other soil properties.
#################################################################################################################

# Remove existing objects in memory and clear memory cache
rm(list = ls());
gc()

# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

############################### Prepare the final table for modelling (regression matrix) ########################
library (sp)
library(raster)

# Load the covariates stack. It was was prepared in the 'data_preparation_covariates' script
load(file = "02-Outputs/covariates_clay.RData")
names(covs)

# Load the processed data for digital soil mapping. This table was prepared in the 'data_preparation_profiles' script
dat <- read.csv("02-Outputs/dat_train_clay.csv")
names(dat)


class(dat)

# Promote to spatialPointsDataFrame and set the coordinate system
coordinates(dat) <- ~ X + Y
proj4string(dat) = CRS("+init=epsg:4326") # WGS84
class(dat)


# extract values from covariates to the points
dat <- extract(x = covs, y = dat, sp = TRUE)
summary(dat)

# Remove NA values
dat<-as.data.frame(dat)
dat <- dat[complete.cases(dat),]

str(dat)

# LandCover and soilmap are categorical variables, 
# they need to be 'factor' type
dat$LandCover <- as.factor(dat$LandCover)
dat$soilmap <- as.factor(dat$soilmap)
str(dat)

# save the final table (regression matrix)
write.csv(dat, "02-Outputs/SOC_RegMatrix_clay.csv", row.names = FALSE)

############################### Multiple Linear Regression ##################################################

# prepare the table for regression, including only clay and the covariates
datdf <- dat[, c("clay", names(covs))]

# Fit a multiple linear regression model between clay and the covariates
model.MLR <- lm(clay ~ ., data = datdf)
summary(model.MLR)

# stepwise variable selection
model.MLR.step <- step(model.MLR, direction="both")

# summary and anova of the new model using stepwise covariates
# selection
summary(model.MLR.step)
anova(model.MLR.step)

# graphical diagnosis of the regression analysis
par(mfrow=c(2,2))
plot(model.MLR.step)
par(mfrow=c(1,1))

# collinearity test using variance inflation factors
library(car)
vif(model.MLR.step)

# problematic covariates should have sqrt(VIF) > 2
sqrt(vif(model.MLR.step))

# Removing a layer from the stepwise model
 model.MLR.step <- update(model.MLR.step, . ~ . - B07CHE3 )
 model.MLR.step <- update(model.MLR.step, . ~ . - B04CHE3 )
## summary  of the new model using stepwise covariates selection
summary(model.MLR.step)


################################################ Regression Kriging #################################################

# To use kriging, all data should be converted to a projected coordinate systems (to meters)

# Project point data
coordinates(dat) <- ~ X + Y
proj4string(dat) = CRS("+init=epsg:4326") # WGS84
dat <- spTransform(dat, CRS("+init=epsg:6204")) # Macedonian state coordinate system

# Project covariates

covs <- projectRaster(covs, crs = CRS("+init=epsg:6204"),
                      method='ngb')

# Promote covariates to spatial grid dataframe (it can take a lot of memory)
covs.sp <- as(covs, "SpatialGridDataFrame")

# LandCover and soilmap are categorical variables, they need to be 'factor' type
covs.sp$LandCover <-as.factor(covs.sp$LandCover) 
covs.sp$soilmap <-as.factor(covs.sp$soilmap) 

#  Define gstat object 
library(gstat)
gstat_data <- gstat(formula = as.formula(model.MLR.step$call$formula), data = dat)

# Compute and explore an experimental semivariogram
vario <- variogram(gstat_data, cutoff=20000, width=500)
plot(vario, plot.nu=FALSE)

# Define initial semivariogram model
vario_mod <- vgm(nugget = 0.1, psill = 0.20, range = 5000, model = "Ste")
plot(vario, vario_mod)

# Define the random numbers table (to get reproducable result)
set.seed(12042019)

# Fit semivariogram model for kriging
vario_mod <- fit.variogram(vario, vario_mod, fit.method=7)
plot(vario, vario_mod)
vario_mod

# Make a prediction across all Macedonia using Regresion Kriging model
pred_gstat <- krige(formula = as.formula(model.MLR.step$call$formula), 
                    locations = dat, 
                    newdata = covs.sp, 
                    model = vario_mod, 
                    debug.level = -1)


RKpred <- raster(pred_gstat)

# Back tranform the coordinate system to WGS 84, using MLR prediction as a template
RKpred <- projectRaster(from = RKpred, to = covs$DEMENV5, method = "ngb")

# Explore and save the result as a tiff file
plot(RKpred)
writeRaster(RKpred,'02-Outputs/Final Maps/MKD_clay_RK.tif', overwrite=TRUE)

# Make an uncertainty estimation as a map of standard deviations
# Standard deviation is the square root of kriging variance
RKsd <- sqrt(raster(pred_gstat, layer='var1.var'))
# If data was not log-transformed, than standard deviation map is in same units as the data (t/ha) and can be reported as uncertainty map
plot(RKsd, col=  topo.colors(255))
# But in our case it was calculated on log-transformed data and cannot be back transformed to t/ha
writeRaster(RKsd,'02-Outputs/Final Maps/MKD_clay_RK_sd.tif', overwrite=TRUE)

