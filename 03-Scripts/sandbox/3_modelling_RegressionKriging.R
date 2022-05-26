#Empty environment 
rm(list=ls())

# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

# Load the covariates stack. 
load(file = "02-Outputs/covariates.RData")
names(covs)
covs <-dropLayer(covs, "LandCover")

# Load the processed data for digital soil mapping. 
dat <- read.csv("02-Outputs/dat_train.csv")
names(dat)

library(sp)

class(dat)

# Promote to spatialPointsDataFrame and set the coordinate system
coordinates(dat) <- ~ X + Y
proj4string(dat) = CRS("+init=epsg:4326") # WGS84
class(dat)

# Check if the points overlay with covariates
plot(covs$B04CHE3)
points(dat)

# extract values from covariates to the soil points
dat <- extract(x = covs, y = dat, sp = TRUE)
summary(dat)


# Remove NA values
dat<-as.data.frame(dat)
dat <- dat[complete.cases(dat),]


str(dat)
#dat$LandCover <- as.factor(dat$LandCover)
dat$soilmap <- as.factor(dat$soilmap)

# prepare the table for regression (only OCSlog and the covariates)
datdf <- dat[, c("OCSlog", names(covs))]
# Fit a multiple linear regression model
model.MLR <- lm(OCSlog ~ ., data = datdf)
summary(model.MLR)


# stepwise variable selection
model.MLR.step <- step(model.MLR, direction="both")

# graphical diagnosis of the regression analysis
par(mfrow=c(2,2))
plot(model.MLR.step)
par(mfrow=c(1,1))

# collinearity test using variance inflation factors
library(car)
vif(model.MLR.step)
# problematic covariates should have 
# sqrt(VIF) > 2
sqrt(vif(model.MLR.step))
# Removing a layer from the stepwise model
model.MLR.step <- update(model.MLR.step, . ~ . - PRSCHE3)

# Project point data
coordinates(dat) <- ~ X + Y
proj4string(dat) = CRS("+init=epsg:4326") # WGS84
dat <- spTransform(dat, CRS("+init=epsg:6204")) # Macedonian state #coordinate system
covs <- projectRaster(covs, crs = CRS("+init=epsg:6204"),
                      method='ngb')
# Promote covariates to spatial grid dataframe (it can take a lot #of memory)
covs.sp <- as(covs, "SpatialGridDataFrame")
covs.sp@grid

# LandCover and soilmap are categorical variables, they need to be #'factor' type
#covs.sp$LandCover <-as.factor(covs.sp$LandCover) 
covs.sp$soilmap <-as.factor(covs.sp$soilmap)


library(gstat)
# Define gstat object 
gstat_data <- gstat(formula = as.formula(model.MLR.step$call$formula), data = dat)

#Compute and explore an experimental #semivariogram
vario <- variogram(gstat_data, cutoff=20000, width=500)
plot(vario, plot.nu=FALSE)


# Define initial semivariogram model
vario_mod <- vgm(nugget = 0.1, psill =
                   0.20, range = 5000, model = "Ste")
plot(vario, vario_mod)


# Define the random numbers table 
#(to get reproducible result)
#set.seed(12042019)
# Fit semivariogram model for kriging
vario_mod <- fit.variogram(vario, 
                           vario_mod, fit.method=7)
plot(vario, vario_mod)
vario_mod

# Make a prediction across all Macedonia using Regression Kriging #model
pred_gstat <- krige(formula = as.formula(model.MLR.step$call$formula), 
                    locations = dat, 
                    newdata = covs.sp, 
                    model = vario_mod, 
                    debug.level = -1)

# Back transform predictions log transformed
RKpred <- exp(raster(pred_gstat))

# Back transform the coordinate system to WGS 84, using MLR #prediction as a template
RKpred <- projectRaster(from = RKpred, to = covs$DEMENV5, method = "ngb")

# Explore and save the result as a tiff file
plot(RKpred)
writeRaster(RKpred,'02-Outputs/Final Maps/MKD_OCS_RK.tif', overwrite=TRUE)

# Make an uncertainty estimation as a map of standard deviations
# Standard deviation is the square root of kriging variance
RKsd_log <- sqrt(raster(pred_gstat, layer='var1.var'))
# Calculate 95% confidence interval on 
# log-transformed data as pred +-2*sd
RK_ci_high<-raster(pred_gstat)+2*RKsd_log
RK_ci_low<-raster(pred_gstat)-2*RKsd_log
# Confidence interval limits can be back transformed to t/ha
RK_ci_high<-exp(RK_ci_high)
RK_ci_low<-exp(RK_ci_low)
plot(RK_ci_high)
plot(RK_ci_low)
# Export limits of confidence interval as measures of uncertainty
writeRaster(RK_ci_high,'02-Outputs/Final Maps/MKD_OCS_RK_ci95_high.tif', overwrite=TRUE)
writeRaster(RK_ci_low,'02-Outputs/Final Maps/MKD_OCS_RK_ci95_low.tif', overwrite=TRUE)
