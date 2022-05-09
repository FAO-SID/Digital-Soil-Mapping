# This script is a part of the FAO/GSP training course on Digital Soil Mapping.
# It is used to prepare soil profile data for digital soil mapping of SOC stocks.
# It can be modified to be used for mapping other soil properties.
#################################################################################################################
# Install all required R packages used in the cookbook
# install.packages(c("aqp", "automap", "car", "caret", "e1071",
#                    "mpspline2", "htmlwidgets", "leaflet", "mapview",
#                    "Metrics", "openair", "plotKML", "psych", 
#                    "quantregForest", "randomForest", "raster",
#                    "rasterVis", "reshape", "rgdal", "RSQLite",
#                    "snow", "soiltexture", "sf", "sp"))
# Remove existing objects in memory and clear memory cache
rm(list = ls());
gc()

library(mpspline2)
# Set working directory
setwd("C:/Users/hp/Documents/FAO/EduSoils/training_material")

#####################################  Import, explore and clean the data ########################################

# Import soil layers (horizons) data from a .csv table
dat_layers <- read.csv(file = "01-Data/horizons.csv")

# Explore the data
str(dat_layers)
summary(dat_layers)

# Import site-level data from a .csv table
dat_sites <- read.csv(file = "01-Data/site-level.csv")

# Explore the data
str(dat_sites)
summary(dat_sites)

# Remove duplicate profiles
test <- duplicated(dat_sites[c(1,3,4)]) # are there duplicates present, based on similar coordinates and profile id?
summary(test)
dat_sites <-  dat_sites[!duplicated(dat_sites[c(1,3,4)]),]

# Merge site-level data with soil layers (horizons)  data
dat <- merge (x=dat_sites, y=dat_layers, by="ProfID")
summary(dat)
names (dat)

# Simplify the table by choosing only the parameters that are needed to model SOC stocks: 
# profile id, coordinates, depth of each layer, SOC, bulk density and coarse fragments.
# Note that SOC is OM/1.724 (OM - organic matter)
dat <- data.frame(id=dat$ProfID,
                  X=dat$X,
                  Y=dat$Y,
                  soil=dat$soiltype,
                  top=dat$DepthFrom,
                  bottom=dat$DepthTo,
                  clay=dat$Clay_perc)

# Explore and clean the SOC data
summary(dat$clay)
dat <- dat[!is.na(dat[,"clay"]),] # remove NA values


# Explore SOC data, identify outliers
summary(dat$clay)
hist(dat$clay, breaks=100)
boxplot(dat$clay, horizontal=TRUE)

# Outliers should be carefully explored and compared with literature values.
# Only if it is clear that outliers represent impossible or highly unlikely values,
# they should be removed as errors.
dat[dat$clay>70, 'soil']
dat[dat$clay>70, c('soil','top' ,'bottom' ,'clay')]


# Carbon content higher than 15% is only typical for organic soil (histosols)
# We will remove all non-histosols with atypically high SOC as outliers




################################ Calculating SOC stocks for standard topsoil depth 0-30cm #######################

dat <-
  dat[, c("id","top","bottom","clay","X",  "Y", "soil")]


#equal area spline clay
clay <-mpspline(dat, 'clay', d = c(0,30))
temp <- data.frame()
for (i in names(clay)){
  
  t2 <- cbind(clay[[i]]$id,clay[[i]]$est_dcm[1])
  
  
  temp <- rbind(temp,t2)
  
}
clay <- temp



##Since we are predicting clay we don't need Profiles without clay data
# keed only unique ids
dat <-dat[,c("id","X",  "Y", "soil") ]
dat <-dat[!duplicated(dat$id),]
#Rename columns of mpspline output
colnames(clay) <- c('id', 'clay')


#Merge the harmonized clay output with the data
dat <- merge(clay, dat[,c("id","X",  "Y", "soil") ], by='id', all.x=T)
# We are mapping clay therefore we exclude rows with NA for the clay column
dat <- dat[complete.cases(dat$clay),]


dat$clay <- as.numeric(dat$clay)


# Save the final table in a .csv file
write.csv(dat, "02-Outputs/dataproc_clay.csv", row.names = FALSE)

############ Splitting the dataset in calibraition (training the model) and validation (testing the model) ##############

library(caret)

# Define the random numbers table (to get reproducible result)
set.seed(42)

# Create random selection of 75% of the data as 'train' dataset and 25% as 'test' dataset
train.ind <- createDataPartition(1:nrow(dat), p = .75, list = FALSE)
train <- dat[ train.ind,]
test  <- dat[-train.ind,]

# Check if both 'train' and 'test datasets' have similar distributions
summary(train$clay)
summary(test$clay)

plot(density (train$clay), col='red',
     main='Statistical distribution of train and test datasets')
lines(density(test$clay), col='blue')
legend('topright', legend=c("train", "test"),
       col=c("red", "blue"), lty=1, cex=1.5)

# Save the 'train' and 'test' datasets as .csv tables

write.csv(train, file="02-Outputs/dat_train_clay.csv", row.names = FALSE)
write.csv(test, file="02-Outputs/dat_test_clay.csv", row.names = FALSE)

