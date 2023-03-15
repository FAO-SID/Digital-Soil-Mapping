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
                  SOC=dat$OM_perc/1.724,
                  BLD=dat$Bulk_density_gcm3,
                  CRF=dat$Stones_perc)

# Explore and clean the SOC data
summary(dat$SOC)
dat <- dat[!is.na(dat[,"SOC"]),] # remove NA values
dat <- dat[dat$SOC>0,] # remove 0 values

# Explore SOC data, identify outliers
summary(dat$SOC)
hist(dat$SOC, breaks=100)
boxplot(dat$SOC, horizontal=TRUE)

# Outliers should be carefully explored and compared with literature values.
# Only if it is clear that outliers represent impossible or highly unlikely values,
# they should be removed as errors.
summary(dat[dat$SOC>15,])

# Carbon content higher than 15% is only typical for organic soil (histosols)
# We will remove all non-histosols with atypically high SOC as outliers
dat <- dat[ !(dat$SOC>15 & dat$soil!='Histosol'),] 

# Explore bulk density data, identify outliers
summary(dat$BLD)
dat$BLD[dat$BLD==0]<-NA   # bulk density cannot be 0, it must be NA
hist(dat$BLD, breaks=50)  
boxplot(dat$BLD, horizontal=TRUE)$out
# BLD values higher than 2.5 g?cm are not typical for fine earth,
#   and most likely correspond to coarse fragments; they should be removed
dat <- dat[!(dat$BLD>=2.5 & !is.na(dat$BLD)),]

# BLD values lower than 1 g?cm3 can only correspond to organic soils (Histosols)
# We should check if all low BLD values correspond to Histosols
summary(dat[dat$BLD<1,])

# Explore and clean coarse fragments data
summary(dat$CRF)
hist(dat$CRF, breaks=50)
boxplot(dat$CRF)

# Remove outliers automatically, using boxplot
out <- boxplot(dat$CRF, horizontal=TRUE)$out
dat <- dat[!(dat$CRF %in% out),]

# The points without CRF, should have 0 values
dat$CRF[is.na(dat$CRF)]<-0


###################################### Calculate BLD using pedotransfer functions ###############################

# Select the best fitting pedotransfer function using subset of data that has measured BLD

BD_test <- dat[is.na(dat$BLD)==FALSE,] 
estimateBD <- function(SOC, method="Saini1996"){
  OM <- SOC * 1.724
  if(method=="Saini1996"){BD <- 1.62 - 0.06 * OM}
  if(method=="Drew1973"){BD <- 1 / (0.6268 + 0.0361 * OM)}
  if(method=="Jeffrey1979"){BD <- 1.482 - 0.6786 * (log(OM))}
  if(method=="Grigal1989"){BD <- 0.669 + 0.941 * exp(1)^(-0.06 * OM)}
  if(method=="Adams1973"){BD <- 100 / (OM /0.244 + (100 - OM)/2.65)}
  if(method=="Honeyset_Ratkowsky1989"){BD <- 1/(0.564 + 0.0556 * OM)}
  return(BD)
}

# Estimate BLD for a subset using the pedotransfer functions

BD_test$Saini <- estimateBD(BD_test$SOC, method="Saini1996")
BD_test$Drew <- estimateBD(BD_test$SOC, method="Drew1973")
BD_test$Jeffrey <- estimateBD(BD_test$SOC, method="Jeffrey1979")
BD_test$Grigal <- estimateBD(BD_test$SOC, method="Grigal1989")
BD_test$Adams <- estimateBD(BD_test$SOC, method="Adams1973")
BD_test$Honeyset_Ratkowsky <- estimateBD(BD_test$SOC, method="Honeyset_Ratkowsky1989")

# Compare results

# Observed values:
summary(BD_test$BLD)

# Predicted values:
summary(BD_test$Saini) 
summary(BD_test$Drew) 
summary(BD_test$Jeffrey) 
summary(BD_test$Grigal)
summary(BD_test$Adams)
summary(BD_test$Honeyset_Ratkowsky) 

# Compare data distributions for observed and predicted BLD
plot(density(BD_test$BLD),type="l",col="black", ylim=c(0,5), lwd=2, main="Bulk Density Pedotransfer Functions")
lines(density(BD_test$Saini),col="green", lwd=2)
lines(density(BD_test$Drew),col="red", lwd=2)
lines(density(BD_test$Jeffrey),col="cyan", lwd=2)
lines(density(BD_test$Grigal),col="orange", lwd=2)
lines(density(BD_test$Adams),col="magenta", lwd=2)
lines(density(BD_test$Honeyset_Ratkowsky),col="blue", lwd=2)
legend("topleft",legend = c("Original", "Saini", "Drew", "Jeffrey", "Grigal", "Adams",
                            "Honeyset_Ratkowsky"), fill=c("black", "green", "red", "cyan", 
                                                          "orange","magenta", "blue"))

# Plot the Selected function again
plot(density(BD_test$BLD),type="l",col="black", ylim=c(0,3.5), lwd=2, main="Bulk Density Selected Function")
lines(density(BD_test$Honeyset_Ratkowsky),col="blue", lwd=2)
legend("topleft",legend = c("Original", "Honeyset_Ratkowsky"), fill=c("black", "blue"))

# Estimate BLD for the missing points with the selected pedotransfer function
dat$BLD[is.na(dat$BLD)] <- estimateBD(dat[is.na(dat$BLD),]$SOC, method="Honeyset_Ratkowsky1989")

# Explore the results
summary(dat$BLD)
plot(density(BD_test$BLD),type="l",col="black", ylim=c(0,3.5), lwd=2, main="Bulk Density Gap-Filling")
lines(density(dat$BLD),col="green", lwd=2)
legend("topleft",legend = c("Original", "Original+Estimated"), fill=c("black", "green"))

################################ Calculating SOC stocks for standard topsoil depth 0-30cm #######################

dat <-
  dat[, c("id","top","bottom","SOC","BLD","CRF","X",  "Y", "soil")]

#equal area spline SOC
SOC <-mpspline(dat, 'SOC', d = c(0,30))
temp <- data.frame()
for (i in names(SOC)){
  
  t2 <- cbind(SOC[[i]]$id,SOC[[i]]$est_dcm[1])
  
  
  temp <- rbind(temp,t2)
  
}
SOC <- temp

#equal area spline BLD
BLD <-mpspline(dat, 'BLD', d = c(0,30))
temp <- data.frame()
for (i in names(BLD)){
  
  t2 <- cbind(BLD[[i]]$id,BLD[[i]]$est_dcm[1])
  
  
  temp <- rbind(temp,t2)
  
}
BLD <- temp

#equal area spline CRF
CRFVOL <-mpspline(dat, 'CRF', d = c(0,30))
temp <- data.frame()
for (i in names(CRFVOL)){
  
  t2 <- cbind(CRFVOL[[i]]$id,CRFVOL[[i]]$est_dcm[1])
  
  
  temp <- rbind(temp,t2)
  
}
CRFVOL <- temp


##Since we are predicting SOC we don't need Profiles without SOC data
# keed only unique ids
dat <-dat[,c("id","X",  "Y", "soil") ]
dat <-dat[!duplicated(dat$id),]
#Rename columns of mpspline output
colnames(SOC) <- c('id', 'SOC')
colnames(BLD) <- c('id', 'BLD')
colnames(CRFVOL) <- c('id', 'CRF')

#Merge the harmonized SOC output with the data
dat <- merge(SOC, dat[,c("id","X",  "Y", "soil") ], by='id', all.x=T)
# We are mapping SOC therefore we exclude rows with NA for the SOC column
dat <- dat[complete.cases(dat$SOC),]
dat <- dat[dat$SOC>0,]

#Finally merge the other attributes and set to numeric
dat <- merge(dat, BLD, by='id', all.x=T )
dat <- merge(dat, CRFVOL, by='id', all.x=T )

dat$SOC <- as.numeric(dat$SOC)
dat$BLD <- as.numeric(dat$BLD)
dat$CRF <- as.numeric(dat$CRF)
summary(dat)




# Estimate Organic Carbon Stock
# SOC must be in g/kg (% * 10)
# BLD in kg/m3
# CRF in percentage


#Create OCSKGM function
OCSKGM <-function (ORCDRC, BLD = 1400, CRFVOL = 0, HSIZE, ORCDRC.sd = 10, 
                   BLD.sd = 100, CRFVOL.sd = 5, se.prop = TRUE) 
{
  if (any(ORCDRC[!is.na(ORCDRC)] < 0) | any(BLD[!is.na(BLD)] < 
                                            0) | any(CRFVOL[!is.na(CRFVOL)] < 0)) {
    warning("Negative values for 'ORCDRC', 'BLD', 'CRFVOL' found")
  }
  OCSKG <- ORCDRC/1000 * HSIZE/100 * BLD * (100 - CRFVOL)/100
  if (se.prop == TRUE) {
    if (any(ORCDRC.sd[!is.na(ORCDRC.sd)] < 0)) {
      ORCDRC.sd = ifelse(is.na(ORCDRC.sd) | ORCDRC.sd < 
                           0, 0, ORCDRC.sd)
      warning("Replacing negative values for 'ORCDRC.sd'")
    }
    if (any(BLD.sd[!is.na(BLD.sd)] < 0)) {
      BLD.sd = ifelse(is.na(BLD.sd) | BLD.sd < 0, 0, BLD.sd)
      warning("Replacing negative values for 'BLD.sd'")
    }
    if (any(CRFVOL.sd[!is.na(CRFVOL.sd)] < 0)) {
      CRFVOL.sd = ifelse(is.na(CRFVOL.sd) | CRFVOL.sd < 
                           0, 0, CRFVOL.sd)
      warning("Replacing negative values for 'CRFVOL.sd'")
    }
    OCSKG.sd <- 1e-07 * HSIZE * sqrt(BLD^2 * (100 - CRFVOL)^2 * 
                                       ORCDRC.sd^2 + ORCDRC^2 * (100 - CRFVOL)^2 * BLD.sd^2 + 
                                       ORCDRC^2 * BLD^2 * CRFVOL.sd^2)
    attr(OCSKG, "measurementError") <- signif(OCSKG.sd, 
                                              3)
    attr(OCSKG, "units") <- "kilograms per square-meter"
  }
  return(OCSKG)
}





OCSKGM <- OCSKGM(ORCDRC = dat$SOC*10, BLD = dat$BLD*1000, CRFVOL = dat$CRF,
                 HSIZE = 30)
# Convert Organic Carbon Stock from kg/m3 to t/ha
dat$OCS <- OCSKGM*10

# Explore calculated SOC stocks
summary(dat$OCS)
hist(dat$OCS, breaks = 50)

#Remove Nas 
dat <- dat[complete.cases(dat$OCS),]

# Check if log-transformation improves the data distribution
hist(log(dat$OCS), breaks = 50)
# Add a new column for log-transformes carbon stocks
dat$OCSlog <- log(dat$OCS)


# Save the final table in a .csv file
write.csv(dat, "02-Outputs/dataproc.csv", row.names = FALSE)

############ Splitting the dataset in calibraition (training the model) and validation (testing the model) ##############

library(caret)

# Define the random numbers table (to get reproducible result)
set.seed(11042019)

# Create random selection of 75% of the data as 'train' dataset and 25% as 'test' dataset
train.ind <- createDataPartition(1:nrow(dat), p = .75, list = FALSE)
train <- dat[ train.ind,]
test  <- dat[-train.ind,]

# Check if both 'train' and 'test datasets' have similar distributions
summary(train$OCS)
summary(test$OCS)

plot(density (train$OCSlog), col='red',
     main='Statistical distribution of train and test datasets')
lines(density(test$OCSlog), col='blue')
legend('topright', legend=c("train", "test"),
       col=c("red", "blue"), lty=1, cex=1.5)

# Save the 'train' and 'test' datasets as .csv tables

write.csv(train, file="02-Outputs/dat_train.csv", row.names = FALSE)
write.csv(test, file="02-Outputs/dat_test.csv", row.names = FALSE)

