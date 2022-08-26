#
# Digital Soil Mapping
# Soil Profile Data
# Cleaning and Processing
#
# GSP-Secretariat
# Contact: Isabel.Luotto@fao.org
#          Marcos.Angelini@fao.org
#_______________________________________________________________________________


rm(list = ls())
gc()

# Content of this script =======================================================
# The goal of this script is to organise the soil data for mapping, including:
# 
# User-defined variables 
# 0 - Set working directory and load necessary packages
# 1 - Import national data 
# 2 - select useful columns
# 3 - Quality check
# 4 - Estimate BD using pedotransfer function
# 5 - Estimate Organic Carbon Stock (ocs) 
# 6 - Harmonize soil layers
# 7 - Plot and save results
#_______________________________________________________________________________

# User-defined variables =======================================================
wd <- 'C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping'
#wd <- "C:/GIT/Digital-Soil-Mapping"

# 0 - Set working directory and load necessary packages ========================
setwd(wd) # change the path accordingly

library(tidyverse) # for data management and reshaping
library(readxl) # for importing excel files
library(mapview) # for seing the profiles in a map
library(sf) # to manage spatial data (shp vectors) 
library(aqp) # for soil profile data
# require("devtools") 
# devtools::install_bitbucket("brendo1001/ithir/pkg") #install ithir package
library(ithir) # for horizon harmonization


# 1 - Import national data =====================================================
# Save your national soil dataset in the data folder /01-Data as a .csv file or 
# as a .xlsx file

## 1.1 - for .xlsx files -------------------------------------------------------
# Import horizon data 
hor <- read_excel("01-Data/soil_data.xlsx", sheet = 2)
# Import site-level data
site <- read_excel("01-Data/soil_data.xlsx", sheet = 1)

## 1.2 - for .csv files --------------------------------------------------------
# Import horizon data 
hor <- read_csv(file = "01-Data/horizon.csv")
# Import site-level data
site <- read_csv(file = "01-Data/site.csv")

# change names of key columns
names(site)
names(site)[1] <- "ProfID"
names(hor)
names(hor)[1] <- "ProfID"
names(hor)[2] <- "HorID"
# scan the data
summary(site)
summary(hor)

# 2 - select useful columns ====================================================
## 2.1 - select columns ---------------------------------------------
hor <- select(hor, ProfID, HorID, top, bottom, humus, bd, 
              crf=coarse_frag, clay, ph=ph_h2o)
hor$soc <- hor$humus/1.724

# remove humus
hor <- select(hor, -humus)

# assume that Coarse fragments = NA means no coarse fragments present
hor <- mutate(hor, 
              crf = ifelse(is.na(crf), 0, crf))


# 3 - Quality check ============================================================

## 3.1 - Check locations -------------------------------------------------------
# https://epsg.io/6204
site %>% 
  st_as_sf(coords = c("x", "y"), crs = 6204) %>% # convert to spatial object
  mapview(zcol = "year", cex = 2, lwd = 0.1) # visualise in an interactive map

# profile P5810 is wrongly located, so let's remove it
site <- filter(site, ProfID != "P5810")


## 3.2 - Convert data into a Soil Profile Collection ---------------------------
depths(hor) <- ProfID ~ top + bottom
site(hor) <- left_join(site(hor), site)
profiles <- hor

profiles

## 3.3 - plot first 20 profiles using pH as color ------------------------------
plotSPC(x = profiles[1:20], name = "clay", color = "clay",
        name.style = "center-center")

## 3.4 - check data integrity --------------------------------------------------
# A valid profile is TRUE if all of the following criteria are false:
#    + depthLogic : boolean, errors related to depth logic
#    + sameDepth : boolean, errors related to same top/bottom depths
#    + missingDepth : boolean, NA in top / bottom depths
#    + overlapOrGap : boolean, gaps or overlap in adjacent horizons
aqp::checkHzDepthLogic(profiles)

# visualize some of these profiles by the pid
subset(profiles, grepl("P0142", ProfID, ignore.case = TRUE))
subset(profiles, grepl("P0494", ProfID, ignore.case = TRUE))
subset(profiles, grepl("P3847", ProfID, ignore.case = TRUE))


## 3.5 - keep only valid profiles ----------------------------------------------
clean_prof <- HzDepthLogicSubset(profiles)
metadata(clean_prof)$removed.profiles
write_rds(clean_prof, "01-Data/soilProfileCollection.rds")

## 3.6 convert soilProfileCollection to a table --------------------------------
dat <- left_join(clean_prof@site, clean_prof@horizons)
dat <- dat <- select(dat, ProfID, HorID, x, y, year, top, bottom, bd:soc )

# 4 - Estimate BD Stock using pedotransfer functions ==========================

## 4.1 - Select a pedotransfer function ----------------------------------------
# create a vector of BD values to test the best fitting pedotransfer function
BD_test <- tibble(SOC = clean_prof@horizons$soc, 
                  BD_test = clean_prof@horizons$bd)
BD_test <-  na.omit(BD_test)

# create the function with all PTF
estimateBD <- function(SOC=NULL, method=NULL){
  OM <- SOC * 1.724
  if(method=="Saini1996"){BD <- 1.62 - 0.06 * OM}
  if(method=="Drew1973"){BD <- 1 / (0.6268 + 0.0361 * OM)}
  if(method=="Jeffrey1979"){BD <- 1.482 - 0.6786 * (log(OM))}
  if(method=="Grigal1989"){BD <- 0.669 + 0.941 * exp(1)^(-0.06 * OM)}
  if(method=="Adams1973"){BD <- 100 / (OM /0.244 + (100 - OM)/2.65)}
  if(method=="Honeyset_Ratkowsky1989"){BD <- 1/(0.564 + 0.0556 * OM)}
  return(BD)
}

## 4.2 - Estimate BLD for a subset using the pedotransfer functions ------------
BD_test$Saini <- estimateBD(BD_test$SOC, method="Saini1996")
BD_test$Drew <- estimateBD(BD_test$SOC, method="Drew1973")
BD_test$Jeffrey <- estimateBD(BD_test$SOC, method="Jeffrey1979")
BD_test$Grigal <- estimateBD(BD_test$SOC, method="Grigal1989")
BD_test$Adams <- estimateBD(BD_test$SOC, method="Adams1973")
BD_test$Honeyset_Ratkowsky <- estimateBD(BD_test$SOC, 
                                         method="Honeyset_Ratkowsky1989")

## 4.3 Compare results ---------------------------------------------------------

# Observed values:
summary(BD_test$BD_test)

# Predicted values:
summary(BD_test$Saini) 
summary(BD_test$Drew) 
summary(BD_test$Jeffrey) 
summary(BD_test$Grigal)
summary(BD_test$Adams)
summary(BD_test$Honeyset_Ratkowsky) 

# Compare data distributions for observed and predicted BLD
plot(density(BD_test$BD_test),type="l",col="black", ylim=c(0,5), 
     lwd=2, main="Bulk Density Pedotransfer Functions")
lines(density(BD_test$Saini),col="green", lwd=2)
lines(density(BD_test$Drew),col="red", lwd=2)
lines(density(BD_test$Jeffrey),col="cyan", lwd=2)
lines(density(BD_test$Grigal),col="orange", lwd=2)
lines(density(BD_test$Adams),col="magenta", lwd=2)
lines(density(BD_test$Honeyset_Ratkowsky),col="blue", lwd=2)
legend("topleft",
       legend = c("Original", "Saini", "Drew", "Jeffrey", "Grigal", "Adams",
                  "Honeyset_Ratkowsky"), 
       fill=c("black", "green", "red", "cyan", "orange","magenta", "blue"))

# Plot the Selected function again
plot(density(BD_test$BD_test),type="l",col="black", ylim=c(0,3.5),
     lwd=2, main="Bulk Density Selected Function")
lines(density(BD_test$Honeyset_Ratkowsky),col="blue", lwd=2)
legend("topleft",legend = c("Original", "Honeyset_Ratkowsky"), 
       fill=c("black", "blue"))


## 4.4 Estimate BD for the missing horizons ------------------------------------
dat$bd[is.na(dat$bd)] <- 
  estimateBD(dat[is.na(dat$bd),]$soc, method="Honeyset_Ratkowsky1989")

# Explore the results
summary(dat$bd)
plot(density(BD_test$BD_test),type="l",col="black", ylim=c(0,3.5), 
     lwd=2, main="Bulk Density Gap-Filling")
lines(density(dat$bd, na.rm = TRUE), col="green", lwd=2)
legend("topleft",legend = c("Original", "Original+Estimated"), 
       fill=c("black", "green"))


## 4.5 - Explore outliers ------------------------------------------------------
# Outliers should be carefully explored and compared with literature values.
# Only if it is clear that outliers represent impossible or highly unlikely 
# values, they should be removed as errors.
# 
# Carbon content higher than 15% is only typical for organic soil (histosols)
# We will remove all atypically high SOC as outliers
summary(dat$soc)
na.omit(dat$ProfID[dat$soc>15])
dat <- dat[dat$ProfID!= "P4601",]

# Explore bulk density data, identify outliers
# remove layers with Bulk Density < 1 g/cm^3
low_bd_profiles <- na.omit(dat$ProfID[dat$bd<1])
dat <- dat[!(dat$ProfID %in% low_bd_profiles),]

# Explore data, identify outliers
x <- pivot_longer(dat, cols = bd:soc, values_to = "value",
                  names_to = "soil_property")
x <- na.omit(x)
ggplot(x, aes(x = soil_property, y = value, fill = soil_property)) +
  geom_boxplot() + 
  facet_wrap(~soil_property, scales = "free")


# 5 - Harmonize soil layers ====================================================
## 5.1 - Set target soil properties and depths ---------------------------------
names(dat)
dat <- select(dat, ProfID, HorID, x, y, top, bottom, crf, bd, soc, clay, ph)

target <- c("crf", "bd", "soc", "clay", "ph")
depths <- t(c(0,30,60,100))

## 5.2 - Create standard layers ------------------------------------------------
d <- unique(select(dat, ProfID, x, y))

for (i in seq_along(target)) {
  vlow <- min(dat[,target[i]][[1]], na.rm = TRUE)
  vhigh <- max(dat[,target[i]][[1]], na.rm = TRUE)
  o <- dat[,c("ProfID", "top", "bottom",target[i])] %>% 
    na.omit() %>%
    as.data.frame(stringsAsFactors = FALSE)
  x <- ithir::ea_spline(obj = o, var.name = target[i], d = depths, 
                        vlow = vlow[[1]], vhigh = vhigh[[1]])$harmonised 
  x[x==-9999] <- NA
  x <- x %>% 
    as_tibble() %>% 
    select(-`soil depth`)
  
  names(x) <- c("ProfID",paste0(target[i],c("_0_30","_30_60","_60_100")))
  d <- d %>% left_join(x)
}
d

# 6 - Estimate Organic Carbon Stock (ocs) ======================================
## 6.1 - Estimate Organic Carbon Stock -----------------------------------------
# SOC must be in g/kg (% * 10)
# BLD in kg/m3
# CRF in percentage

# 0 - 30 cm
ORCDRC <- d$soc_0_30*10
HSIZE <- 30
BLD <- d$bd_0_30*1000
CRFVOL <- d$crf_0_30

OCSKG_0_30 <- ORCDRC/1000 * HSIZE/100 * BLD * (100 - CRFVOL)/100

# Convert Organic Carbon Stock from kg/m3 to t/ha
d$ocs_0_30 <- OCSKG_0_30*10

# 30 - 60 cm
ORCDRC <- d$soc_30_60*10
HSIZE <- 30
BLD <- d$bd_30_60*1000
CRFVOL <- d$crf_30_60

OCSKG_30_60 <- ORCDRC/1000 * HSIZE/100 * BLD * (100 - CRFVOL)/100

# Convert Organic Carbon Stock from kg/m3 to t/ha
d$ocs_30_60 <- OCSKG_30_60*10

# 60 - 100 cm
ORCDRC <- d$soc_60_100*10
HSIZE <- 40
BLD <- d$bd_60_100*1000
CRFVOL <- d$crf_60_100

OCSKG_60_100 <- ORCDRC/1000 * HSIZE/100 * BLD * (100 - CRFVOL)/100

# Convert Organic Carbon Stock from kg/m3 to t/ha
d$ocs_60_100 <- OCSKG_60_100*10


# 7 - Plot  and save results ===================================================

x <- pivot_longer(d, cols = crf_0_30:ocs_60_100, values_to = "value",
                  names_sep = "_", 
                  names_to = c("soil_property", "top", "bottom"))
x <- mutate(x, depth = paste(top, "-" , bottom))
x <- na.omit(x)
ggplot(x, aes(x = depth, y = value, fill = soil_property)) +
  geom_boxplot() + 
  facet_wrap(~soil_property, scales = "free")

# remove BD and CF
d <- select(d, ProfID:y, soc_0_30:ocs_60_100)
# save data
write_csv(d, "02-Outputs/soil_data.csv")
