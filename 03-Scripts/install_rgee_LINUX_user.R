
# Script: Script with installation procedures to use the RGEE ------------------------------------------------------------
# Author: Ricardo Dal'Agnol da Silva (ricds@hotmail.com)
# Date Created: 2021-10-26
# R version 4.1.1 (2021-08-10)
#

# clean environment
rm(list = ls()); gc()

# general libraries
#install.packages("pacman")
library(pacman)



# GEE account -------------------------------------------------------------

## you need a GEE account
## log in the https://code.earthengine.google.com/ and register for one

# LINUX USERS - install and setup rgee --------------------------------------------------
# Install remotes package if it is not installed yet
# install.packages("remotes")
# install rgee
remotes::install_github("r-spatial/rgee")

# Additionally, rgee depends on the Python packages: 
# numpy and ee. To install them, users can follow any of these three methods:
# 1.Use ee_install (Highly recommended for users with 
# no experience with Python environments)
# https://github.com/r-spatial/rgee#installation

# follow the instructions
rgee::ee_clean_pyenv()
rgee::ee_install()
library(rgee)
ee_check()

# if the current version is not updated, run the following line
rgee::ee_install_upgrade()
library(rgee)
ee_check()
# rgee_environment_dir = '~/.virtualenvs/rgee/'
# the result must be:
# ◉  Python version
# ✔ [Ok] ~/.virtualenvs/rgee/bin/python v3.8
# ◉  Python packages:
# ✔ [Ok] numpy
# ✔ [Ok] earthengine-api

# pre-requirements for R --------------------------------------------------
## R: version at least 3.6 (this is the version that I tested so far and works)
# Link: https://cran.r-project.org/bin/windows/base/

## RStudio: a recent version is recommended.
## Older versions do not show the GEE images in Viewer correctly.
# Link: https://www.rstudio.com/products/rstudio/download/

## RTools: needed to build some packages in R from source
# Link: https://cran.r-project.org/bin/windows/Rtools/
## after installing the Rtools, make sure to run this command line below inside RStudio:
# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")



# R packages -------------------------------------------------------------

## if you installed everything above, you can now install the packages inside R

# install/load general packages used in the scripts
p_load(raster,
       rgdal,
       rgeos,
       gdalUtils,
       sp,
       sf,
       leaflet,
       mapview,
       caret)

# if gdalUtils does not install correctly, use the following line
# devtools:::install_github("gearslaboratory/gdalUtils")

# now some more specific packages related to using the rgee
p_load(rgee, geojsonio, remotes, reticulate, devtools, googledrive)

## sometimes at this point you are required to restart R or the computer before proceeding
## try restarting if the installation do not finish properly and run the installation again after restart

# Initialize the Python Environment
# to clean credentials: ee_clean_credentials()
rgee::ee_Initialize(drive = T)
# The result should be:
# ── rgee 1.1.3 ──────────────────────────────────── earthengine-api 0.1.302 ── 
# ✔ user: not_defined
# ✔ Google Drive credentials:  FOUND
# ✔ Initializing Google Earth Engine:  DONE!
# ✔ Earth Engine account: users/<YOUR_USERNAME>


