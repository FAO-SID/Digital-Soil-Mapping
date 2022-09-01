# Introduction to R

# 0. Playground

# learn important keyboard shortcuts
# Ctrl + enter for running code
# tab after writing the first three characters of the function name
# F1 to access the help

# explore the use of <-, $, [], ==, !=, c(), :, data.frame(), list(), as.factor()

a <- 10:15
a[2]
a[2:3]
b <- c("1", "a", a )
length(b)
df <- data.frame(column_a = 1:8, column_b = b)

df[,1]
df$column_b
as.numeric(df$column_b)
plot(df)

df[1:3,]
df[,1]

as.factor(b)

d <- list(a, b, df)
d
names(d)
names(d) <- c("numeric_vector", "character_vector", "dataframe")
d
d[[1]]
d$numeric_vector

a == b
a != b

# 1. Set working directory using the function 
setwd("C:/GIT/Digital-Soil-Mapping/")

# 2. Install and load readxl, tidyverse, and data.table packages using the functions
install.packages("tidyverse")
install.packages("readxl")
install.packages("data.table")
library(tidyverse)
library(readxl)
library(data.table)

# 3. Import an spreadsheet
# 3.1 Read the soil_data.xlsx file, spreadsheet 2, using the function 
read_excel(path = "01-Data/soil_data.xlsx", sheet = 2)

# 3.2 Read the csv file with the native function 
read.csv("01-Data/horizon.csv")

# 3.3 Read the csv file with the tidyverse function 
read_csv("01-Data/horizon.csv")

# 3.4 Read the csv file with the data.table function 
fread("01-Data/horizon.csv")

# 3.5 Assign the dataframe to an object called dat
dat <- read_csv("01-Data/horizon.csv")

# 4. Tidyverse functions =======================================================
## 4.1 Select pid, hip, top, bottom, ph_h2o, clay from dat ---------------------
dat_1 <- dat %>% 
  select(pid, hid, top, bottom, ph_h2o, clay)

## 4.2 Filter: pick observations by their values -------------------------------
# filter observations with clay > 50%

dat_2 <- dat_1 %>% 
  filter(clay > 50)

dat_2
## 4.3 Mutate: create a new variable -------------------------------------------
# thickness = top - bottom

dat_3 <- dat_2 %>% 
  mutate(thickness = bottom - top)

## 4.4 Group_by and summarise --------------------------------------------------
# group by variable pid
# summarise taking the mean of pH and clay

dat_4 <- dat_3 %>% 
  group_by(pid) %>% 
  summarise(mean_ph = mean(ph_h2o),
            mean_clay = mean(clay))

## 4.5 Reshape the table using pivot_longer ------------------------------------
# use dat_3
# put the names of the variables ph_h2o, clay and thickness in the column 
# variable and keep the rest of the table. Save in dat_5 

dat_5 <- dat_3 %>% 
  pivot_longer(ph_h2o:thickness, names_to = "soil_property", values_to = "value")

## 4.6 Join the table sites.csv with dat_3 -------------------------------------
# Load site.csv (in 01-Data folder) 
# Join its columns with dat_3 keeping all the rows of dat_3
# save the result as dat_6

sites <- read_csv("01-Data/site.csv") 
dat_6 <- dat_3 %>% 
  left_join(sites)
# or
dat_6 <- sites %>% 
  right_join(dat_3)

# 5. Data visualization with ggplot2 ===========================================
## 5.1 1D plot: histograms -----------------------------------------------------
# histograms of clay and ph_h2o

ggplot(dat_3, aes(x=clay)) + geom_histogram()

## 5.2 2D plot: scatterplot ----------------------------------------------------
# Scatterplot bottom vs. ph_h2o

ggplot(dat_3, aes(x = bottom, y = ph_h2o)) + 
  geom_point() 

# add a fitting line

ggplot(dat_3, aes(x = bottom, y = ph_h2o)) + 
  geom_point() +
  geom_smooth(method = "lm" )

## 5.3 3D plot: scatterplot ----------------------------------------------------
# Scatterplot bottom vs. ph_h2o, add clay as color and size inside the 
# function aes()

ggplot(dat_3, aes(x = bottom, y = ph_h2o, color = clay, size = clay)) + 
  geom_point() 

## 5.4 2D plot + facets --------------------------------------------------------
# Make a histogram of pH for each year
# use dat_6

ggplot(dat_6, aes(x = ph_h2o, fill = year)) +
  geom_histogram()+
  facet_wrap(~year)

# 6. Geospatial data with terra ================================================
## Load packages (install them if needed)
library(terra)

