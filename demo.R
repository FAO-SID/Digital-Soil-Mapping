# Isabel Luotto 13/06/22

# This is a comment
rm(list=ls())
gc()
# setting the work environment ----
setwd('C:/Users/hp/Documents/GitHub/Digital-Soil-Mapping/')

# load packages

# Create first object
myfirstobject = 42
myfirstobject <- 42
myfirstobject

# create a vector that goes from 5 to 100 ':'
(a <- 5:100) #In this case we're creating a short #vector counting from 5 to 10
#Let's multiply x with all the numbers between 5 #and 10
x =142
b2 <- x*a
b #type b to see in the console what this #simple vector contains 
#sqrt of myfirstobject?
sqrt(x)

#Set seed is a function to make a 
#random sample reproducible
?set.seed #we can get some info in 
#the help pane by using "?"

set.seed(13622)

(vec<- rnorm(1000,2,1))

?rnorm #if you leave something between the commas blank
#it goes to default
b<-rnorm(120, 2,)
hist(b)
b
boxplot(b)


data <- c(2,3) 
Data #For R data and Data are two different things
#Error: object 'Data' not found
library(raster)#It's library not Library
#Error in Library(raster) : could not find function "Library"

numbers <- c(-5,6)
#Error in numbers[-5, 6] : incorrect number

(dat <- 1:42)
dat
dat <- dat[- 4,5]#[]are indices  
dat <- dat[-c(4,5)];dat#this is a comment
dat <- 2:42 
dat
dat[c(4,5)] <- c(2,3)
vec <- c('cat','mouse','dog','giraffe')
vec[c(2,3)]
vec[2] <-'bird'

dat <- 2:42 
dat <- dat[-c(8,41)];dat 

# Example with the iris dataset ----
#load in the iris example dataset
dat <- data.frame(iris)

View(dat)#Opens up a window showing you the whole dataset. 
#You can also open it from the global environment by clicking on it
head(dat,10) #first 6 observations. tail(dat) gives you the last 6
names(dat) #Tells you the names and position of each variable

str(dat)#how many variables, observations, the class/type of data
summary(dat)#Tells you about the distribution of the data 
#Minimum, Mean,Median,Maximum of each variable

#Visually explore specific columns
hist(dat$Sepal.Width)#The $ sign allows you
#to select specific columns
boxplot(dat$Petal.Length)

plot(dat$Petal.Length, dat$Sepal.Length)

#paste0 no space e.g. viriginica2
# paste with space e.g. viriginica 2
dat$Species2 <- paste(dat$Species, 2)#with the $ sign you can #create a new column. The paste function allows you to add #sequences to characters. This is useful when creating an ID

str(dat)
dat$Species2 <- as.factor(dat$Species2)
str(dat)

names(dat)#check the name and position of each variable in the #console

#names(dat)[6] <- 'new_name'
names(dat) <- c('sl','sw','pl','pw','s','s2')
#names(dat)[c(2,3)] <- c('sw2','pl2')

mean(dat$pl)

# Write data into R ----
write.csv(dat, file= "iris.csv")#Put the 
#name of the file in ""
#If you want to save the file somewhere else other than the wd you #just have to specify where
write.csv(dat, file= "folder_path/iris2.csv")

# Export dataframe as excel sheet using the writexl package
#install.packages('readxl') # to read excel sheets
#install.packages("writexl")# to write excel sheets
library(readxl)
library(writexl)

write_xlsx(dat, 'iris.xlsx')

# Read data into R ----
iris2 <- read.csv('01-Data/horizons.csv')
#NOT RECOMMENDED
#iris2 <- read.csv(file.choose())

# Conditional selection in dataframes ----
dat[1:6, 5]
head(dat)
#by leaving the section after the comma blank we're selecting all the columns
dat <- dat[ , -6]
head(dat)#let's get rid of the last column

#Create a random subset with 85% of the data

smp_size <- floor(0.85 * nrow(dat))#define the subsample size
train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
train <- dat[train_ind, ] 
test <- dat[-train_ind, ]

write.csv(train, file = "train_iris.csv")
write.csv(test, file = "test_iris.csv")


# Conditional selection based on logical operators ----
subset <- dat[dat$pl> 5,]

# AND operator all conditions need to be TRUE
subset = dat[dat$pl >5 & dat$s == "virginica"& dat$pw >= 2, c('pl', 's')]

# Or operator not all conditions need to be TRUE
subset2 <- dat[dat$pl >5 | dat$sl > 3, ]

# example with the operator %in%
levels(dat$s)
dat[dat$s %in% c('setosa','virginica'), ]

# ifelse
dat$s2 <- ifelse(dat$s == 'virginica','v',as.character(dat$s))
# double ifelse
dat$s3 <- ifelse(dat$s == 'virginica',
                 'v',ifelse(dat$s =='setosa','s',
                  as.character(dat$s)))

# ifelse with %in%
dat$s2 <- ifelse(dat$s %in% c('setosa','virginica'),
                 'pretty',as.character(dat$s))

# Remove outliers using a boxplot ----
dat[c(2,4,7,9), 'pw'] <- c(-10,20,15,-10)
bx<- boxplot(dat$pw)
bx$out

out <-boxplot(dat$pw, range = 1.5, plot= FALSE)$out
out
#Exclude extreme observations
dat <- dat[-which(dat$pw %in% out),] 

boxplot(dat$pw)
View(dat)

# Remove missing values (NAs) ----
dat <- data.frame(iris)
dat <- setNames(dat, c("SL","SW","PL","PW","S"))
dat$PL2 <- dat$PL#first let's create a column with NAs
dat[c(9,3,12,6),6]=NA 

summary(dat)#the summary function is good 
#to check if there are any NAs

dat <- dat[!is.na(dat$PL2),]
# Preferred option (also the only one when using)
# data.tables instead of dataframes (data.table package)
dat <- dat[complete.cases(dat),]#the complete
#cases function selects all the rows without 
#NAs

dat <- data.frame(iris) 
#dat$Species <- as.factor(dat$Species)
par(bty="u", family = "mono")#we want an L shaped graph with #font "mono" 
boxplot(dat$Petal.Length ~ dat$Species, 
        main = "Petal length by Species",
        col.main = "#009999",
        ylab = "Petal length [cm]",xlab= '')
