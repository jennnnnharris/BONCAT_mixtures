# nod analysis
# j harris 6/19/2025


#load packages#
library(tidyverse)
library(lubridate)
library(readxl)

#set working directory
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology/nod_pictures/cropped/splits/threshold/Result")


#import data frames that start with summary 

#make a list of files
files<-list.files( pattern = "summary_blue*")
files

# use function to import files. 
dfs<-lapply(files,read.csv)

# join the dataframes
df<-do.call(rbind, dfs)
df

# reformat treatment code
# remove the blue
trt<-df$Slice
trt<-sub("blue_", "", trt)
trt<-sub(".png", "", trt)
trt<-sub("-1", "", trt)
trt<-sub("N", "", trt)
df$Slice<-trt
df


# split into a few columns
x<-strsplit(trt, "_")
x<-do.call(rbind, x)
df<-cbind(x,df)
df

# rename columns
# first part is treatment
# secound is Nitrogen
# last is rep
n<-colnames(df)
n<-c("Treatment", "Nitrogen", "Rep", "Trt_ID", n[5:length(n)])
colnames(df) <- n
df


# reformat the treatment code to match the long name for the DF
# currenntly:    L_1_2
# goal: L+N2
trt<-df$Trt_ID
#change first + to +
trt<-sub("_0_", "", trt)
trt
# change teh string  "1" to +N
trt<-sub("_1_", "+N", trt)
trt
trt
df$Trt_ID <- trt
head(df)


# add weight data

# import nod rhizo weights
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
wt<-read_xlsx("nod_rhizo_weights.xlsx", sheet = 2)
head(wt)

# bind with nodule df
df1<-left_join(wt, df)
df1

# save out put

# save as a new data file for analysis
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")

write.csv(df1, "noduledata.csv")



