# summarise biomass to pot level
# Dec 10 2025

#load libraries
library(readxl)
library(tidyverse)

#### import biomass data and process #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read_excel("biomass_species.xlsx") # biomass data

df<-df %>% group_by(Trt_ID, Treatment, N, Rep, Brassicae, Legume, Grass ) %>%
  summarise(
    Stem.Biomass.g = sum(Stem.Biomass.g),
    Root.Biomass.g = sum(Root.Biomass.g),
    Bulk.Root.g = sum(Bulk.Root.g),
    Total.Root.g = sum(Total.Root.g)
  )

df  


write.csv(df, biomass_potlevel.csv)
