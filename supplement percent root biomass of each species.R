

#### supplement: percent root biomass of each species ####

# load libraries and cols
library(readxl)
library(tidyverse)
IBM <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
  "#DC267F", # magenta pink GB
  "#FE6100", # bright orange LB
  "#FFB000", # golden yellow LG
  "#865338" # medium mocha brown LGB
)
mono_cols <- 
  c( #IBM colors
    "navy", # dark royal blue L
    "#648FFF", # french blue G
    "#785EF0" # light purple B
  )
# import data frame 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read_excel("biomass_species.xlsx")

# tidy data
head(df)

# make plot
ggplot(df, aes(fill=Species, y=Root.Biomass.g, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# make total.root.species.g column
total<-df %>% group_by(Trt_ID) %>%
  summarise(
    Total.no.bulk = sum(Root.Biomass.g))

# add col to df
df<-left_join(df, total)

# make percent col
df<-df%>% mutate(
  percent= (Root.Biomass.g/Total.no.bulk)*100)

# make plot percent
ggplot(df, aes(fill=Species, y=percent, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#add unknown bulk to df

bulk<-df %>%
  dplyr::select(Treatment, N, Rep, Trt_ID,  Bulk.Root.g) %>%
  group_by(Trt_ID, Treatment, N, Rep) %>%
  summarise(
    Root.Biomass.g= sum(Bulk.Root.g)) %>%
  mutate(
    Species="bulk",
    Brassicae= 0,
    Legume=  0,
    Grass= 0
  )


df1 <- df %>% dplyr:: select( Treatment, N, Rep, Trt_ID, Species,  Brassicae, Legume, Grass, Root.Biomass.g )
df1<-full_join(df1,bulk)

#add total+bulk to df 

total<-df %>%
  dplyr::select(Trt_ID, Total.Root.g) %>%
  group_by(Trt_ID) %>%
  summarise(
    Total.Root.g= sum(Total.Root.g))

df1<-left_join(df1, total)

# calculate percent
df1<-df1 %>% mutate(
  percent = Root.Biomass.g/Total.Root.g
)

df1$Species<- factor(df1$Species, levels= c("bulk", "legume", "grass", "brassica"))
# make plot percent
ggplot(df1, aes(fill=Species, y=percent, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


### bulk assigned to crops ####

# import data frame 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read_excel("biomass_species.xlsx")

# tidy data
head(df)

# make plot
ggplot(df, aes(fill=Species, y=Root.Biomass.g, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# make total.root.species.g column
total<-df %>% group_by(Trt_ID) %>%
  summarise(
    Total.withbulk = sum(Total.Root.g))

# add col to df
df<-left_join(df, total)

# make percent col
df<-df%>% mutate(
  percent= (Total.Root.g/Total.withbulk)*100)

df$Species<- factor(df$Species, levels= c("legume", "grass", "brassica"))

# make plot percent
mono_cols <- 
  c( #IBM colors
    "#0a3170", # dark royal blue L
    "#cdddf7", # french blue G
    "#785EF0" # light purple B
  )

# recode name
df1<-df%>% filter(Treatment!="B" & Treatment!="L" & Treatment!="G") 

df1$label<-gsub("LGB", "" ,df1$Trt_ID)
df1$label
df1$label<-gsub("GB", "" ,df1$label)
df1$label<-gsub("LG", "" ,df1$label)
df1$label<-gsub("LB", "" ,df1$label)




setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement")
svg("percent.biomass.svg", width = 8, height =4.5 )
df1%>%
  ggplot(aes(fill=Species, y=percent, x=label)) + 
  theme_bw(base_size = 12)+
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values=mono_cols)+
  facet_grid(~Treatment, space="free", scales="free")+
  labs(x= "Sample",
       y="proportion dry biomass (g)")

dev.off()
write.csv(df, "percent.biomass.csv")

