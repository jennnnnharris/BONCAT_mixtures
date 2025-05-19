# Mixtures Boncat
# created: Mau 2025
# last edited: May 2025
# author: Jennifer Harris
# data: plant biomass


rm(list=ls())
library(readxl)
library(tidyverse)
library(lubridate)

#import colors
mycols7<-c( "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF",  "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
mycols8<-c("grey", "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF",  "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")



# import data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")

df<-read_excel("Rice_greenhouse_ccexp_biomass.xlsx", sheet = 1)

head(df)

# make n species column
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<-rowSums(select(df, Brassicae, Legume, Grass))

head(df)

# summarize at pot level

df<-df %>% group_by(Pot.Number, Treatment, Rep, Brassicae, Legume, Grass, n_species ) %>% summarise(Root.Biomass = sum(Total.Root.g), Shoot.Biomass = sum(Stem.Biomass.g), )
df$Root.to.Shoot <- df$Root.Biomass / df$Shoot.Biomass


# make treatment and day factors
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))


# plot for n species

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
#svg(file="biomass.svg",width = 7, height=4)

p1<-df  %>%
  ggplot(aes(x=n_species, y=Root.Biomass )) +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= blues[2])+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)

p2<-df  %>%
  ggplot(aes(x=n_species, y=Shoot.Biomass )) +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= blues[3])+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))

p3<-df  %>%
    ggplot(aes(x=n_species, y=Root.to.Shoot )) +
    geom_jitter(width = .2, size=1 )+
    geom_smooth(method = lm, color= blues[4])+
    theme_classic(base_size = 14)+
    theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
          plot.title = element_text(hjust = 0.5))

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="biomass.species.svg",width = 10, height=4)
require(gridExtra)
#windows(8,4)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()


# plot for each treatment

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
#svg(file="biomass.svg",width = 7, height=4)

p1<-df  %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  #ylab("percent active")+
  facet_grid( ~n_species, scales = "free", space = "free")
  #ggtitle("number of species")
  #geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
  #geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)

dev.off()  
#


p2<-df  %>%
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  #ylab("")+
  facet_grid( ~n_species, scales = "free", space = "free")
  #ggtitle("number of species")
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)

dev.off()  

p3<-df  %>%
  ggplot(aes(x=Treatment, y=Root.to.Shoot, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  #ylab("")+
  facet_grid( ~n_species, scales = "free", space = "free")
  #ggtitle("number of species")
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)

dev.off()  


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="biomass.svg",width = 10, height=4)
require(gridExtra)
#windows(8,4)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



#### normalize to monoculures
#### kind of like reaction norms

#example
#LB = L/ 2 + B/2 

df1<-df %>% group_by(Treatment) %>% summarise(Root.Biomass = mean(Root.Biomass), Shoot.Biomass = mean(Shoot.Biomass), Root.to.Shoot = mean(Root.to.Shoot), )
head(df1)

GB<-filter(df1, Treatment=="G")[,2:4]/3 + filter(df1, Treatment=="B")[,2:4]/2 
LB<-filter(df1, Treatment=="L")[,2:4]/2 + filter(df1, Treatment=="B")[,2:4]/2 
LG<-filter(df1, Treatment=="L")[,2:4]/2 + filter(df1, Treatment=="G")[,2:4]/2 
LGB<-filter(df1, Treatment=="L")[,2:4]/3 + filter(df1, Treatment=="G")[,2:4]/3 + filter(df1, Treatment=="B")[,2:4]/3 
val<- rbind(GB, LB, LG, LGB)
Treatment <- c("GB.predict", "LB.predict", "LG.predict", "LGB.predict")
df1<-cbind(Treatment, val)
head(df1)
 
## ad to df
df1<-full_join(df1, (df%>%filter(n_species!="1")))

mycols<-c("#365C83FF", "grey" ,"#384351FF", "grey", "#4D8F8BFF", "grey", "#CDD6ADFF", "grey")
as.factor(df1$Treatment)

##plot
p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  ylab("Shoot Biomass")
 


p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
 
p3<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.to.Shoot, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="overyielding.biomass.svg",width = 10, height=4)
require(gridExtra)
#windows(8,4)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()
