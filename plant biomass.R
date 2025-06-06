# Mixtures Boncat
# created: Mau 2025
# last edited: May 2025
# author: Jennifer Harris
# data: plant biomass


rm(list=ls())
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)

#import colors
mycols7<-c( "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF",  "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
mycols8<-c("grey", "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF",  "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")


# write functions
get.predict<-function(df, sp1, sp2, sp3, trait) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(all_of(trait))
    #print(sp1.df)
    sp1.predict<-sp1.df/2
    print(sp1.predict)
    sp2.df<-filter(df, Treatment==sp2) %>% select(all_of(trait))
    sp2.predict<-sp2.df /2
    print(sp2.predict)
    predict <- sp1.predict + sp2.predict
    print(predict)
    #return(predict)
  }
  else {
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(all_of(trait))
  sp1.predict<-sp1.df/2
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$trait/2
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$trait/2
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  }
  
}

get.predict(df, "G", "B", trait="Shoot.Biomass")
t<-"Shoot.Biomass"
df %>% select(all_of(t))

get.root.predict<-function(df, sp1, sp2, sp3) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(Root.Biomass)
    sp1.predict<-sp1.df$Root.Biomass/3
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$Root.Biomass/3
    predict <- sp1.predict + sp2.predict
    print(predict)
    return(predict)
  }
  else
    
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(Root.Biomass)
  sp1.predict<-sp1.df$Root.Biomass/2
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$Root.Biomass/2
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$Root.Biomass/2
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  
}


# import data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read.csv("Rice_greenhouse_ccexp_biomass_block.csv")
head(df)

# data wrangling
# make n species column
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<-rowSums(select(df, Brassicae, Legume, Grass))
head(df)

# summarize at pot level
df<-df %>% group_by(Pot.Number, Treatment, Rep, Brassicae, N, Legume, Grass, n_species, Block ) %>% summarise(Root.Biomass = sum(Total.Root.g), Shoot.Biomass = sum(Stem.Biomass.g), )
df$Root.to.Shoot <- df$Root.Biomass / df$Shoot.Biomass


# make treatment and day factors
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))


# plot for n species
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
#svg(file="biomass.svg",width = 7, height=4)

p1<-df  %>%
  ggplot(aes(x=n_species, y=Shoot.Biomass )) +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= blues[2])+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)
p1

p2<-  df  %>%
  ggplot(aes(x=n_species, y=Root.Biomass )) +
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
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0, size=14))+
  facet_grid( ~n_species, scales = "free", space = "free")
  #ggtitle("A Root biomass")
  #geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
  #geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)
p1

p2<-df  %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0, size= 14))+
  facet_grid( ~n_species, scales = "free", space = "free")
p2



p3<-df  %>%
  ggplot(aes(x=Treatment, y=Root.to.Shoot, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0, size=14))+
  ggtitle("C root:shoot")+
  facet_grid( ~n_species, scales = "free", space = "free")
  

p3


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="biomass.svg",width = 10, height=4)
require(gridExtra)
#windows(9,4)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

####linear model#####

#number of species
#lm
m1<-lm(Shoot.Biomass~n_species*N*Block ,data=df)
summary(m1)
plot(m1)
# number of species
#Adjusted R-squared:  0.002231 
#F-statistic: 1.092 on 2 and 80 DF,  p-value: 0.3406

#Legume
m1<-lm(Shoot.Biomass~ Legume*Block*N ,data=df)
summary(m1)
plot(m1)

#Leg +Brass
m1<-lm(Shoot.Biomass~Brassicae*Legume*Block*N ,data=df)
summary(m1)
plot(m1)

#Leg +Grass
m1<-lm(Shoot.Biomass~Grass*Brassicae*Legume*Block*N ,data=df)
summary(m1)
plot(m1)

#Grass + Legume
m1<-lm(Shoot.Biomass~Grass*Legume*Block*N ,data=df)
summary(m1)
plot(m1)

###roots
#plots
df  %>%
  ggplot(aes(x=n_species, y=Root.Biomass )) +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= blues[3])+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~N)

#n species
m1<-lm(Shoot.Biomass~n_species*N*Block ,data=df)
summary(m1)

#Legume
m1<-lm(Shoot.Biomass~ Legume*Block*N ,data=df)
summary(m1)
plot(m1)


###root2shoot
#plots
df  %>%
  ggplot(aes(x=n_species, y=Root.to.Shoot )) +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= blues[3])+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Block)

#n species
m1<-lm(Shoot.Biomass~n_species*N*Block ,data=df)
summary(m1)

#Legume
m1<-lm(Shoot.Biomass~ Legume*Block*N ,data=df)
summary(m1)
plot(m1)


####improvement over prediction #######

## we expect that a plant makes the same amount of biomass in monoculures vs mixtures
#example
#LB biomass = L monoculture/ 2 + B monocultre/2 

#get shoot predictions
df<-df %>% group_by(Rep, N)

GB<-get.predict(df, "G", "L" , "G", Shoot.Biomass) 


LB<-get.shoot.predict(df, "L", "B") 
LG<-get.shoot.predict(df, "L", "G") 
LGB<-get.shoot.predict(df, "L", "G", "B")  
Shoot.Biomass<- round(as.numeric(c(GB, LB, LG, LGB)), 2)
Shoot.Biomass
Treatment<-c(rep("GB.predict", n_groups(df)), rep("LB.predict", n_groups(df)), rep("LG.predict", n_groups(df)), rep("LGB.predict", n_groups(df)) )
Treatment
predict <- as.data.frame(cbind(Treatment, Shoot.Biomass))
predict$Shoot.Biomass<-as.numeric(predict$Shoot.Biomass)
predict

#get root predictions

GB<-get.root.predict(df, "G", "B") 
LB<-get.root.predict(df, "L", "B") 
LG<-get.root.predict(df, "L", "G") 
LGB<-get.root.predict(df, "L", "G", "B")  
Root.Biomass<- round(as.numeric(c(GB, LB, LG, LGB)), 2)

predict<- cbind(predict, Root.Biomass)
predict$Root.Biomass<-as.numeric(predict$Root.Biomass)
predict
 

## add to df
df1<-df %>% filter(n_species!="1")
mixtures<-df %>% filter(n_species!="1")
47*2

df1<-full_join(df1, predict) 
df1<-df1 %>% ungroup()

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
 
p1

p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
p2
 

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="overyielding.biomass.svg",width = 10, height=4)
require(gridExtra)
windows(8,4)
grid.arrange(p1, p2, ncol=2)
dev.off()

####tukey test###

df1
anova()
TukeyHSD()