# Nodule
# J Harris


# clear workspace
rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

#import packages
library(tidyverse)
library(lme4)
require(gridExtra)

# import colors
#import colors
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")
mycols7<-c( "#4B2D4BFF", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
mycols4 <- c("#4B2D4BFF",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "LB", "LG", "LGB"))



# import data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
metadat <- read.csv("metadata.csv")
df <- read.csv("plant_physiology/noduledata.csv", row.names = 1)
df

# combine data and metadata
df<-left_join(df, metadat)
head(df)


##### 1. Does number of species effect nodule traits? #####

# basic data viz
hist(df$Count)
plot(df$n_species, df$Count)
plot(df$n_species, df$Total.Area)
plot(df$n_species, df$Nodule_weight)

# prettier plot

p2<-df  %>%
  ggplot(aes(x=n_species, y=Total.Area )) +
  geom_jitter(width = .2, size=2 )+
  geom_smooth(method = lm, color= "grey5")+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  ylab("Nodule Area")
p2


#save image
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Figure1_nspecies")
svg(file="nodule.species.svg",width = 3, height=3)
p2
#grid.arrange(p1, p2, p3, ncol=1)
dev.off()

# with colors
p2<-df  %>%
  ggplot(aes(x=n_species, y=Total.Area, colour = Treatment )) +
  geom_jitter(width = .2, size=2 )+
  geom_smooth(method = lm, color= "grey10")+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = mycols4)+
  ylab("Nodule Area")
p2


#save image
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Figure1_nspecies")
svg(file="col.nodule.species.svg",width = 3, height=3)
p2
#grid.arrange(p1, p2, p3, ncol=1)
dev.off()


# linear model
m1<-lm(Count~n_species ,data=df)
summary(m1)
plot(m1)
# not significant

m1<-lm(Total.Area~n_species ,data=df)
summary(m1)
plot(m1)
# not significant

m1<-lm(Nodule_weight~n_species*N*Block ,data=df)
summary(m1)
plot(m1)
# not significant



########### 2. Does cover crop treatment effect nodule traits? ###########


# basic data viz
hist(df$Count)
hist(df$Total.Area)
boxplot(df$Count~df$Treatment)
boxplot(df$Total.Area~df$Treatment)

# prettier plot
fig2path <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Figure2_treatment"
setwd(fig2path)

svg(file="nodule.trt.svg",width = 4, height=4)
df  %>% filter(Treatment!="NA") %>%
  ggplot(aes(x=Treatment, y=Total.Area, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols4) +
  scale_fill_manual(values = mycols4)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0, size=14))
dev.off()

# count
df  %>% filter(Treatment!="NA") %>%
  ggplot(aes(x=Treatment, y=Count, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0, size=14))


# weight
df  %>% filter(Treatment!="NA") %>%
  ggplot(aes(x=Treatment, y=Nodule_weight, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0, size=14))



# anova count
df$Treatment   <- factor(df$Treatment, levels= c( "L", "LB", "LG", "LGB"))
m1 <- aov(Count ~  Treatment+Block+N, data =df)
summary(m1)
plot(m1)

# SA
m1 <- aov(Total.Area~  Treatment+Block+N, data =df)
summary(m1)
plot(m1)
# significant

# WT
m1 <- aov(Nodule_weight ~ Treatment +Block+N, data =df)
summary(m1)
plot(m1)


# 3. Do mixtures preform better than predicted?


