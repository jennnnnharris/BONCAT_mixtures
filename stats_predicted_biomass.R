# stats for predicted verse measured biomass 

# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
#library(lme4)
#library(nlme)
#library(emmeans)
#library(multcomp)
#library(dplyr)


###############write functions for predictions###############
get.predict<-function(df, sp1, sp2, trait, sp3) {
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
    print(trait)
    sp1.df<- filter(df, Treatment==sp1) %>% select(all_of(trait))
    sp1.predict<-sp1.df/3
    #print(sp1.df)
    #print(sp1.predict)
    sp2.df<-filter(df, Treatment==sp2) %>% select(all_of(trait))
    sp2.predict<-sp2.df/3
    print(sp2.df)
    print(sp2.predict)
    sp3.df<-filter(df, Treatment==sp3) %>% select(all_of(trait))
    sp3.predict<-sp3.df/3
    #print(sp3.df)
    #print(sp3.predict)
    predict <- sp1.predict + sp2.predict + sp3.predict
    print(predict)
    return(predict)
  }
  
}


get.shoot.predict<-function(df, sp1, sp2, sp3) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(Shoot.Biomass)
    sp1.predict<-sp1.df$Shoot.Biomass/2
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$Shoot.Biomass/2
    predict <- sp1.predict + sp2.predict
    print(predict)
    return(predict)
  }
  else
    
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(Shoot.Biomass)
  sp1.predict<-sp1.df$Shoot.Biomass/3
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$Shoot.Biomass/3
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$Shoot.Biomass/3
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  
}
get.root.predict<-function(df, sp1, sp2, sp3) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(Root.Biomass)
    sp1.predict<-sp1.df$Root.Biomass/2
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$Root.Biomass/2
    predict <- sp1.predict + sp2.predict
    print(predict)
    return(predict)
  }
  else
    
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(Root.Biomass)
  sp1.predict<-sp1.df$Root.Biomass/3
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$Root.Biomass/3
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$Root.Biomass/3
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  
}


#####predictions from monocultures for biomass #######

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
#library(lme4)
#library(nlme)

#### import biomass data and process
biomasspath <- "C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology"
setwd(biomasspath)
df <- read.csv("biomass_potlevel.csv", row.names = 1) # biomass data

# process to make n species column, summarize at pot level, make treatment a factor.
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<- df %>% select(c(Brassicae, Legume, Grass )) %>% rowSums()
#df<-df %>% group_by(Trt_ID, Treatment, Rep, Brassicae, Legume, Grass, N, n_species ) %>% summarise(Root.Biomass = sum(Total.Root.g), Shoot.Biomass = sum(Stem.Biomass.g), )
# rename cols
# rename cols
df$Shoot.Biomass<-df$Stem.Biomass.g
df$Stem.Biomass.g=NULL
df$Root.Biomass = df$Total.Root.g
df$Root.Biomass.g=NULL
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# rename
df <- df
head(df)

## we expect that a plant makes the same amount of biomass in monoculures vs mixtures
#example
#LB biomass = L monoculture/  2 + B monocultre/2 

#get shoot predictions
df<-df %>% group_by(Rep, N)
df

GB<-get.shoot.predict(df, "G", "B") 
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

# add N
dim(predict)
predict$Nitrogen_label<-rep(rep(c("Nitrogen +", "Nitrogen -"), each=6), 4)

## add to df

df1<-full_join(df, predict) 
df1<-df1 %>% ungroup()
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB.predict", "GB", "LB.predict",  "LB",   "LG.predict", "LG", 
                                                "LGB.predict"   , "LGB" ))

#### plot predictions from monocultures for biomass

mycols <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
  "grey",  
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338" # medium mocha brown LGB
)

df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
  )+
    labs(title = "A",
       x="",
       y="Root biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #

df1  %>% 
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",
  )+
  labs(title = "B",
       x="",
       y="shoot biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #



######### anova#####
# gott include nitrogen + or minus in the predicted values
# Analysis of variance percent N from BNF
library(multcomp)
library(emmeans)
# include nitrogen in this model

# root biomass
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Root.Biomass~ Treatment*Nitrogen_label, data=df2)
summary(m1)
anova(m1)
anov = aov(m1)
summary(anov)


emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)



#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- aov(Root.Biomass~ Treatment*N, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

# shoots
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)


