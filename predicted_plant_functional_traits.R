# Figure 3 prediction


# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)


IBM <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
  "#DC267F", # magenta pink GB
  "#FE6100", # bright orange LB
  "#FFB000", # golden yellow LG
  "#865338" # medium mocha brown LGB
)
legume_cols <- c( #IBM colors
  "navy", # dark royal blue L
  "#FE6100", # bright orange LB
  "#FFB000", # golden yellow LG
  "#865338" # medium mocha brown LGB
)

# load paths
biomasspath <- "C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology"


#### import biomass data and process #####
setwd(biomasspath)
df <- read.csv("biomass_potlevel.csv") # biomass data

# process to make n species column, summarize at pot level, make treatment a factor.
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<- df %>% select(c(Brassicae, Legume, Grass )) %>% rowSums()
df<-df %>% group_by(Trt_ID, Treatment, Rep, Brassicae, Legume, Grass, N, n_species ) %>% summarise(Root.Biomass = sum(Total.Root.g), Shoot.Biomass = sum(Stem.Biomass.g), )
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# give long composition name
# make long_name var
# df$long_name <- df$Treatment
# df$long_name<-gsub("L", "Legume", df$long_name)
# df$long_name<-gsub("G", "Grass", df$long_name)
# df$long_name<-gsub("B", "Brassica", df$long_name)
# df$long_name<-gsub("LegumeBrassica", "Legume_Brassica", df$long_name)
# df$long_name<-gsub("LegumeGrass", "Legume_Grass", df$long_name)
# df$long_name<-gsub("GrassBrassica", "Grass_Brassica", df$long_name)
# df$long_name 
# df$long_name<- factor(df$long_name, levels = c("Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))

# rename
df <- df
head(df)




###############write functions###############
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


## add to df
df1<-df %>% filter(n_species!="1")
mixtures<-df %>% filter(n_species!="1")
47*2

df1<-full_join(df1, predict) 
df1<-df1 %>% ungroup()
as.factor(df1$Treatment)

#### plot predictions from monocultures for biomass
#old
# new
#mycols<-c("#FEB5A2FF","grey" , "#9D7660FF", "grey", "#D7B5A6FF", "grey", "#3896C4FF" , "grey")
mycols <- c( #IBM colors
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)

#setwd(fig3path)
#svg(file="biomass.root.predict.svg",width = 2.8, height=3)
label <- df1$Treatment
label <- gsub("LGB.predict", "*" ,label )
label <- gsub("LB.predict", "*" ,label )
label <- gsub("LG.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label

p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  
  geom_text(y=6, label =label , nudge_x = -.8, size=8)


p2

label <- df1$Treatment
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "*" ,label )


#setwd(fig3path)
#svg(file="biomass.shoot.predict.svg",width = 4, height=3)
p3<-df1  %>% 
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  #geom_text(y=8, label =label , nudge_x = -.8, size=7)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))

p3

# put the tow plots together

require(gridExtra)
grid.arrange(p2, p3, ncol=2)



######### example calculation  #####


mycols7<-c("#466F9DFF", "#91B3D7FF",  "white", "white", "white", "#D7B5A6FF", "#3896C4FF" )

p1<-df  %>% filter(n_species!="NA") %>%
  filter(Treatment!= "LGB") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
    xlab("") 
p1


# caculated half for L
sp1.df<- filter(df, Treatment=="L") %>% select(Root.Biomass)
sp1.predict<-sp1.df$Root.Biomass/2
Root.Biomass<-sp1.predict
Treatment <- rep("half_L", length(Root.Biomass))
halfl<-data.frame(Treatment, Root.Biomass)
df1<-full_join(df, halfl)
df1<-df1 %>% filter(Treatment!="LGB" & Treatment!="LB" & Treatment!="B")
df1$Treatment<-factor(df1$Treatment, levels=c("L", "G", "half_L", "GB", "LG"))


# caculated half for G
sp1.df<- filter(df, Treatment=="G") %>% select(Root.Biomass)
sp1.predict<-sp1.df$Root.Biomass/2
Root.Biomass<-sp1.predict
Treatment <- rep("half_G", length(Root.Biomass))
halfg<-data.frame(Treatment, Root.Biomass)
df1<-full_join(df1, halfg)
df1<-df1 %>% filter(Treatment!="GB")
df1$Treatment<-factor(df1$Treatment, levels=c("L", "G", "half_L", "half_G", "LG"))
as.factor(df1$Treatment)



# caculated LG 
sp1<- filter(df1, Treatment=="half_L") %>% select(Root.Biomass)
sp2<- filter(df1, Treatment=="half_G") %>% select(Root.Biomass)

LG.predict <- sp1$Root.Biomass + sp2$Root.Biomass
Root.Biomass<-LG.predict
Treatment <- rep("LG.predict", length(Root.Biomass))
predict<-data.frame(Treatment, Root.Biomass)
df1<-full_join(df1, predict)
df1$Treatment<-factor(df1$Treatment, levels=c("L", "G", "half_L", "half_G", "LG.predict", "LG"))
as.factor(df1$Treatment)

### plot###
mycols7<-c( "#715b8a", "#4D8F8BFF", "grey", "grey", "grey", "#AD5A6BFF") 

mycols <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "grey", # light purple B
  "grey", # magenta pink GB
  "grey", # bright orange LB
  "#FFB000" # golden yellow LG
)

p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
  xlab("") 
p1




######### anova#####


#GB
#filter
df2<-dfb1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LB
df2<-dfb1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LG
df2<-dfb1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LGB
df2<-dfb1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

# shoots
#GB
#filter
df2<-dfb1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LB
df2<-dfb1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LG
df2<-dfb1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LGB
df2<-dfb1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)




######### make df of difference between predicted and not predicted #######
#get shoot predictions
df<-df %>% group_by(Rep, N)
df

GB<-get.shoot.predict(df, "G", "B") 
LB<-get.shoot.predict(df, "L", "B") 
LG<-get.shoot.predict(df, "L", "G") 
LGB<-get.shoot.predict(df, "L", "G", "B")  
predict.Shoot.Biomass<- round(as.numeric(c(GB, LB, LG, LGB)), 2)
Treatment<-c(rep("GB", n_groups(df)), rep("LB", n_groups(df)), rep("LG", n_groups(df)), rep("LGB", n_groups(df)) )

predict <- as.data.frame(cbind(Treatment, predict.Shoot.Biomass))
predict$predict.Shoot.Biomass<-as.numeric(predict$predict.Shoot.Biomass)
predict

#get root predictions

GB<-get.root.predict(df, "G", "B") 
LB<-get.root.predict(df, "L", "B") 
LG<-get.root.predict(df, "L", "G") 
LGB<-get.root.predict(df, "L", "G", "B")  
predict.Root.Biomass<- round(as.numeric(c(GB, LB, LG, LGB)), 2)

predict<- cbind(predict, predict.Root.Biomass)
predict$predict.Root.Biomass<-as.numeric(predict$predict.Root.Biomass)
predict

N<-rep(c(rep(1, 6), rep(0,6)), 4)
Rep<-rep(1:6, 8)

predict$Rep <- Rep
predict$N <- N

## add to df
df1<-df %>% filter(n_species!="1")

dim(df1)
dim(predict)
predict
df2<-full_join(df1, predict) 
df2<-df2 %>% ungroup()


# make df of difference between predicted and not predicted
# actual - predicted 


df2<-df2 %>%
mutate( root.difference=  Root.Biomass-predict.Root.Biomass,
        shoot.difference = Shoot.Biomass-predict.Shoot.Biomass,
        Nitrogen=N)


write.csv(df2, "predicted.biomass.csv", row.names = FALSE)



######predictions weed seed########
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
# group by
df<- df %>% group_by(Rep, Nitrogen)
df

# pigdf
LG<-get.predict(df, "L", "G", "pigweed_prop_nongerm")
LB<-get.predict(df, "L", "B", "pigweed_prop_nongerm")
GB<-get.predict(df, "G", "B", "pigweed_prop_nongerm")
LGB<-get.predict(df, "L", "G", "pigweed_prop_nongerm", "B")
LG

predict<-rbind(GB, LB, LG, LGB)
predict$pigweed_prop_nongerm<-round(as.numeric(predict$pigweed_prop_nongerm), 3)
predict$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict


# foxtail
LG<-get.predict(df, "L", "G", "foxtail_prop_nongerm")
LB<-get.predict(df, "L", "B", "foxtail_prop_nongerm")
GB<-get.predict(df, "G", "B", "foxtail_prop_nongerm")
LGB<-get.predict(df, "L", "G", "foxtail_prop_nongerm", "B")
LG

predict1<-rbind(GB, LB, LG, LGB)
predict1$foxtail_prop_nongerm<-round(as.numeric(predict1$foxtail_prop_nongerm), 3)
predict1$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict1
predict<-full_join(predict, predict1)

## add to df
df1<-df %>% filter(n_species!="1")

df1<-full_join(df1, predict) 
df1<-df1 %>% ungroup()
as.factor(df1$Treatment)

#### plot predictions from monocultures for biomass
#mycols<-c("#FEB5A2FF","grey" , "#9D7660FF", "grey", "#D7B5A6FF", "grey", "#3896C4FF" , "grey")
mycols <- c( #IBM colors
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)


p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
p2

p3<-df1  %>% 
  ggplot(aes(x=Treatment, y=pigweed_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  #geom_text(y=8, label =label , nudge_x = -.8, size=7)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))

p3

# put the two plots together
require(gridExtra)
grid.arrange(p2, p3, ncol=2)

#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(foxtail_prop_nongerm~ Treatment, data=df2)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(foxtail_prop_nongerm~ Treatment, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(foxtail_prop_nongerm~ Treatment, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(foxtail_prop_nongerm~ Treatment, data=df2)
anova(m1)

# pigweed
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(pigweed_prop_nongerm~ Treatment, data=df2)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(pigweed_prop_nongerm~ Treatment, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(pigweed_prop_nongerm~ Treatment, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(pigweed_prop_nongerm~ Treatment, data=df2)
anova(m1)

##########weed prop non germinated #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
# group by
df<- df %>% group_by(Rep, Nitrogen)
df

# pigweed non germ 
LG<-get.predict(df, "L", "G", "pigweed_num_nongerm")
LB<-get.predict(df, "L", "B", "pigweed_num_nongerm")
GB<-get.predict(df, "G", "B", "pigweed_num_nongerm")
LGB<-get.predict(df, "L", "G", "pigweed_num_nongerm", "B")
predict<-rbind(GB, LB, LG, LGB)
predict$pigweed_num_nongerm<-round(as.numeric(predict$pigweed_num_nongerm), 0)
predict$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict

# foxtail
LG<-get.predict(df, "L", "G", "foxtail_num_nongerm")
LB<-get.predict(df, "L", "B", "foxtail_num_nongerm")
GB<-get.predict(df, "G", "B", "foxtail_num_nongerm")
LGB<-get.predict(df, "L", "G", "foxtail_num_nongerm", "B")
predict1<-rbind(GB, LB, LG, LGB)
predict1$foxtail_num_nongerm<-round(as.numeric(predict1$foxtail_num_nongerm), 0)
predict1$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict1
predict<-full_join(predict, predict1)


# pigweed 2
LG<-get.predict(df, "L", "G", "pigweed_num_germ")
LB<-get.predict(df, "L", "B", "pigweed_num_germ")
GB<-get.predict(df, "G", "B", "pigweed_num_germ")
LGB<-get.predict(df, "L", "G", "pigweed_num_germ", "B")
predict2<-rbind(GB, LB, LG, LGB)
head(predict2)
predict2$pigweed_num_germ<-round(as.numeric(predict2$pigweed_num_germ), 0)
predict2$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict<-full_join(predict, predict2)
predict


# foxtail
LG<-get.predict(df, "L", "G", "foxtail_num_germ")
LB<-get.predict(df, "L", "B", "foxtail_num_germ")
GB<-get.predict(df, "G", "B", "foxtail_num_germ")
LGB<-get.predict(df, "L", "G", "foxtail_num_germ", "B")
predict3<-rbind(GB, LB, LG, LGB)
predict3$foxtail_num_germ<-round(as.numeric(predict3$foxtail_num_germ), 0)
predict3$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict3
predict<-full_join(predict, predict3)
head(predict)

## add to df
df1<-df %>% filter(n_species!="1")

df1<-full_join(df1, predict) 
df1<-df1 %>% ungroup()
as.factor(df1$Treatment)
head(df1)

# anova binomial model pigweed
#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

# anova binomial model foxtail
#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)


# there are difference between the predicted verse measured for weed seed decay for some treatments



# predicted vs measured for activity ######

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
# group by
df<- df %>% group_by(Rep, Nitrogen)
df

# pigdf
LG<-get.predict(df, "L", "G", "pigweed_prop_nongerm")
LB<-get.predict(df, "L", "B", "pigweed_prop_nongerm")
GB<-get.predict(df, "G", "B", "pigweed_prop_nongerm")
LGB<-get.predict(df, "L", "G", "pigweed_prop_nongerm", "B")
LG

predict<-rbind(GB, LB, LG, LGB)
predict$pigweed_prop_nongerm<-round(as.numeric(predict$pigweed_prop_nongerm), 3)
predict$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict


# foxtail
LG<-get.predict(df, "L", "G", "foxtail_prop_nongerm")
LB<-get.predict(df, "L", "B", "foxtail_prop_nongerm")
GB<-get.predict(df, "G", "B", "foxtail_prop_nongerm")
LGB<-get.predict(df, "L", "G", "foxtail_prop_nongerm", "B")
LG

predict1<-rbind(GB, LB, LG, LGB)
predict1$foxtail_prop_nongerm<-round(as.numeric(predict1$foxtail_prop_nongerm), 3)
predict1$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict1
predict<-full_join(predict, predict1)

## add to df
df1<-df %>% filter(n_species!="1")

df1<-full_join(df1, predict) 
df1<-df1 %>% ungroup()
as.factor(df1$Treatment)



















#### plot
mycols <- c( #IBM colors
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)
p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
p2
p3<-df1  %>% 
  ggplot(aes(x=Treatment, y=pigweed_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 16)+
  #geom_text(y=8, label =label , nudge_x = -.8, size=7)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))
p3
# put the two plots together
require(gridExtra)
grid.arrange(p2, p3, ncol=2)




