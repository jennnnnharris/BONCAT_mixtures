# Figure 3 prediction


# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
#library(lme4)
#library(nlme)
library(emmeans)
#library(multcomp)
#library(dplyr)


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

#plant functional traits without prediction 
#biomass####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read.csv("biomass_potlevel.csv") # biomass data

# process to make n species column, summarize at pot level, make treatment a factor.
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<- df %>% select(c(Brassicae, Legume, Grass )) %>% rowSums()
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))


head(df)

p1<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass.g, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=.7, width=.1)+
  
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=12))+
  labs(title = "A",
       x="",
       y="Root biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #
p1

p2<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Stem.Biomass.g, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=.7, width=.1)+
  
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0, size=12),
        legend.position = "none")+
  labs(title = "B",
       x="",
       y="shoot biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #
p2
require(gridExtra)
grid.arrange(p1, p2, ncol=2)


#root biomass N
#dfN0<- df %>% filter(N==0)
df1<-df %>% filter(n_species==1)
m1<-aov(Root.Biomass.g ~ Treatment, data = df1)
summary(m1) # difference between treatments
library(emmeans)
emm_object <- emmeans(m1, specs = ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 

summary(pairwise_comparison)
cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  Letters = letters) 
print(cld_result)
     


#shoot biomass
df1<-df %>% filter(n_species==1)

m1<-aov(Stem.Biomass.g ~ Treatment*N, data = df1)
summary(m1) # difference between treatments
emm_object <- emmeans(m1, specs = ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
library(multcomp)
summary(pairwise_comparison)
cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  Letters = letters) 
print(cld_result)

#G  L   B    GB  LB    LGB  LG
#bc cd  a    ab  bc    cd   d 





####### weed seed decay #######
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

# load data 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read.csv("weed_seed_decay.csv", row.names = 1 )
head(df)

# set factor
df$Treatment   <- factor(df$Treatment, levels= c("S", "L", "G", "B", "GB", "LB", "LG", "LGB"))
unique(df$Treatment)


# foxtail
#L G B LB GB LG LGB
# a a a ab ab ab b
df1<-df  %>%   filter(Treatment!="S") 
df1$Treatment   <- factor(df1$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

lab = as.character(df1$Treatment)
unique(lab)
lab<-gsub("LGB", "X", lab)
lab<-gsub("GB", "AX", lab)
lab<-gsub("LG", "AX", lab)
lab<-gsub("LB", "AX", lab)
lab<-gsub("B", "A", lab)
lab<-gsub("G", "A", lab)
lab<-gsub("L", "A", lab)

lab<-gsub("X", "B", lab)
unique(lab)
length(lab)


p1<- df1 %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_nongerm, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "none")+
  scale_fill_manual(values = IBM)+
  geom_text(y=.21, label = lab, size=3)+
  labs(title = "A",
       x="",
       y="non germinating foxtail (%)")
p1





# pigweed
# L   G   B   GB    LB  LG  LGB
# bc abc  d   ab    a   abc cd

lab = as.character(df1$Treatment)
unique(lab)
lab<-gsub("LGB", "CD", lab)
lab<-gsub("GB", "AX", lab)
lab<-gsub("LG", "AXC", lab)
lab<-gsub("LB", "A", lab)
lab<-gsub("B", "D", lab)
lab<-gsub("G", "ABC", lab)
lab<-gsub("L", "XC", lab)
lab<-gsub("X", "B", lab)

p2<-df1 %>%
  ggplot(aes(x=Treatment, y=pigweed_prop_nongerm, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "none")+
  scale_fill_manual(values = IBM)+
  geom_text(y=.73, label = lab, size=3)+
  labs(title = "B",
       x="",
       y="non germinating pigweed (%)")
p2


require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
svg("weeddecay.svg", width = 5, height = 2.5)
grid.arrange(p1, p2, ncol=2)
dev.off()




#### model just monocultures ####

# model - binomial model with percent data##
# make vector of successes and failures
df1<-df %>% filter(df$n_species==1) %>% filter(Treatment!="S")

# model foxtail
y<-cbind(df1$foxtail_num_germ, df1$foxtail_num_nongerm)
m1<-glm(data= df1, y~Treatment, family = binomial)
summary(m1) # no treatmen effect 
library(emmeans)
emm_object <- emmeans(m1, specs = ~ Treatment)
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
library(multcomp)
cld_result <- cld(emm_object, adjust = "tukey", alpha = 0.05, Letters = letters) 
print(cld_result)




# model pigweed
y<-cbind(df1$pigweed_num_germ, df1$pigweed_num_nongerm)
m1<-glm(data= df1, y~Treatment, family = binomial)
summary(m1) # no treatmen effect 
library(emmeans)
emm_object <- emmeans(m1, specs = ~ Treatment)
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
library(multcomp)
cld_result <- cld(emm_object, adjust = "tukey", alpha = 0.05, Letters = letters) 
print(cld_result)




############## N fix figures ####################

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df.leg<- read.csv("Nfix.csv")
head(df.leg) 

# by treatment 
label <- df.leg$Treatment
label
label <- gsub("LGB", "C" ,label )
label<- gsub("LG", "B" ,label )
label <- gsub("LB", "B" ,label )
label<- gsub("L", "A" ,label )


p1<- ggplot(df.leg, aes(x=Treatment, y=perc.Ndfa, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = legume_cols)+
  # geom_text(y=92, label = label, size=4)+
  labs(title = "E",
       x="",
       y= "Nitrogen from Fixation (%)") +
  scale_shape_manual(values = c(17, 16)) +
  facet_grid(~Nitrogen_label)#


p1

# Analysis of variance 
# full model 
one.way.Nadd <- aov(n_fix_per_legume ~ Treatment*Nitrogen, data = df.leg)
summary(one.way.Nadd) # difference between treatments

df<-df.leg %>% filter(Nitrogen==0)
one.way.Nadd <- aov(n_fix_per_legume ~ Treatment, data = df)
tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
print(tukey.result.Nadd) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine

# by treatment 
label <- df.leg$Treatment
label

label <- gsub("LGB", "B" ,label )
label<- gsub("LG", "AB" ,label )
label <- gsub("LB", "AB" ,label )
label<- gsub("L", "A" ,label )

p2<- ggplot(df.leg, aes(x=Treatment, y=n_fix_per_legume, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = legume_cols)+
  # geom_text(y=.07, label = label, size=4)+
  labs(title = "F",
       x="",
       y= "N fixed (mg per legume)")  +
  scale_shape_manual(values = c(17, 16)) +
  facet_grid(~Nitrogen_label)#

p2


# Analysis of variance 
df<-df.leg %>% filter(Nitrogen==1)
one.way.Nadd <- aov(n_fix_per_legume ~ Treatment, data = df)
summary(one.way.Nadd) # difference between treatments
tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
print(tukey.result.Nadd) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine\

df<-df.leg %>% filter(Nitrogen==0)
one.way.Nadd <- aov(n_fix_per_legume ~ Treatment, data = df)
summary(one.way.Nadd) # difference between treatments
tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
print(tukey.result.Nadd) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
svg("nfix.svg", height = 3, width = 8)
grid.arrange(p1, p2, ncol=2)
dev.off()


# nfixed per mg of soil
ggplot(df.leg, aes(x=treatment, y=totalN.mg.g.1, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  ylab('Total Nitrogen Fixed per gram of soil') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = mycols4)+
  facet_grid( ~spp.number, scales = "free", space = "free")







# Analysis of variance 
# full model 
one.way.Nadd <- aov(perc.Ndfa ~ Treatment*Nitrogen, data = df.leg)
summary(one.way.Nadd) # difference between treatments

df<-df.leg %>% filter(Nitrogen==0)
one.way.Nadd <- aov(perc.Ndfa ~ Treatment, data = df)
print(tukey.result.Nadd) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine


df<-df.leg %>% filter(Nitrogen==1)
one.way.Nadd <- aov(perc.Ndfa ~ Treatment, data = df)
tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
print(tukey.result.Nadd) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine


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
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB", "GB.predict", "LB", "LB.predict",  "LG", 
                                 "LG.predict", "LGB", "LGB.predict" ))

#### plot predictions from monocultures for biomass

mycols <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
    "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)


label <- df1$Treatment
label <- gsub("LGB.predict", "*" ,label )
label <- gsub("LB.predict", "*" ,label )
label <- gsub("LG.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label

p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
       )+
  
  geom_text(y=6, label =label , nudge_x = -.5, size=8)+
  labs(title = "A",
       x="",
       y="Root biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #



p1



p2<-df1  %>% 
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

p2

# put the plots together

require(gridExtra)
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
#svg("biomasspredict.svg", height = 6, width = 4.5)
grid.arrange(p1, p2, ncol=1)
#dev.off()


######### anova#####

# root biomass
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
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
 #df1<-df %>% filter(n_species!="1")

df1<-full_join(df, predict) 
df1<-df1 %>% ungroup()
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB", "GB.predict", "LB", "LB.predict",  "LG", 
                                                "LG.predict", "LGB", "LGB.predict" ))

#### plot predictions from monocultures for weed seed

mycols <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)



# make labels 
label <- df1$Treatment
label
label <- gsub("LGB.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label

# plot
p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none")+
  labs(title = "C",
       x="",
       y="non germinating foxtail (%)")+
  geom_text(y=.25, label =label , nudge_x = -.8, size=8)
  
  
p1

# make labels 
label <- df1$Treatment
label <- gsub("LB.predict", "*" ,label )
label <- gsub("GB.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label


# plot
p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=pigweed_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  geom_text(y=.8, label =label , nudge_x = -.8, size=8)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",)+
  labs(title = "D",
       x="",
       y="non germinating pigweed (%)")

p2

# put the two plots together
require(gridExtra)
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
#svg("weeddecay.predict.svg", width = 4.5, height = 6)
grid.arrange(p1, p2, ncol=1)
#dev.off()

##########stats weed prop non germinated #####
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

