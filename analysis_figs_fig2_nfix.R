# N fix 
# clear workspace and restart R
rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)


# load colors
legume_cols <- c( #IBM colors
  "navy", # dark royal blue L
  "#FE6100", # bright orange LB
  "#FFB000", # golden yellow LG
  "#865338" # medium mocha brown LGB
)



############## percent nitrogen from BNF ######

#load data
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df.leg<- read.csv("Nfix.csv")
head(df.leg) 


# add nitrogen label
  df.leg<-df.leg%>%
    mutate(Nitrogen_label = recode(nitrogen.added,
                             "N" = "Nitrogen -",
                             "Y" = "Nitrogen +")) 
  
  

# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_experiment_planning.xlsx")
head(block)  
df.leg<-left_join(df.leg, block)
  
# plot
p1<- ggplot(df.leg, aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = legume_cols)+
   labs(#title = "E",
       x="",
       y= "Nitrogen from Fixation (%)") +
  scale_shape_manual(values = c(17, 16)) +
  facet_grid(~Nitrogen_label)
p1


####### stats percent N from BNF ############

# load libraries 
library(lme4)
library(nlme)
library(multcomp)
library(emmeans)

# overall model perc ndfa
m1<- aov(perc.Ndfa ~ treatment*nitrogen.added+block, data = df.leg)
summary(m1)

# 2. Get estimated marginal means grouped by Nitrogen level
emm_object <- emmeans(m1, ~ treatment | nitrogen.added)

# 3. Perform pairwise comparisons within each Nitrogen level
pairwise_comparison <- pairs(emm_object, adjust = "tukey")
summary(pairwise_comparison)

# Requires multcomp / multcompView packages
cld_results <- cld(emm_object, Letters = letters, adjust = "tukey")
print(cld_results)


########################n fix per legume #######

p2<- ggplot(df.leg, aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = legume_cols)+
   labs(#title = "F",
       x="",
       y= "N fixed (g per legume)")  +
  scale_shape_manual(values = c(17, 16)) +
  facet_grid(~Nitrogen_label)#

p2
 
# Analysis of variance 
# overall model 
m1<- aov(n_fix_per_legume ~ treatment*nitrogen.added+block, data = df.leg)
summary(m1)

# 2. Get estimated marginal means grouped by Nitrogen level
emm_object <- emmeans(m1, ~ treatment | nitrogen.added)

# 3. Perform pairwise comparisons within each Nitrogen level
pairwise_comparison <- pairs(emm_object, adjust = "tukey")
summary(pairwise_comparison)

# Requires multcomp / multcompView packages
cld_results <- cld(emm_object, Letters = letters, adjust = "tukey")
print(cld_results)




###plot ###
 setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
 svg("nfix.svg", height = 4.5, width = 3.5)
  require(gridExtra)
  grid.arrange(p1, p2, ncol=1)
 dev.off()
 
 

