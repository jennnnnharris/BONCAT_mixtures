# N fix and biomass weed seed
# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)

#old cols:


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



############## N fix figures ####################

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


# Analysis of variance percent N from BNF
library(multcomp)
library(emmeans)
# overall model 
m1<- lm(perc.Ndfa ~ treatment*nitrogen.added+block, data = df.leg)
summary(m1)
anov = aov(m1)
summary(anov)

# Nitrogen -
df<-df.leg %>% filter(nitrogen.added=="N")
m1 <- aov(perc.Ndfa ~ treatment, data = df)
summary(m1) # difference between treatments
tukey <- TukeyHSD(m1)
print(tukey) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine
# emmeans 
emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

library(multcompView)
cld <- multcompLetters4(m1, tukey)
print(cld)
#LGB  LG  LB   L 
#"a" "a" "a" "b"

df<-df.leg %>% filter(nitrogen.added=="Y")
m1 <- aov(perc.Ndfa ~ treatment, data = df)
summary(m1) # difference between treatments
tukey <- TukeyHSD(m1)
print(tukey) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine
# emmeans 
emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

library(multcompView)
cld <- multcompLetters4(m1, tukey)
print(cld)

########################################n fix per legume #

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
 
library(multcomp)
library(emmeans)
# Analysis of variance 
#overall 
m1<- lm(n_fix_per_legume ~ treatment*nitrogen.added+block, data = df.leg)
summary(m1)
anov = aov(m1)
summary(anov)

# nitrogen -
df<-df.leg %>% filter(nitrogen.added=="N")
m1 <- aov(n_fix_per_legume ~ treatment+block, data = df)
summary(m1) # difference between treatments

emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  Letters = letters) 
print(cld_result)







# with nitrogen +
df<-df.leg %>% filter(nitrogen.added=="Y")
m1 <- aov(n_fix_per_legume ~ treatment+block, data = df)
summary(m1) # difference between treatments

emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)


cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  Letters = letters) 
print(cld_result)





 setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
 svg("nfix.svg", height = 4.5, width = 3.5)
  require(gridExtra)
  grid.arrange(p1, p2, ncol=1)
 dev.off()
 
 

