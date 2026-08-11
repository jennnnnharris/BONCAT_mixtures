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



# calculate nitrogen from fixation #############


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Nfix")
df <- read.csv("raw_nitrogen_fixation.csv", header=T, stringsAsFactors = F) # fix data


### Take average of δ15N triticale monoculture across all 6 pots for 
# δ15Nref (standard background δ15N reference) 
# Select for only monoculture and triticale
# reference = Tricale monculture with nitrogen, which should have 
# pretty much zero nitrogen from fixation. 
tnref <- df[df$spp.number %in% "1", ] 
tnref1 <- tnref[tnref$treatment %in% "G", ]
tnrefN <- tnref1[tnref1$nitrogen.added %in% "Y", ]
# Where duplicate take average so one reference not over weighted
tnrefN.avg <- tnrefN %>%
  group_by(pot.number, rep) %>%
  summarise(avg.δ15N = mean(δ15Nvs.At.Air))
# Take average across all reps
δ15Nref <- mean(tnrefN.avg$avg.δ15N)

## Repeat with noN addition samples
tnrefnoN <- tnref1[tnref1$nitrogen.added %in% "N", ]
# Where duplicate take average so one reference not over weighted
tnrefnoN.avg <- tnrefnoN %>%
  group_by(pot.number, rep) %>%
  summarise(avg.δ15N = mean(δ15Nvs.At.Air))
# Take average across all reps
δ15Nref.noNadd <- mean(tnrefnoN.avg$avg.δ15N)

## Divide the df into nitrogen added and no nitrogen added
dfNadd <- df[df$nitrogen.added %in% "Y", ]
dfnoNadd <- df[df$nitrogen.added %in% "N", ]

#Calculate %Ndfa with for legumes only ##
## N addition first
dfNadd.leg <- dfNadd[dfNadd$spp.type %in% "legume", ]
# Equation: %Ndfa = 100[(δ15Nref - δ15Nleg) / (δ15Nref - B)]
dfNadd.leg$perc.Ndfa <- 100*((δ15Nref - dfNadd.leg$δ15Nvs.At.Air) / (δ15Nref - dfNadd.leg$B))


# Calculate the average for each treatment
dfNadd.leg.trtavg <- dfNadd.leg %>%
  group_by(treatment) %>%
  summarise(avg.perc.Ndfa = mean(perc.Ndfa))




# Calculating how much of total plant nitrogen from BNF in g/pot 

# BNF (g pot-1) = biomass (g plot-1) • plant N concentration (%)/100 • %Ndfa/100
dfNadd.leg$BNF.g.pot.1 <- dfNadd.leg$stem.biomass.g * (dfNadd.leg$perc.N/100) * (dfNadd.leg$perc.Ndfa/100)

# How much total N in each pot in g
dfNadd.leg$totalN.g.pot.1 <- NA  
dfNadd.leg$totalN.g.pot.1 <- dfNadd.leg$stem.biomass.g * (dfNadd.leg$perc.N/100) 
# Soil N retention in g.pot
dfNadd.leg$soilNret.g.pot.1 <- NA  
dfNadd.leg$soilNret.g.pot.1 <- dfNadd.leg$totalN.g.pot.1 - dfNadd.leg$BNF.g.plot.1

## No N addition 
dfnoNadd.leg <- dfnoNadd[dfnoNadd$spp.type %in% "legume", ]
# Equation: %Ndfa = 100[(δ15Nref - δ15Nleg) / (δ15Nref - B)]
dfnoNadd.leg$perc.Ndfa <- 100*((δ15Nref - dfnoNadd.leg$δ15Nvs.At.Air) / (δ15Nref - dfnoNadd.leg$B))

# Calculate the average for each treatment
dfnoNadd.leg.trt.avg <- dfnoNadd.leg %>%
  group_by(treatment) %>%
  summarise(avg.perc.Ndfa = mean(perc.Ndfa))
# average %Ndfa the same in LGB with and without N addition, all others higher without N added

# Calculating how much of total plant nitrogen from BNF in g/pot 
# BNF (g pot-1) = biomass (g plot-1) • plant N concentration (%)/100 • %Ndfa/100
dfnoNadd.leg$BNF.g.pot.1 <- dfnoNadd.leg$stem.biomass.g * (dfnoNadd.leg$perc.N/100) * (dfnoNadd.leg$perc.Ndfa/100)
# How much total N in each pot in g
dfnoNadd.leg$totalN.g.pot.1 <- NA  
dfnoNadd.leg$totalN.g.pot.1 <- dfnoNadd.leg$stem.biomass.g * (dfnoNadd.leg$perc.N/100) 
# Soil N retention in g.pot
dfnoNadd.leg$soilNret.g.pot.1 <- NA  
dfnoNadd.leg$soilNret.g.pot.1 <- dfnoNadd.leg$totalN.g.pot.1 - dfnoNadd.leg$BNF.g.plot.1


# now we have a few data dataframes
head(dfnoNadd.leg)
head(dfNadd.leg)

#put together N+ and N- ##
dfNadd.leg
dfnoNadd.leg
df.leg<-rbind(dfnoNadd.leg, dfNadd.leg)
#df.leg$treatment <- factor(df.leg$treatment)

# calcluated the amount fo N fixed per plant
df.leg<-df.leg%>%
  mutate(n_legume = recode(spp.number,
                           "1" = 6,
                           "2" = 3,
                           "3" =2)) %>%
  mutate(totalN.g.pot.1.n_legume = totalN.g.pot.1/n_legume,
         n_fix_per_legume = BNF.g.pot.1/n_legume )

df.leg

head(df.leg)
#optional: write df of process data
#write.csv(df.leg, "Nfix.csv")





############## percent nitrogen from BNF ######


head(df.leg) 


# add nitrogen label
  df.leg<-df.leg%>%
    mutate(Nitrogen_label = recode(nitrogen.added,
                             "N" = "Nitrogen -",
                             "Y" = "Nitrogen +")) 
  
  

# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_blockinfo.xlsx")
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
 
 

