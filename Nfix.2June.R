### June 2nd, 2025
### Emma Rice
### Running r version 4.3.2


### Purpose: Calculating percent nitrogen from fixation using the naturl abundance method
#            Looking at samples with nitrogen added to pots and not added
#            Calcuating %Ndfa in legumes under each treatment
#            And how much 15N is shared with grass and brassica when in mixtures with legumes under different treatments


# clear workspace
rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

### Set the working directory
#setwd("/Users/emmarice/Documents/PhD Research/Greenhouse Exp 2022/FixedN")
path <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Analysis/N Fixation from Emma"
setwd(path)

# set colors
#import colors
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")
mycols7<-c( "#4B2D4BFF", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
mycols4 <- c("#4B2D4BFF",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "LB", "LG", "LGB"))
fig2path <-  "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig2_nspecies"
fig3path <-  "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig3_predict"

### Import Field Mixture Data 
df <- read.csv("merged.plate.ghbiomass.sheet.csv", header=T, stringsAsFactors = F)

### Load necessary packges 
library(dplyr)
library(ggplot2)
library(tidyr)
library(nlme)

#### Getggplot2#### Get set up for calculations ####
### Take average of δ15N triticale monoculture across all 6 pots for δ15Nref (standard background δ15N reference) 
# Select for only monoculture and triticale
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

#### Calculate %Ndfa with for legumes only ####
## N addition first
dfNadd.leg <- dfNadd[dfNadd$spp.type %in% "legume", ]
# Equation: %Ndfa = 100[(δ15Nref - δ15Nleg) / (δ15Nref - B)]
dfNadd.leg$perc.Ndfa <- 100*((δ15Nref - dfNadd.leg$δ15Nvs.At.Air) / (δ15Nref - dfNadd.leg$B))

# Calculate the average for each treatment
dfNadd.leg.trtavg <- dfNadd.leg %>%
  group_by(treatment) %>%
  summarise(avg.perc.Ndfa = mean(perc.Ndfa))

# Calculating how much of total plant nitrogen from BNF in g/pot (need to update this with soil volume info)
# BNF (g pot-1) = biomass (g plot-1) • plant N concentration (%)/100 • %Ndfa/100
dfNadd.leg$BNF.g.plot.1 <- dfNadd.leg$stem.biomass.g * (dfNadd.leg$perc.N/100) * (dfNadd.leg$perc.Ndfa/100)
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

# Calculating how much of total plant nitrogen from BNF in g/pot (need to update this with soil volume info)
# BNF (g pot-1) = biomass (g plot-1) • plant N concentration (%)/100 • %Ndfa/100
dfnoNadd.leg$BNF.g.plot.1 <- dfnoNadd.leg$stem.biomass.g * (dfnoNadd.leg$perc.N/100) * (dfnoNadd.leg$perc.Ndfa/100)
# How much total N in each pot in g
dfnoNadd.leg$totalN.g.pot.1 <- NA  
dfnoNadd.leg$totalN.g.pot.1 <- dfnoNadd.leg$stem.biomass.g * (dfnoNadd.leg$perc.N/100) 
# Soil N retention in g.pot
dfnoNadd.leg$soilNret.g.pot.1 <- NA  
dfnoNadd.leg$soilNret.g.pot.1 <- dfnoNadd.leg$totalN.g.pot.1 - dfnoNadd.leg$BNF.g.plot.1

#### Calculate biomass N in mg N /g soil ####
# No added 
# converting from g/pot to mg N / g soil
# soil mass in g = 3 liters soil * 1.3 g/mL (approximate soil density) * 1000
dfNadd.leg$BNF.mg.g.1 <- NA  
dfNadd.leg$BNF.mg.g.1 <- (dfNadd.leg$BNF.g.plot.1*1000)/3900
dfNadd.leg$soilNret.mg.g.1 <- NA  
dfNadd.leg$soilNret.mg.g.1 <- (dfNadd.leg$soilNret.g.pot.1*1000)/3900
dfNadd.leg$totalN.mg.g.1 <- NA  
dfNadd.leg$totalN.mg.g.1 <- (dfNadd.leg$totalN.g.pot.1*1000)/3900

# No N
dfnoNadd.leg$BNF.mg.g.1 <- NA  
dfnoNadd.leg$BNF.mg.g.1 <- (dfnoNadd.leg$BNF.g.plot.1*1000)/3900
dfnoNadd.leg$soilNret.mg.g.1 <- NA  
dfnoNadd.leg$soilNret.mg.g.1 <- (dfnoNadd.leg$soilNret.g.pot.1*1000)/3900
dfnoNadd.leg$totalN.mg.g.1 <- NA  
dfnoNadd.leg$totalN.mg.g.1 <- (dfnoNadd.leg$totalN.g.pot.1*1000)/3900



#### Plots - legume N added  ####
dfNadd.leg$treatment <- factor(dfNadd.leg$treatment)

########jenny plot by treatment ####### 

#%Ndfa boxplots
setwd(fig2path)


# put together N+ and N-
dfNadd.leg
dfnoNadd.leg
df.leg<-rbind(dfnoNadd.leg, dfNadd.leg)

df.leg$treatment <- factor(df.leg$treatment)


svg(file="nfix.treatment.svg",width = 2.3, height=2.3)
ggplot(df.leg, aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  ylab('Nitrogen from Fixation (%)') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = mycols4)+
  facet_grid( ~spp.number, scales = "free", space = "free")
  
dev.off() 

# Analysis of variance 
one.way.Nadd <- aov(perc.Ndfa ~ treatment, data = dfNadd.leg)
summary(one.way.Nadd) # difference between treatments
tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
print(tukey.result.Nadd) # All difference except LG-LB

plot(one.way.Nadd) #homoscedasticity looks fine

### jenny plot by n species ###############
dfNadd.leg$treatment <- factor(dfNadd.leg$treatment)

# put together N+ and N-
dfNadd.leg
dfnoNadd.leg
df.leg<-rbind(dfnoNadd.leg, dfNadd.leg)


# plot with colors
ggplot(df.leg, aes(x=spp.number, y=perc.Ndfa, colour = treatment)) + 
  geom_jitter(size=2) +
  ylab('Nitrogen from Fixation (%)') +
  xlab("Number of species ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  scale_color_manual(values=mycols4)

# put together N+ and N-
dfNadd.leg
dfnoNadd.leg
df.leg<-rbind(dfnoNadd.leg, dfNadd.leg)


setwd(fig2path)
svg(file="nfix.species.svg",width = 2.4, height=2.4)
ggplot(df.leg, aes(x=spp.number, y=perc.Ndfa)) + 
  ylab('Nitrogen from Fixation (%)') +
  xlab("Number of species") +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= "grey10")+
  theme_classic(base_size = 12)
  #theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
  #      plot.title = element_text(hjust = 0.5))

dev.off()


svg(file="col.nfix.species.svg",width = 3, height=3)
ggplot(dfNadd.leg, aes(x=spp.number, y=perc.Ndfa, colour = treatment)) + 
  ylab('Nitrogen from Fixation (%)') +
  xlab("Number of species") +
  geom_jitter(width = .2, size=2 )+
  geom_smooth(method = lm, color= "grey10")+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=mycols4)


dev.off()

# run linear model

# linear model
m1<-lm( perc.Ndfa~spp.number , data= dfNadd.leg )
summary(m1)
plot(m1)

#end#

# Looking at %Ndfa with legume proportion in pot
# Fit the model
model.Nadd<-lm(perc.Ndfa ~ prop.legume + treatment, data=dfNadd.leg)
vif(model.Nadd, type="predictor") # VIF high, leave for now because probably wont use this, look at points without line of fit

# Meet linear assumptions?
plot(model.Nadd)
E<-resid(model.Nadd)
hist(E, xlab="residuals", main="") # okay
#test
shapiro.test(resid(model.Nadd)) # no, not normally distributed, 0.003661


# Create a new data frame for predictions
new.df.Nadd <- expand.grid(
  prop.legume = seq(min(dfNadd.leg$prop.legume), max(dfNadd.leg$prop.legume), length.out = 100),
  treatment = unique(dfNadd.leg$treatment)
)

# Predict using the model
new.df.Nadd$predicted.Ndfa <- predict(model.Nadd, newdata = new.df.Nadd)

# visualize
ggplot(dfNadd.leg, aes(x = prop.legume, y = perc.Ndfa, fill = treatment)) + 
#  geom_line(data = new.df.Nadd, aes(x = prop.legume, y = predicted.Ndfa, color = treatment, group = treatment), size = 1, show.legend = FALSE) +
#  scale_color_manual(values = c("L" = "#018571", 
#                                "LB" = "#80cdc1", 
#                                "LG" = "#dfc27d", 
#                                "LGB" = "#a6611a")) +
  geom_point(size = 5, shape = 21, color = "black") + 
  scale_fill_manual(values = c("L" = "#018571", 
                               "LB" = "#80cdc1", 
                               "LG" = "#dfc27d", 
                               "LGB" = "#a6611a")) +
  ylab("%Ndfa") + 
  xlab("Proportion Legume") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 12)) 
# Shows that the less legume prop in mixture the more N derived from fixation
# less clear than noN addition

## biomass N from fixation (mg N/g soil)
ggplot(dfNadd.leg, aes(x=treatment, y=BNF.mg.g.1, fill=treatment)) + 
  geom_boxplot() +
  ylab('Nitrogen from Fixation \n(mg N/g soil)') +
  xlab("Legume Treatments") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 12), 
        legend.position = "none") + 
  scale_fill_manual(values = c("L" = "#018571", 
                               "LB" = "#80cdc1", 
                               "LG" = "#dfc27d", 
                               "LGB" = "#a6611a")
  )

# Analysis of variance 
one.way.bnf.Nadd <- aov(BNF.mg.g.1 ~ treatment, data = dfNadd.leg)
summary(one.way.bnf.Nadd) # difference between treatments
tukey.result.bnf.Nadd <- TukeyHSD(one.way.bnf.Nadd)
print(tukey.result.bnf.Nadd) # all difference except LG-LB and LGB-LG

plot(one.way.bnf.Nadd) #homoscedasticity looks fine


# Plot with shoot biomass N with bars divided into Fixed N and Soil N sources
df.long.N <- subset(dfNadd.leg, select = c(treatment, BNF.mg.g.1, soilNret.mg.g.1))

df.long1.N <- df.long.N %>%
  gather(key = "variable", value = "value", BNF.mg.g.1, soilNret.mg.g.1)

# side by side bars
# Summarize data: Calculate the mean across replicates for each treatment and variable
df.summary.N <- df.long1.N %>%
  dplyr::group_by(treatment, variable) %>%
  dplyr::summarize(
    mean.value = mean(value, na.rm = TRUE),
    se.value = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )

# Visual
ggplot(df.summary.N, aes(x = treatment, y = mean.value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black") +
  geom_errorbar(aes(ymin = mean.value - se.value, ymax = mean.value + se.value), 
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(
    values = c("BNF.mg.g.1" = "#80cdc1", "soilNret.mg.g.1" = "#dfc27d"),
    labels = c("Fixed N", "Soil N uptake")
  ) +
  labs(
    x = "Treatment",  
    y = "Aboveground biomass Nitrogen \n(mg N / g soil)", 
    fill = "Nitrogen source"     
  ) + 
  theme_bw() +
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.text = element_text(size = 12, colour = "black")) + 
  theme(axis.title = element_text(size = 12)) + 
  theme(legend.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 10)) 
  # In mixes less soil N for legumes because likely taken up by other spp




#### Plots - legume no N added  ####
dfnoNadd.leg$treatment <- factor(dfnoNadd.leg$treatment)
# %Ndfa boxplots
ggplot(dfnoNadd.leg, aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
  geom_boxplot() +
  ylab('Nitrogen from Fixation (%)') +
  xlab("Legume Treatments") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c("L" = "#018571", 
                               "LB" = "#80cdc1", 
                               "LG" = "#dfc27d", 
                               "LGB" = "#a6611a")
  )

  # Analysis of variance 
  one.way <- aov(perc.Ndfa ~ treatment, data = dfnoNadd.leg)
  summary(one.way) # difference between treatments
  tukey.result <- TukeyHSD(one.way)
  print(tukey.result)
  
  plot(two.way) #homoscedasticity looks fine
  
  
# Looking at %Ndfa with legume proportion in pot
  # Fit the model
  model<-lm(perc.Ndfa ~ prop.legume + treatment, data=dfnoNadd.leg)
  vif(model, type="predictor")
  
  # Meet linear assumptions?
  plot(model)
  E<-resid(model)
  hist(E, xlab="residuals", main="") # not great
  #test
  shapiro.test(resid(model)) # yes, normally distributed! 0.126
  
  
  # Create a new data frame for predictions
  new.df <- expand.grid(
    prop.legume = seq(min(dfnoNadd.leg$prop.legume), max(dfnoNadd.leg$prop.legume), length.out = 100),
    treatment = unique(dfnoNadd.leg$treatment)
  )
  
  
  # Predict using the model
  new.df$predicted.Ndfa <- predict(model, newdata = new.df)
  
  # visualize
  ggplot(dfnoNadd.leg, aes(x = prop.legume, y = perc.Ndfa, fill = treatment)) + 
  #  geom_line(data = new.df, aes(x = prop.legume, y = predicted.Ndfa, color = treatment, group = treatment), size = 1, show.legend = FALSE) +
  #  scale_color_manual(values = c("L" = "#018571", 
  #                               "LB" = "#80cdc1", 
  #                               "LG" = "#dfc27d", 
  #                               "LGB" = "#a6611a")) +
    geom_point(size = 5, shape = 21, color = "black") + 
    scale_fill_manual(values = c("L" = "#018571", 
                                 "LB" = "#80cdc1", 
                                 "LG" = "#dfc27d", 
                                 "LGB" = "#a6611a")) +
    ylab("%Ndfa") + 
    xlab("Proportion Legume") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 11, colour = "black"), 
          axis.title = element_text(size = 12)) 
  # Shows that the less legume prop in mixture the more N derived from fixation
  
## biomass N from fixation (kg ha-1)
  ggplot(dfnoNadd.leg, aes(x=treatment, y=BNF.mg.g.1, fill=treatment)) + 
    geom_boxplot() +
    ylab('Nitrogen from Fixation \n(mg N/g soil)') +
    xlab("Legume Treatments") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 11, colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(size = 12), 
          legend.position = "none") + 
     scale_fill_manual(values = c("L" = "#018571", 
                               "LB" = "#80cdc1", 
                               "LG" = "#dfc27d", 
                               "LGB" = "#a6611a")
     )
  
# Analysis of variance 
one.way.bnf <- aov(BNF.mg.g.1 ~ treatment, data = dfnoNadd.leg)
summary(one.way.bnf) # difference between treatments
tukey.result.bnf <- TukeyHSD(one.way.bnf)
print(tukey.result.bnf)
  
plot(one.way.bnf) #homoscedasticity looks fine
  
  
# Plot with shoot biomass N with bars divided into Fixed N and Soil N sources
df.long.noN <- subset(dfnoNadd.leg, select = c(treatment, BNF.mg.g.1, soilNret.mg.g.1))

df.long1.noN <- df.long.noN %>%
  gather(key = "variable", value = "value", BNF.mg.g.1, soilNret.mg.g.1)

# side by side bars
# Summarize data: Calculate the mean across replicates for each treatment and variable
df.summary.N <- df.long1.noN %>%
  dplyr::group_by(treatment, variable) %>%
  dplyr::summarize(
    mean.value = mean(value, na.rm = TRUE),
    se.value = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )

# Visual
ggplot(df.summary.N, aes(x = treatment, y = mean.value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black") +
  geom_errorbar(aes(ymin = mean.value - se.value, ymax = mean.value + se.value), 
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(
    values = c("BNF.mg.g.1" = "#80cdc1", "soilNret.mg.g.1" = "#dfc27d"),
    labels = c("Fixed N", "Soil N uptake")
    ) +
  labs(
    x = "Treatment",  
    y = "Aboveground biomass N \n(mg N/g soil)", 
    fill = "Nitrogen source"     
  ) + 
  theme_bw() +
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.text = element_text(size = 12, colour = "black")) + 
  theme(axis.title = element_text(size = 12)) + 
  theme(legend.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 10)) 





#### Calculate shared fixed N with canola and triticale when in mixture ####
## Triticale first
## N added
dfNadd.g <- dfNadd[dfNadd$spp.type %in% "grass", ]
dfNadd.gmix <- dfNadd.g[!dfNadd.g$spp.number %in% "1", ] # Only want to look in mixes

# Equation: % trit N from transfer = 100[(δ15Ntrit monoculture – δ15Ntrit mix)/(δ15Ntrit monoculture – C)
# C is the lowest δ15N grass in mixture or 0, whichever is lowest, here C = 0 lowest
dfNadd.gmix$perc.Ndfa <- 100*((δ15Nref - dfNadd.gmix$δ15Nvs.At.Air)/(δ15Nref - 0))

## No N added
dfnoNadd.g <- dfnoNadd[dfnoNadd$spp.type %in% "grass", ]
dfnoNadd.gmix <- dfnoNadd.g[!dfnoNadd.g$spp.number %in% "1", ] # Only want to look in mixes

# Equation: % trit N from transfer = 100[(δ15Ntrit monoculture – δ15Ntrit mix)/(δ15Ntrit monoculture – C)
# C is the lowest δ15N grass in mixture or 0, whichever is lowest, here C = 0 lowest
dfnoNadd.gmix$perc.Ndfa <- 100*((δ15Nref.noNadd - dfnoNadd.gmix$δ15Nvs.At.Air)/(δ15Nref.noNadd - 0))

# In both instances (N added and no N added) the %Ndfa is very low of negative which is impossible. 
#   Not sure I trust this measurement, we could make all negatives 0, but I would lean toward not including it at all.


## Canola
## N added
### Take average of δ15N canola monoculture across all 6 pots for δ15Nref.b
# Select for only monoculture and canola
bnref <- df[df$spp.number %in% "1", ] 
bnref1 <- bnref[bnref$treatment %in% "B", ]
bnrefN <- bnref1[bnref1$nitrogen.added %in% "Y", ]
# Where duplicate take average so one reference not over weighted
bnrefN.avg <- bnrefN %>%
  group_by(pot.number, rep) %>%
  summarise(avg.δ15N = mean(δ15Nvs.At.Air))
# Take average across all reps
δ15Nref.b <- mean(bnrefN.avg$avg.δ15N)

# Get mixture df set up for canola
dfNadd.b <- dfNadd[dfNadd$spp.type %in% "brassica", ]
dfNadd.bmix <- dfNadd.b[!dfNadd.b$spp.number %in% "1", ] # Only want to look in mixes

# Equation: % trit N from transfer = 100[(δ15Ntrit monoculture – δ15Ntrit mix)/(δ15Ntrit monoculture – C)
# C is the lowest δ15N grass in mixture or 0, whichever is lowest, here C = 0 lowest
dfNadd.bmix$perc.Ndfa <- 100*((δ15Nref.b - dfNadd.bmix$δ15Nvs.At.Air)/(δ15Nref.b - 0))

## No N added
### Take average of δ15N canola monoculture across all 6 pots for δ15Nref.b
# Select for only monoculture and canola
bnrefnoN <- bnref1[bnref1$nitrogen.added %in% "N", ]
# Where duplicate take average so one reference not over weighted
bnrefN.avg.noN <- bnrefnoN %>%
  group_by(pot.number, rep) %>%
  summarise(avg.δ15N = mean(δ15Nvs.At.Air))
# Take average across all reps
δ15Nref.b.noNadd <- mean(bnrefN.avg.noN$avg.δ15N)

# Get mixture df set up for canola
dfnoNadd.b <- dfnoNadd[dfnoNadd$spp.type %in% "brassica", ]
dfnoNadd.bmix <- dfnoNadd.b[!dfnoNadd.b$spp.number %in% "1", ] # Only want to look in mixes

# Equation: % canola N from transfer = 100[(δ15Ncanola monoculture – δ15Ncanola mix)/(δ15Ncanola monoculture – C)
# C is the lowest δ15N grass in mixture or 0, whichever is lowest, here C = 0 lowest
dfnoNadd.bmix$perc.Ndfa <- 100*((δ15Nref.b.noNadd - dfnoNadd.bmix$δ15Nvs.At.Air)/(δ15Nref.b.noNadd - 0))

# again mostly negative

