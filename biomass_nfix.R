# Figure 2 N fix and biomass 


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



#### import biomass data and process #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read.csv("Rice_greenhouse_ccexp_biomass_block.csv") # biomass data

# process to make n species column, summarize at pot level, make treatment a factor.
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<-rowSums(select(df, Brassicae, Legume, Grass))
df<-df %>% group_by(Pot.Number, Treatment, Rep, Brassicae, N, Legume, Grass, n_species, Block ) %>% summarise(Root.Biomass = sum(Total.Root.g), Shoot.Biomass = sum(Stem.Biomass.g), )
df$Root.to.Shoot <- df$Root.Biomass / df$Shoot.Biomass
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# give long composition name
# make composition var
df$composition <- df$Treatment
df$composition<-gsub("L", "Legume", df$composition)
df$composition<-gsub("G", "Grass", df$composition)
df$composition<-gsub("B", "Brassica", df$composition)
df$composition<-gsub("LegumeBrassica", "Legume_Brassica", df$composition)
df$composition<-gsub("LegumeGrass", "Legume_Grass", df$composition)
df$composition<-gsub("GrassBrassica", "Grass_Brassica", df$composition)
df$composition 
df$composition<- factor(df$composition, levels = c("Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))

# rename
dfb <- df
head(dfb)

###biomass figures ##########

# n species###
p1<-df  %>%
  ggplot(aes(x=n_species, y=Shoot.Biomass, colour = Treatment )) +
  geom_jitter(width = .1, size=1 )+
  geom_smooth(method = lm, color= "grey5")+
  theme_classic(base_size = 16)+
  scale_color_manual(values = IBM)+
  xlab("Number of Species") 
p1

p2<-df  %>%
  ggplot(aes(x=n_species, y=Root.Biomass, colour = Treatment )) +
  geom_jitter(width = .1, size=1 )+
  geom_smooth(method = lm, color= "grey5")+
  theme_classic(base_size = 16)+
  scale_color_manual(values = IBM)+
  xlab("Number of Species")
p2



###bar plot for each treatment##


p1<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
  facet_grid( ~n_species, scales = "free", space = "free")+
  xlab("composition")+
  ylab("root biomass (g dry weight)")
p1




p2<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0, size=14),
        legend.position = "none")+
  facet_grid( ~n_species, scales = "free", space = "free")+
  xlab("composition")+
  ylab("shoot biomass (g dry weight)")

p2


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig2_nspecies")
svg(file="biomass.all.trt.svg",width = 5.2, height=2.5)
require(gridExtra)
#windows(2,6)
grid.arrange(p1, p2, ncol=2)
dev.off()


#### import nfix data and process####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Analysis/N Fixation from Emma")
df <- read.csv("merged.plate.ghbiomass.sheet.csv", header=T, stringsAsFactors = F) # fix data


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

# Calculate biomass N in mg N /g soil #
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


# now we have a few data dataframes
head(dfnoNadd.leg)
head(dfNadd.leg)

#put together N+ and N- ##
dfNadd.leg
dfnoNadd.leg
df.leg<-rbind(dfnoNadd.leg, dfNadd.leg)
df.leg$treatment <- factor(df.leg$treatment)

# calcluated the amount fo N fixed per plant
 df.leg<-df.leg%>%
  mutate(n_legume = recode(spp.number,
                                  "1" = 6,
                                  "2" = 3,
                                  "3" =2)) %>%
   mutate(n_fix_per_legume = totalN.g.pot.1/n_legume)

 # make composition var
 df.leg$composition <- df.leg$treatment
 df.leg$composition<-gsub("L", "Legume", df.leg$composition)
 df.leg$composition<-gsub("G", "Grass", df.leg$composition)
 df.leg$composition<-gsub("B", "Brassica", df.leg$composition)
 df.leg$composition<-gsub("LegumeBrassica", "Legume_Brassica", df.leg$composition)
 df.leg$composition<-gsub("LegumeGrass", "Legume_Grass", df.leg$composition)
 df.leg$composition<-gsub("GrassBrassica", "Grass_Brassica", df.leg$composition)
 df.leg$composition 
 df.leg$composition<- factor(df.leg$composition, levels = c( "Legume", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))
 
 
 
 

############## N fix figures ####################

 
# by treatment 
  label <- df.leg$treatment
 label
  label <- gsub("LGB", "C" ,label )
  label<- gsub("LG", "B" ,label )
  label <- gsub("LB", "B" ,label )
  label<- gsub("L", "A" ,label )

  
p1<- ggplot(df.leg, aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
   geom_jitter(size=.5)+
   ylab('Nitrogen from Fixation (%)') +
   xlab("Treatment") +
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = legume_cols)+
   facet_grid( ~spp.number, scales = "free", space = "free")+
    xlab("Treatment")+
  geom_text(y=90, label = label, size=5)
p1


# Analysis of variance 
one.way.Nadd <- aov(perc.Ndfa ~ treatment, data = dfNadd.leg)
summary(one.way.Nadd) # difference between treatments
tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
print(tukey.result.Nadd) # All difference except LG-LB
#plot(one.way.Nadd) #homoscedasticity looks fine

#plot perc by nspecies
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig2_nspecies")
 #svg(file="nfix.percent.treatment.svg",width = 2.3, height=2.3) 
 p1
 #dev.off() 
 
p2<- ggplot(df.leg, aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
   geom_jitter(size=.5)+
   #ylab('Nitrogen from Fixation (%)') +
   xlab("Treatment") +
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = legume_cols)+
   facet_grid( ~spp.number, scales = "free", space = "free")+
   ylab( "N fixed per legume g")
 
p2
 
 
 # Analysis of variance 
 one.way.Nadd <- aov(perc.Ndfa ~ treatment, data = dfNadd.leg)
 summary(one.way.Nadd) # difference between treatments
 tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
 print(tukey.result.Nadd) # All difference except LG-LB
 #plot(one.way.Nadd) #homoscedasticity looks fine
 
 
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



