# Figure 2 N fix and biomass 


# clear workspace and restart R
rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)

#old cols:
#mycols7<-c( "#715b8a", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
#mycols4 <- c("#715b8a",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )

#current cols"
mycols7<-c("#466F9DFF", "#91B3D7FF",  "#ED444AFF", "#FEB5A2FF", "#9D7660FF", "#D7B5A6FF", "#3896C4FF" )
mycols4<-c("#466F9DFF", "#9D7660FF", "#D7B5A6FF", "#3896C4FF" )
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "LB", "LG", "LGB"))

# load paths
biomasspath <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology"
nfixpath <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Analysis/N Fixation from Emma"
fig2path <-  "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig2_nspecies"
#fig3path <-  "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig3_predict"


#### import biomass data and process #####
setwd(biomasspath)
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

setwd(fig2path)

#no colors
p1<-df  %>%
  ggplot(aes(x=n_species, y=Shoot.Biomass )) +
  geom_jitter(width = .1, size=1 )+
  #geom_jitter(size=1)+
  geom_smooth(method = lm, color=  "grey10" )+
  theme_classic(base_size = 12)+
  theme(lot.title = element_text(hjust = 0.5))+
  xlab("Number of Species") 
#scale_color_manual(values = mycols7)
#geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
p1

p2<-df  %>%
  ggplot(aes(x=n_species, y=Root.Biomass )) +
  geom_jitter(width = .1, size=1 )+
  #geom_jitter(size=1)+
  geom_smooth(method = lm, color=  "grey10")+
  theme_classic(base_size = 12)+
  theme(lot.title = element_text(hjust = 0.5))+
  xlab("Number of Species") 
#scale_color_manual(values = mycols7)
#geom_text(aes(,y=4, label = ifelse(df$Treatment=="LB", "", "")), size=10)
p2


setwd(fig2path)
svg(file="biomass.species.svg",width = 5, height=2.5)
require(gridExtra)
#windows(2,6)
grid.arrange(p1, p2, ncol=2)
dev.off()


##with colors#######################3
p1<-df  %>%
  ggplot(aes(x=n_species, y=Shoot.Biomass, colour = Treatment )) +
  geom_jitter(width = .1, size=1 )+
  geom_smooth(method = lm, color= "grey5")+
  theme_classic(base_size = 16)+
  scale_color_manual(values = mycols7)+
  xlab("Number of Species") 
p1

p2<-df  %>%
  ggplot(aes(x=n_species, y=Root.Biomass, colour = Treatment )) +
  geom_jitter(width = .1, size=1 )+
  geom_smooth(method = lm, color= "grey5")+
  theme_classic(base_size = 16)+
  scale_color_manual(values = mycols7)+
  xlab("Number of Species")
p2


setwd(fig2path)
svg(file="biomass.species.col.svg",width = 10, height=4)
require(gridExtra)
#windows(2,6)
grid.arrange(p1, p2, ncol=2)
dev.off()
###bar plot for each treatment##


p1<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
  facet_grid( ~n_species, scales = "free", space = "free")+
  xlab("Treatment") 
p1
setwd(fig2path)
svg(file="biomass.root.trt.svg",width = 2.4, height=2.5)
p1
dev.off()



setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_liana_kayla")
p2<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0, size=14))+
  facet_grid( ~n_species, scales = "free", space = "free")+
  xlab("Treatment")

p2



#setwd(fig2path)
svg(file="biomass.shoot.trt.svg",width = 3.3, height=2.5)
p2

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

#number of species
#lm
m1<-lm(Root.Biomass~n_species ,data=df)
summary(m1)
plot(m1)
# number of species
#Adjusted R-squared:  0.002231 
#F-statistic: 1.092 on 2 and 80 DF,  p-value: 0.3406


###roots

#n species
m1<-lm(Shoot.Biomass~n_species*N*Block ,data=df)
summary(m1)

#Legume
m1<-lm(Shoot.Biomass~ Legume*Block*N ,data=df)
summary(m1)
plot(m1)



#### import nfix data and process####
setwd(nfixpath)
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
  
 #set wd
setwd(fig2path)
 
 #plot nspecies perc
 #setwd(fig2path)
 #svg(file="nfix.perc.species.svg",width = 2.4, height=2.4)
 #ggplot(df.leg, aes(x=spp.number, y=perc.Ndfa)) + 
#   ylab('Nitrogen from Fixation (%)') +
#   xlab("Number of species") +
#   geom_jitter(width = .2, size=1 )+
#   geom_smooth(method = lm, color= "grey10")+
#   theme_classic(base_size = 12)
  #dev.off()

  # plot nspecies per legume
  #setwd(fig2path)
  #svg(file="nfix.perc.species.svg",width = 2.4, height=2.4)
 # ggplot(df.leg, aes(x=spp.number, y=n_fix_per_legume)) + 
#    ylab('grams N fixed per legume') +
#    xlab("Number of species") +
#    geom_jitter(width = .2, size=1 )+
#    geom_smooth(method = lm, color= "grey10")+
#    theme_classic(base_size = 12)
  #dev.off()

 #per legume
# ggplot(df.leg, aes(x=spp.number, y=n_fix_per_legume, colour = treatment)) + 
#   ylab('grams N fixed per legume') +
#   xlab("Number of species") +
#   geom_jitter(width = .2, size=2 )+
#   geom_smooth(method = lm, color= "grey10")+
#   theme_classic(base_size = 12)+
#   theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
#         plot.title = element_text(hjust = 0.5))+
#   scale_color_manual(values=mycols4)
 
 # linear model
 m1<-lm( perc.Ndfa~spp.number , data= df.leg )
 summary(m1)
 #plot(m1)
 
 
 
p1<- ggplot(df.leg, aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
   geom_jitter(size=.5)+
   ylab('Nitrogen from Fixation (%)') +
   xlab("Treatment") +
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = mycols4)+
   facet_grid( ~spp.number, scales = "free", space = "free")
p1
#plot perc by nspecies
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_liana_kayla")
 svg(file="nfix.percent.treatment.svg",width = 2.3, height=2.3) 
 p1
 dev.off() 
 
 #svg(file="nfix.perleg.treatment.svg",width = 2.3, height=2.3)
 ggplot(df.leg, aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
   geom_jitter(size=.5)+
   #ylab('Nitrogen from Fixation (%)') +
   xlab("Treatment") +
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = mycols4)+
   facet_grid( ~spp.number, scales = "free", space = "free")
 
 #dev.off() 
 
 
 
 # Analysis of variance 
 one.way.Nadd <- aov(perc.Ndfa ~ treatment, data = dfNadd.leg)
 summary(one.way.Nadd) # difference between treatments
 tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
 print(tukey.result.Nadd) # All difference except LG-LB
 #plot(one.way.Nadd) #homoscedasticity looks fine
 
 
 #### figure 2 all #######
 
 
 
 p1<- ggplot(df.leg, aes(x=spp.number, y=perc.Ndfa, colour = treatment)) + 
   ylab('Nitrogen from Fixation (%)') +
   xlab("Number of species") +
   geom_jitter(width = .1, size=1.5 )+
   geom_smooth(method = lm, color= "grey10")+
   theme_classic(base_size = 12)+
   scale_color_manual(values=mycols4)+
   annotate("text", x =2.4, y = 65, label = "p<0.001
            Rsquared=.43 ", color = "black", size = 3.4)
 
 p1
 
 
 
 p2<-dfb  %>%
   ggplot(aes(x=n_species, y=Shoot.Biomass, colour = Treatment )) +
   geom_jitter(width = .1, size=1.5 )+
   geom_smooth(method = lm, color= "grey5")+
   theme_classic(base_size = 12)+
   scale_color_manual(values = mycols7)+
   xlab("Number of Species")+
   annotate("text", x =2.4, y = 4, label = "p=0.461", color = "black", size = 3.4)
 
p2
 
 m1<-lm( Shoot.Biomass~n_species , data= dfb )
 summary(m1)
 
 
 p3<-dfb  %>%
   ggplot(aes(x=n_species, y=Root.Biomass, colour = Treatment )) +
   geom_jitter(width = .1, size=1.5 )+
   geom_smooth(method = lm, color= "grey5")+
   theme_classic(base_size = 12)+
   scale_color_manual(values = mycols7)+
   xlab("Number of Species")+
   annotate("text", x =2.4, y = 2, label = "p<0.001
            Rsquared=.14 ", color = "black", size = 3.4)
 p3 
 
 m1<-lm( Root.Biomass~n_species , data= dfb )
 summary(m1)
 
 

 
 setwd(fig2path)
 svg(file="fig2_linegraph.svg",width = 12, height=3)
 require(gridExtra)
 #windows(10,3)
 grid.arrange(p1, p2, p3, ncol=3)
 dev.off()
 
### some other N fixed figures ??######
ggplot(df.leg, aes(x=treatment, y=totalN.mg.g.1, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  ylab('Total Nitrogen Fixed per gram of soil') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = mycols4)+
  facet_grid( ~spp.number, scales = "free", space = "free")

#### figure N fixed ###
ggplot(df.leg, aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  #ylab('Total Nitrogen Fixed per gram of soil') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = mycols4)+
  facet_grid( ~spp.number, scales = "free", space = "free")


#####figure percent N ###
ggplot(df.leg, aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  ylab('Nitrogen from Fixation (%)') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = mycols4)+
  facet_grid( ~spp.number, scales = "free", space = "free")


