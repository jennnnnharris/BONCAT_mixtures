# Figure 3 prediction


# clear workspace and restart R
rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)

#import colors
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")
mycols7<-c( "#4B2D4BFF", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
mycols4 <- c("#4B2D4BFF",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "LB", "LG", "LGB"))

# load paths
biomasspath <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology"
nfixpath <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Analysis/N Fixation from Emma"
fig2path <-  "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig2_nspecies"
fig3path <-  "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig3_predict"


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
dfb <- df
head(dfb)

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


############## Fig 2  N fix ####################
  
 #set wd
setwd(fig2path)
 
 #plot nspecies perc
 setwd(fig2path)
 svg(file="nfix.perc.species.svg",width = 2.4, height=2.4)
 ggplot(df.leg, aes(x=spp.number, y=perc.Ndfa)) + 
   ylab('Nitrogen from Fixation (%)') +
   xlab("Number of species") +
   geom_jitter(width = .2, size=1 )+
   geom_smooth(method = lm, color= "grey10")+
   theme_classic(base_size = 12)
  dev.off()

  # plot nspecies per legume
  setwd(fig2path)
  svg(file="nfix.perc.species.svg",width = 2.4, height=2.4)
  ggplot(df.leg, aes(x=spp.number, y=n_fix_per_legume)) + 
    ylab('grams N fixed per legume') +
    xlab("Number of species") +
    geom_jitter(width = .2, size=1 )+
    geom_smooth(method = lm, color= "grey10")+
    theme_classic(base_size = 12)
  dev.off()

  
  # plot nspecies with colors
 ggplot(df.leg, aes(x=spp.number, y=perc.Ndfa, colour = treatment)) + 
   ylab('Nitrogen from Fixation (%)') +
   xlab("Number of species") +
   geom_jitter(width = .2, size=2 )+
   geom_smooth(method = lm, color= "grey10")+
   theme_classic(base_size = 12)+
   theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
         plot.title = element_text(hjust = 0.5))+
   scale_color_manual(values=mycols4)

 #per legume
 ggplot(df.leg, aes(x=spp.number, y=n_fix_per_legume, colour = treatment)) + 
   ylab('grams N fixed per legume') +
   xlab("Number of species") +
   geom_jitter(width = .2, size=2 )+
   geom_smooth(method = lm, color= "grey10")+
   theme_classic(base_size = 12)+
   theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
         plot.title = element_text(hjust = 0.5))+
   scale_color_manual(values=mycols4)
 
 # linear model
 m1<-lm( perc.Ndfa~spp.number , data= df.leg )
 summary(m1)
 #plot(m1)
 
 
 #plot perc by nspecies
 svg(file="nfix.prec.treatment.svg",width = 2.3, height=2.3)
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
 
 svg(file="nfix.perleg.treatment.svg",width = 2.3, height=2.3)
 ggplot(df.leg, aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
   geom_jitter(size=.5)+
   #ylab('Nitrogen from Fixation (%)') +
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
 #plot(one.way.Nadd) #homoscedasticity looks fine
 
 
 
 
 
### figure N fixed######
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



# Fig 3 write functions
get.predict<-function(df, sp1, sp2, sp3, trait) {
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
    sp1.df<- filter(df, Treatment==sp1) %>% select(all_of(trait))
    sp1.predict<-sp1.df/2
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$trait/2
    sp3.df<-filter(df, Treatment==sp3)
    sp3.predict<-sp3.df$trait/2
    predict <- sp1.predict + sp2.predict + sp3.predict
    print(predict)
    return(predict)
  }
  
}
get.shoot.predict<-function(df, sp1, sp2, sp3) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(Shoot.Biomass)
    sp1.predict<-sp1.df$Shoot.Biomass/3
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$Shoot.Biomass/3
    predict <- sp1.predict + sp2.predict
    print(predict)
    return(predict)
  }
  else
    
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(Shoot.Biomass)
  sp1.predict<-sp1.df$Shoot.Biomass/2
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$Shoot.Biomass/2
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$Shoot.Biomass/2
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  
}
get.root.predict<-function(df, sp1, sp2, sp3) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(Root.Biomass)
    sp1.predict<-sp1.df$Root.Biomass/3
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$Root.Biomass/3
    predict <- sp1.predict + sp2.predict
    print(predict)
    return(predict)
  }
  else
    
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(Root.Biomass)
  sp1.predict<-sp1.df$Root.Biomass/2
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$Root.Biomass/2
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$Root.Biomass/2
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  
}


######### Fig 3 jenny predictions from monocultures for biomass #######


## we expect that a plant makes the same amount of biomass in monoculures vs mixtures
#example
#LB biomass = L monoculture/  2 + B monocultre/2 

#get shoot predictions
dfb<-dfb %>% group_by(Rep, N)
dfb



GB<-get.shoot.predict(dfb, "G", "B") 
LB<-get.shoot.predict(dfb, "L", "B") 
LG<-get.shoot.predict(dfb, "L", "G") 
LGB<-get.shoot.predict(dfb, "L", "G", "B")  
Shoot.Biomass<- round(as.numeric(c(GB, LB, LG, LGB)), 2)
Shoot.Biomass
Treatment<-c(rep("GB.predict", n_groups(dfb)), rep("LB.predict", n_groups(dfb)), rep("LG.predict", n_groups(dfb)), rep("LGB.predict", n_groups(dfb)) )
Treatment
predict <- as.data.frame(cbind(Treatment, Shoot.Biomass))
predict$Shoot.Biomass<-as.numeric(predict$Shoot.Biomass)
predict

#get root predictions

GB<-get.root.predict(dfb, "G", "B") 
LB<-get.root.predict(dfb, "L", "B") 
LG<-get.root.predict(dfb, "L", "G") 
LGB<-get.root.predict(dfb, "L", "G", "B")  
Root.Biomass<- round(as.numeric(c(GB, LB, LG, LGB)), 2)

predict<- cbind(predict, Root.Biomass)
predict$Root.Biomass<-as.numeric(predict$Root.Biomass)
predict


## add to dfb
dfb1<-dfb %>% filter(n_species!="1")
mixtures<-dfb %>% filter(n_species!="1")
47*2

dfb1<-full_join(dfb1, predict) 
dfb1<-dfb1 %>% ungroup()


as.factor(dfb1$Treatment)


######### Fig 3 plot jenny predictions from monocultures for biomass################

mycols<-c("#365C83FF", "grey" ,"#AD5A6BFF", "grey", "#E3C1CBFF", "grey", "#384351FF", "grey")

setwd(fig3path)
svg(file="biomass.root.predict.svg",width = 2.8, height=3)
dfb1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.8 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))
dev.off()




setwd(fig3path)
svg(file="biomass.shoot.predict.svg",width = 4, height=3)
dfb1  %>% 
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.8 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))
dev.off()



######### Fig 3 jenny predictions from monocultures for Nfix######
 
mycols4 <- c("grey", "#AD5A6BFF", "#E3C1CBFF", "#384351FF")

setwd(fig3path)
svg(file="nfix.perplant.predict.svg",width = 2.5, height=2.5)
df.leg %>%
ggplot(aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
  geom_jitter(width = .2, size=.8 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols4) +
  scale_fill_manual(values = mycols4)+
  ylab("grams of nitrogen fixed per legume") +
  theme_classic(base_size = 12)+
  theme(legend.position = "none", 
  axis.text.x = element_text(angle=60, hjust=1),
      plot.title = element_text(hjust = 0.5))
  
dev.off()

  
#####figure percent N ###
setwd(fig3path)
svg(file="nfix.percent.predict.svg",width = 2.5, height=2.5)
df.leg %>%
ggplot( aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(width = .2, size=.8 )+
  ylab('Nitrogen from Fixation (%)') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = mycols4)
dev.off()



##### fig 3 anova#####

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





