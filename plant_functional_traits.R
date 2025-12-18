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
df$long_name <- df$Treatment
df$long_name<-gsub("L", "Legume", df$long_name)
df$long_name<-gsub("G", "Grass", df$long_name)
df$long_name<-gsub("B", "Brassica", df$long_name)
df$long_name<-gsub("LegumeBrassica", "Legume_Brassica", df$long_name)
df$long_name<-gsub("LegumeGrass", "Legume_Grass", df$long_name)
df$long_name<-gsub("GrassBrassica", "Grass_Brassica", df$long_name)
df$long_name 
df$long_name<- factor(df$long_name, levels = c("Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))

# rename
dfb <- df
head(dfb)

###biomass figures ##########


p1<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
  ylab("root biomass (g dry weight)")
p1




p2<-df  %>% filter(n_species!="NA") %>%
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0, size=14),
        legend.position = "none")+
  ylab("shoot biomass (g dry weight)")

p2



grid.arrange(p1, p2, ncol=2)



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
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
#write.csv(df.leg, "Nfix.csv")

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
 
 
 # by treatment 
 label <- df.leg$treatment
 label
 
 label <- gsub("LGB", "B" ,label )
 label<- gsub("LG", "AB" ,label )
 label <- gsub("LB", "AB" ,label )
 label<- gsub("L", "A" ,label )
 
p2<- ggplot(df.leg, aes(x=treatment, y=n_fix_per_legume, fill=treatment)) + 
   geom_boxplot(alpha=.7, outlier.shape = NA)+
   geom_jitter(size=.5)+
   #ylab('Nitrogen from Fixation (%)') +
   xlab("Treatment") +
   theme_classic(base_size = 12) +
   theme(legend.position = "none")+
   scale_fill_manual(values = legume_cols)+
   facet_grid( ~spp.number, scales = "free", space = "free")+
   ylab( "N fixed per legume g") +
 geom_text(y=.07, label = label, size=5)
p2
 

 # Analysis of variance 
 one.way.Nadd <- aov(n_fix_per_legume ~ treatment, data = df.leg)
 summary(one.way.Nadd) # difference between treatments
 tukey.result.Nadd <- TukeyHSD(one.way.Nadd)
 print(tukey.result.Nadd) # All difference except LG-LB
 #plot(one.way.Nadd) #homoscedasticity looks fine
 
 
 #scale_shape_discrete() 
 require(gridExtra)
 grid.arrange(p1, p2, ncol=2)
 
 
  
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



######## weed seed decay #######
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
 IBM2 <- c( #IBM colors
   "grey",
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
df <- read_csv("GH Mix Germination Data(Sheet1).csv")
head(df)

# set factor
df$treatment   <- factor(df$treatment, levels= c("S", "L", "G", "B", "GB", "LB", "LG", "LGB"))
unique(df$treatment)

# look at daeta distrubution 
hist(df$prop_nongerm)

# foxtail
#L G B LB GB LG LGB
# a a a ab ab ab b
df1<-df  %>%   filter(treatment!="S") %>% filter(species=="foxtail")
lab = as.character(df1$treatment)
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


p1<- df %>% filter(species=="foxtail")%>%  filter(treatment!="S") %>%
  ggplot(aes(x=treatment, y=prop_nongerm, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = IBM)+
  ylab("percent non germinating foxtail seeds")+
  geom_text(y=.21, label = lab, size=5)
p1





# pigweed
# L   G   B   GB    LB  LG  LGB
# bc abc  d   ab    a   abc cd
df1<-df  %>%   filter(treatment!="S") %>% filter(species=="pigweed")
lab = as.character(df1$treatment)
unique(lab)
lab<-gsub("LGB", "CD", lab)
lab<-gsub("GB", "AX", lab)
lab<-gsub("LG", "AXC", lab)
lab<-gsub("LB", "A", lab)
lab<-gsub("B", "D", lab)
lab<-gsub("G", "ABC", lab)
lab<-gsub("L", "XC", lab)
lab<-gsub("X", "B", lab)

p2<-df %>% filter(species=="pigweed")%>%  filter(treatment!="S") %>%
  ggplot(aes(x=treatment, y=prop_nongerm, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(legend.position = "none")+
  scale_fill_manual(values = IBM)+
  ylab("percent non germinating pigweed seeds")+
  geom_text(y=.73, label = lab, size=5)
p2


require(gridExtra)
grid.arrange(p1, p2, ncol=2)





# model - binomial model with percent data##
# make vector of successes and failures
y<-cbind(df$num_nongerm, df$num_germ)
m1<-glm(data= df, y~(treatment+block+species)^2, family = binomial)
summary(m1)
# significant effect of block and species

# model without soil
df <- df %>% filter(treatment!="S")
y<-cbind(df$num_nongerm, df$num_germ)
m1<-glm(data= df, y~(treatment+block+species)^2, family = binomial)
summary(m1)
# significant effect of block and species and species treatment interaction


# model foxtail
df1 <- df %>% filter(species=="foxtail")
y<-cbind(df1$num_nongerm, df1$num_germ)
m1<-glm(data= df1, y~treatment+block, family = binomial)
summary(m1)
# significant effect of block and species and species treatment interaction
# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
# This will perform the pairwise tests on the log-odds scale, 
# apply the sidek adjustment, and assign letters based on the results.
library(multcomp)
cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  # The Letters argument is optional, but common for CLDs
                  Letters = letters) 

print(cld_result)



# model pigweed
df2 <- df %>% filter(species=="pigweed")
y<-cbind(df2$num_nongerm, df2$num_germ)
m1<-glm(data= df2, y~treatment+block, family = binomial)
summary(m1)

# significant effect of block and species and species treatment interaction
# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
  emm_object <- emmeans(m1, specs = ~ treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
# This will perform the pairwise tests on the log-odds scale, 
# apply the sidek adjustment, and assign letters based on the results.
library(multcomp)
cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  # The Letters argument is optional, but common for CLDs
                  Letters = letters) 

print(cld_result)

# L   G   B   GB    LB  LG  LGB
# bc abc  d   ab    a   abc cd



#### percent root biomass of each species ####

# load libraries and cols
library(readxl)
library(tidyverse)
IBM <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
  "#DC267F", # magenta pink GB
  "#FE6100", # bright orange LB
  "#FFB000", # golden yellow LG
  "#865338" # medium mocha brown LGB
)
mono_cols <- 
  c( #IBM colors
    "navy", # dark royal blue L
    "#648FFF", # french blue G
    "#785EF0", # light purple B
    )
# import data frame 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read_excel("biomass_species.xlsx")

# tidy data
head(df)

# make plot
ggplot(df, aes(fill=Species, y=Root.Biomass.g, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# make total.root.species.g column
total<-df %>% group_by(Trt_ID) %>%
  summarise(
  Total.no.bulk = sum(Root.Biomass.g))

# add col to df
df<-left_join(df, total)

# make percent col
df<-df%>% mutate(
  percent= (Root.Biomass.g/Total.no.bulk)*100)

# make plot percent
ggplot(df, aes(fill=Species, y=percent, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#add unknown bulk to df
  
bulk<-df %>%
  select(Treatment, N, Rep, Trt_ID,  Bulk.Root.g) %>%
  group_by(Trt_ID, Treatment, N, Rep) %>%
  summarise(
    Root.Biomass.g= sum(Bulk.Root.g)) %>%
  mutate(
    Species="bulk",
    Brassicae= 0,
    Legume=  0,
    Grass= 0
  )


df1 <- df %>% select( Treatment, N, Rep, Trt_ID, Species,  Brassicae, Legume, Grass, Root.Biomass.g )
df1<-full_join(df1,bulk)

#add total+bulk to df 

total<-df %>%
  select(Trt_ID, Total.Root.g) %>%
  group_by(Trt_ID) %>%
  summarise(
    Total.Root.g= sum(Total.Root.g))

df1<-left_join(df1, total)

# calculate percent
df1<-df1 %>% mutate(
  percent = Root.Biomass.g/Total.Root.g
)

df1$Species<- factor(df1$Species, levels= c("bulk", "legume", "grass", "brassica"))
# make plot percent
ggplot(df1, aes(fill=Species, y=percent, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


### block with bulk assigned to crops ####

# import data frame 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read_excel("biomass_species.xlsx")

# tidy data
head(df)

# make plot
ggplot(df, aes(fill=Species, y=Root.Biomass.g, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# make total.root.species.g column
total<-df %>% group_by(Trt_ID) %>%
  summarise(
    Total.withbulk = sum(Total.Root.g))

# add col to df
df<-left_join(df, total)

# make percent col
df<-df%>% mutate(
  percent= (Total.Root.g/Total.withbulk)*100)

df$Species<- factor(df$Species, levels= c("legume", "grass", "brassica"))

# make plot percent
mono_cols <- 
  c( #IBM colors
    "navy", # dark royal blue L
    "#98b2fa", # french blue G
    "#785EF0" # light purple B
  )
ggplot(df, aes(fill=Species, y=percent, x=Trt_ID)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values=mono_cols)

write.csv(df, "percent.biomass.csv")

