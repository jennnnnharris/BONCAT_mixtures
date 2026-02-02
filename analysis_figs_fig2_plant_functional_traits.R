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



#### import biomass data #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read.csv("biomass_potlevel.csv") # biomass data

# process to make n species column, summarize at pot level, make treatment a factor.
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<- df %>% select(c(Brassicae, Legume, Grass )) %>% rowSums()
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))


###biomass figures ##########
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
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
svg("biomass.svg", height = 2.5, width = 5)
grid.arrange(p1, p2, ncol=2)
dev.off()


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
df<-df.leg %>% filter(Nitrogen==0)
one.way.Nadd <- aov(n_fix_per_legume ~ Treatment, data = df)
summary(one.way.Nadd) # difference between treatments
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
  dplyr::select(Treatment, N, Rep, Trt_ID,  Bulk.Root.g) %>%
  group_by(Trt_ID, Treatment, N, Rep) %>%
  summarise(
    Root.Biomass.g= sum(Bulk.Root.g)) %>%
  mutate(
    Species="bulk",
    Brassicae= 0,
    Legume=  0,
    Grass= 0
  )


df1 <- df %>% dplyr:: select( Treatment, N, Rep, Trt_ID, Species,  Brassicae, Legume, Grass, Root.Biomass.g )
df1<-full_join(df1,bulk)

#add total+bulk to df 

total<-df %>%
  dplyr::select(Trt_ID, Total.Root.g) %>%
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
    "#0a3170", # dark royal blue L
    "#cdddf7", # french blue G
    "#785EF0" # light purple B
  )

windows(6,4)
df%>% filter(Treatment!="B" & Treatment!="L" & Treatment!="G") %>%
ggplot(aes(fill=Species, y=percent, x=Trt_ID)) + 
  theme_bw(base_size = 12)+
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values=mono_cols)+
  labs(x= "Treatment ID")

write.csv(df, "percent.biomass.csv")

