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
#library(emmeans)
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

# add block info 
# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_experiment_planning.xlsx")
head(block)  
df<-left_join(df, block)
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


# load stats libraries
library(multcomp)
library(emmeans)
library(lme4)
library(lmerTest)
library(car)

####### shoot biomass overall model
# 1. Fit your ANOVA model
# Note: Since 'block' is an additive block factor, don't include it in emmeans specs
m1 <- aov(Stem.Biomass.g ~ Treatment * N + block, data = df)
summary(m1) # difference between treatments and nitrogen interaction

# 2. Get estimated marginal means grouped by Nitrogen level
emm_object <- emmeans(m1, ~ Treatment | N)

# 3. Perform pairwise comparisons within each Nitrogen level
pairwise_comparison <- pairs(emm_object, adjust = "tukey")
summary(pairwise_comparison)

# Requires multcomp / multcompView packages
cld_results <- cld(emm_object, Letters = letters, adjust = "tukey")
print(cld_results)



######root biomass overall model 
m1<-aov(Root.Biomass.g ~ Treatment*N+block, data = df)
summary(m1) # difference between treatments and nitrogen interaction block is sig


# 2. Get estimated marginal means grouped by Nitrogen level
emm_object <- emmeans(m1, ~ Treatment | N)

# 3. Perform pairwise comparisons within each Nitrogen level
pairwise_comparison <- pairs(emm_object, adjust = "tukey")
summary(pairwise_comparison)

# Requires multcomp / multcompView packages
cld_results <- cld(emm_object, Letters = letters, adjust = "tukey")
print(cld_results)





# 
# 
# 
# 
# ###maybe cut this below>
# # nitrogen -
# dfN0<- df %>% filter(N==0)
# m1<-aov(Root.Biomass.g ~ Treatment+block, data = dfN0)
# summary(m1) 
# 
# emm_object <- emmeans(m1, specs = ~ Treatment)
# # Perform all pairwise comparisons with Tukey adjustment
# pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
# 
# summary(pairwise_comparison)
# cld_result <- cld(emm_object, 
#                   adjust = "tukey", 
#                   alpha = 0.05,
#                   Letters = letters) 
# print(cld_result)
# 
# #nitrogen +
# dfN1<- df %>% filter(N==1)
# m1<-aov(Root.Biomass.g ~ Treatment+block, data = dfN1)
# summary(m1) 
# # Perform all pairwise comparisons with Tukey adjustment
# emm_object <- emmeans(m1, specs = ~ Treatment)
# pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
# summary(pairwise_comparison)
# 
# cld_result <- cld(emm_object, 
#                   adjust = "tukey", 
#                   alpha = 0.05,
#                   Letters = letters) 
# print(cld_result)
# 

#########################prediction#############################
# clear workspace and restart R
rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)

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
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB.predict", "GB", "LB.predict",  "LB",   "LG.predict", "LG", 
                             "LGB.predict"   , "LGB" ))

#### plot predictions from monocultures for biomass

mycols <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
    "grey",  
  "#DC267F", # magenta pink GB
    "grey",
  "#FE6100", # bright orange LB
    "grey",
  "#FFB000", # golden yellow LG
    "grey",
  "#865338" # medium mocha brown LGB
)

# 
# label <- df1$Treatment
# label <- gsub("LGB.predict", "*" ,label )
# label <- gsub("LB.predict", "*" ,label )
# label <- gsub("LG.predict", "*" ,label )
# label<- gsub("L", "" ,label )
# label<- gsub("G", "" ,label )
# label <- gsub("B", "" ,label )
# label <- gsub(".predict", "" ,label )
# label

p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
       )+
  
  #geom_text(y=6, label =label , nudge_x = .5, size=8)+
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
# shoots
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
df2
m1<- lm(Shoot.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Shoot.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Shoot.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Shoot.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)

############ root biomass ###
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Root.Biomass~ Treatment*Nitrogen_label, data=df2)
summary(m1)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Root.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Root.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Root.Biomass~ Treatment*Nitrogen_label, data=df2)
anova(m1)







#### supplement: percent root biomass of each species ####

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
    "#785EF0" # light purple B
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


### bulk assigned to crops ####

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

# recode name
df1<-df%>% filter(Treatment!="B" & Treatment!="L" & Treatment!="G") 

df1$label<-gsub("LGB", "" ,df1$Trt_ID)
df1$label
df1$label<-gsub("GB", "" ,df1$label)
df1$label<-gsub("LG", "" ,df1$label)
df1$label<-gsub("LB", "" ,df1$label)




setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement")
svg("percent.biomass.svg", width = 8, height =4.5 )
df1%>%
  ggplot(aes(fill=Species, y=percent, x=label)) + 
  theme_bw(base_size = 12)+
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values=mono_cols)+
  facet_grid(~Treatment, space="free", scales="free")+
  labs(x= "Sample",
       y="proportion dry biomass (g)")

dev.off()
write.csv(df, "percent.biomass.csv")

