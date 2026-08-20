# Jennifer Harris
# biomass 
# Aug 10 2026

# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

#load libraries
library(tidyverse)
library(lubridate)

mycols <- c( #IBM colors
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

#biomass without prediction 

#### import biomass data and process #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
df <- read.csv("biomass_species.csv") # biomass data
df

df<-df %>% group_by(Trt_ID, Treatment, N, Rep, Brassicae, Legume, Grass ) %>%
  summarise(
    Stem.Biomass.g = sum(Stem.Biomass.g),
    Root.Biomass.g = sum(Root.Biomass.g),
    Bulk.Root.g = sum(Bulk.Root.g),
    Total.Root.g = sum(Total.Root.g)
  )

df  


# process to make n species column, summarize at pot level, make treatment a factor.

df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df<-df%>% ungroup() 
#df$n_species<- df %>% select(c(Brassicae, Legume, Grass )) %>% rowSums()
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
df

# added nitrogen column 
df <- df %>%
  mutate(Nitrogen_label = if_else(N == 1, "Nitrogen +", "Nitrogen -"))

# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
block <- read.csv("metadata_blockinfo.csv")
head(block)  
df<-left_join(df, block)
head(df)

p1<-df  %>% filter(Treatment=="soil") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass.g, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=.7, width=.1)+
  
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=12))+
  labs(title = "A",
       x="",
       y="Root biomass (g)")+
  scale_shape_manual(values = c(17, 16)) 
p1

p2<-df  %>% filter(Treatment=="soil") %>%
  ggplot(aes(x=Treatment, y=Stem.Biomass.g, fill = Treatment)) +
  geom_jitter(aes(shape=Nitrogen_label), size=.7, width=.1)+
  
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0, size=12),
        legend.position = "none")+
  labs(title = "B",
       x="",
       y="Shoot biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #
p2
require(gridExtra)
grid.arrange(p1, p2, ncol=2)


####### shoot biomass overall model


# 1. Fit your ANOVA model
# Note: Since 'block' is an additive block factor, don't include it in emmeans specs
m1 <- aov(Stem.Biomass.g ~ Treatment * N + block, data = df)
summary(m1) # difference between treatments and nitrogen interaction

# 2. check for normality and homogeneity of variance
ks.test(df$Stem.Biomass.g, "pnorm", mean = mean(df$Stem.Biomass.g), sd = sd(df$Stem.Biomass.g))
qqnorm(residuals(m1))
qqline(residuals(m1), col = "red") # Adds reference line

# 2. Get estimated marginal means grouped by Nitrogen level
emm_object <- emmeans::emmeans(m1, ~ Treatment | N)

# 3. Perform pairwise comparisons within each Nitrogen level
pairwise_comparison <- pairs(emm_object, adjust = "tukey")
summary(pairwise_comparison)

# Requires multcomp / multcompView packages
cld_results <- multcomp::cld(emm_object, Letters = letters, adjust = "tukey")
print(cld_results)



######root biomass overall model 
m1<-aov(Root.Biomass.g ~ Treatment*N+block, data = df)
summary(m1) # difference between treatments and nitrogen interaction block is sig
# 2. check for normality and homogeneity of variance
ks.test(df$Root.Biomass.g, "pnorm", mean = mean(df$Root.Biomass.g), sd = sd(df$Root.Biomass.g))
qqnorm(residuals(m1))
qqline(residuals(m1), col = "red") # Adds reference line


# 2. Get estimated marginal means grouped by Nitrogen level
emm_object <- emmeans::emmeans(m1, ~ Treatment | N)

# 3. Perform pairwise comparisons within each Nitrogen level
pairwise_comparison <- pairs(emm_object, adjust = "tukey")
summary(pairwise_comparison)

# Requires multcomp / multcompView packages
cld_results <- multcomp::cld(emm_object, Letters = letters, adjust = "tukey")
print(cld_results)



#########################prediction#############################
# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

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

#### import biomass data and process
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
df <- read.csv("biomass_species.csv") # biomass data
df

df<-df %>% group_by(Trt_ID, Treatment, N, Rep, Brassicae, Legume, Grass ) %>%
  summarise(
    Stem.Biomass.g = sum(Stem.Biomass.g),
    Root.Biomass.g = sum(Root.Biomass.g),
    Bulk.Root.g = sum(Bulk.Root.g),
    Total.Root.g = sum(Total.Root.g)
  )

df  


# process to make n species column, summarize at pot level, make treatment a factor.

df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df<-df%>% ungroup() 
df$n_species<- df %>% select(c(Brassicae, Legume, Grass )) %>% rowSums()
df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
df

# added nitrogen column 
df <- df %>%
  mutate(Nitrogen_label = if_else(N == 1, "Nitrogen +", "Nitrogen -"))

# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
block <- read.csv("metadata_blockinfo.csv")
head(block)   
df<-left_join(df, block)
head(df)




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
Treatment<-c(rep("GB_expected", n_groups(df)), rep("LB_expected", n_groups(df)),
             rep("LG_expected", n_groups(df)), rep("LGB_expected", n_groups(df)) )
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
# add rep
predict$Rep<-  rep(1:6,8)



## add to df

df1<-full_join(df, predict) 
df1<-df1 %>% ungroup()
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB_expected", "GB", "LB_expected",  "LB",   "LG_expected", "LG", 
                             "LGB_expected"   , "LGB" ))

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
       y="Shoot biomass (g)")+
  scale_shape_manual(values = c(17, 16)) #

p2

# put the plots together

require(gridExtra)
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/")
#svg("biomasspredict.svg", height = 7, width = 4.5)
grid.arrange(p1, p2, ncol=1)
#dev.off()


######### anova#####
# shoots
df1$Treatment
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB_expected")
df2
m1<- lm(Shoot.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)



#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB_expected")
m1<- lm(Shoot.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG_expected")
m1<- lm(Shoot.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB_expected")
m1<- lm(Shoot.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)

############ root biomass ###
#GB
#filter
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB_expected")
m1<- lm(Root.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
summary(m1)
anova(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB_expected")
m1<- lm(Root.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG_expected")
m1<- lm(Root.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB_expected")
m1<- lm(Root.Biomass~ Rep+Treatment*Nitrogen_label, data=df2)
anova(m1)







#### supplement: percent root biomass of each species ####

# load libraries and cols
library(readxl)
library(tidyverse)

# import data frame 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
df <- read.csv("biomass_species.csv") # biomass data
df

# check df
head(df)

# make total.root.species.g column
total<-df %>% group_by(Trt_ID) %>%
  summarise(
    Total.withoutbulk = sum(Root.Biomass.g))

# add col to df
df<-left_join(df, total)

# make percent col
df<-df%>% mutate(
  percent= (Root.Biomass.g/Total.withoutbulk)*100)
df$Species<- factor(df$Species, levels= c("legume", "grass", "brassica"))

# save data frame
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement")
write.csv(df, "biomass_percent.csv")

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
df$Species<- factor(df$Species, levels= c("legume", "grass", "brassica"))

# make plot
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement")
#svg("percent.biomass.svg", width = 8, height =4.5 )
df1%>%
  ggplot(aes(fill=Species, y=percent, x=label)) + 
  theme_bw(base_size = 12)+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=mono_cols)+
  facet_grid(~Treatment, space="free", scales="free")+
  labs(x= "Sample",
       y="proportion dry biomass (g)")

#dev.off()



