# Figure 3 prediction


# clear workspace and restart R
rm(list=ls())
#rstudioapi::restartSession(clean = TRUE)

#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)

#import colors
#install.packages("viridis")
#library(viridis)

mycols7<-c( "#715b8a", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
mycols7vivid<-c( "#715b8a", "#4D8F8BFF", "#b1de64", "#365C83FF", "#bd0262", "#c77597",  "#384351FF")
mycols4 <- c("#715b8a",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )#

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

# load paths
biomasspath <- "C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology"
#nfixpath <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Analysis/N Fixation from Emma"


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




############### Fig 3 write functions###############
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
    sp1.predict<-sp1.df/3
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$trait/3
    sp3.df<-filter(df, Treatment==sp3)
    sp3.predict<-sp3.df$trait/3
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
#old
# new
#mycols<-c("#FEB5A2FF","grey" , "#9D7660FF", "grey", "#D7B5A6FF", "grey", "#3896C4FF" , "grey")
mycols <- c( #IBM colors
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)

#setwd(fig3path)
#svg(file="biomass.root.predict.svg",width = 2.8, height=3)
label <- dfb1$Treatment
label <- gsub("LGB.predict", "*" ,label )
label <- gsub("LB.predict", "*" ,label )
label <- gsub("LG.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label

p2<-dfb1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.8 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  
  geom_text(y=6, label =label , nudge_x = -.8, size=7)


p2

label <- dfb1$Treatment
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "*" ,label )


#setwd(fig3path)
#svg(file="biomass.shoot.predict.svg",width = 4, height=3)
p3<-dfb1  %>% 
  ggplot(aes(x=Treatment, y=Shoot.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.8 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  #geom_text(y=8, label =label , nudge_x = -.8, size=7)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",
        plot.title = element_text(hjust = 0.5))

p3

# put the tow plots together

require(gridExtra)
grid.arrange(p2, p3, ncol=2)



### example calculation  #####


mycols7<-c("#466F9DFF", "#91B3D7FF",  "white", "white", "white", "#D7B5A6FF", "#3896C4FF" )

p1<-df  %>% filter(n_species!="NA") %>%
  filter(Treatment!= "LGB") %>%
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
    xlab("") 
p1


# caculated half for L
sp1.df<- filter(df, Treatment=="L") %>% select(Root.Biomass)
sp1.predict<-sp1.df$Root.Biomass/2
Root.Biomass<-sp1.predict
Treatment <- rep("half_L", length(Root.Biomass))
halfl<-data.frame(Treatment, Root.Biomass)
df1<-full_join(df, halfl)
df1<-df1 %>% filter(Treatment!="LGB" & Treatment!="LB" & Treatment!="B")
df1$Treatment<-factor(df1$Treatment, levels=c("L", "G", "half_L", "GB", "LG"))


# caculated half for G
sp1.df<- filter(df, Treatment=="G") %>% select(Root.Biomass)
sp1.predict<-sp1.df$Root.Biomass/2
Root.Biomass<-sp1.predict
Treatment <- rep("half_G", length(Root.Biomass))
halfg<-data.frame(Treatment, Root.Biomass)
df1<-full_join(df1, halfg)
df1<-df1 %>% filter(Treatment!="GB")
df1$Treatment<-factor(df1$Treatment, levels=c("L", "G", "half_L", "half_G", "LG"))
as.factor(df1$Treatment)



# caculated LG 
sp1<- filter(df1, Treatment=="half_L") %>% select(Root.Biomass)
sp2<- filter(df1, Treatment=="half_G") %>% select(Root.Biomass)

LG.predict <- sp1$Root.Biomass + sp2$Root.Biomass
Root.Biomass<-LG.predict
Treatment <- rep("LG.predict", length(Root.Biomass))
predict<-data.frame(Treatment, Root.Biomass)
df1<-full_join(df1, predict)
df1$Treatment<-factor(df1$Treatment, levels=c("L", "G", "half_L", "half_G", "LG.predict", "LG"))
as.factor(df1$Treatment)

#####example plot#######
mycols7<-c( "#715b8a", "#4D8F8BFF", "grey", "grey", "grey", "#AD5A6BFF") 

mycols <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "grey", # light purple B
  "grey", # magenta pink GB
  "grey", # bright orange LB
  "#FFB000" # golden yellow LG
)

p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=Root.Biomass, fill = Treatment)) +
  geom_jitter(width = .2, size=.5 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0, size=14))+
  xlab("") 
p1

#####figure percent N ###
#svg(file="nfix.percent.predict.svg",width = 2.5, height=2.5)
# rename
df1<-df.leg %>% filter(treatment=="L")
df1$treatment<-gsub("L", "L.predict", df1$treatment )
m<-mean(df1$perc.Ndfa)

df2<-df.leg %>% filter(treatment!="L")
df.leg1<-rbind(df1, df2)
# set cols
mycols4 <- c("grey", "#35b779", "#90d743", "#fde725")

lab <- df.leg1$treatment
lab<- gsub("L.predict", "" , lab )
lab<- gsub("LGB", "*" , lab  )
lab <- gsub("LB", "*" , lab )
lab <- gsub("LG", "*" , lab )
lab


p1<-df.leg1 %>%
  ggplot( aes(x=treatment, y=perc.Ndfa, fill=treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(width = .2, size=.8 )+
  ylab('Nitrogen from Fixation (%)') +
  xlab("Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = mycols4)+
  geom_hline(yintercept = m, color = "grey", linewidth = 1.5, linetype = "dashed")+
  annotate("text", x =3.5, y = 70, label = "   prediction from
           monoculture", color = "black", size = 3.4)+
  geom_text(y=90, label =lab , nudge_x = .1 , size=7)



p1

setwd(fig3path)
svg(file="fig3_predict.svg",width = 10, height=3)
require(gridExtra)
#windows(10,3)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



##### anova#####


#GB
#filter
df2<-dfb1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LB
df2<-dfb1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LG
df2<-dfb1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

#LGB
df2<-dfb1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Root.Biomass~ Treatment, data=df2)
anova(m1)

# shoots
#GB
#filter
df2<-dfb1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LB
df2<-dfb1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LG
df2<-dfb1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)

#LGB
df2<-dfb1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<- lm(Shoot.Biomass~ Treatment, data=df2)
anova(m1)





