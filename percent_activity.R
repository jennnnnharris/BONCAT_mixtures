# Mixtures Boncat
# created: March 2023
# last edited: April 25
# author: Jennifer Harris


rm(list=ls())
library(readxl)
library(tidyverse)
library(lubridate)

#import colors
#mycols8<- c( "grey", "#1F78B4",  "#eb05db", "#33A02C", "#6A3D9A", "#1a635a","#FF7F00")
#mycols7<-c("#E3C1CBFF", "#AD5A6BFF", "#C993A2FF", "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
mycols7<-c( "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF", , "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
mycols8<-c("grey", "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF",  "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
#mycols7 <-c("#4B2D4BFF", "#3C3C5AFF", "#4B6987FF", "#789696FF", "#968787FF", "#D2C3C3FF", "#875A2DFF", "#873C3CFF")

# import data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")

df<-read_excel("Flow_cyto_master.xlsx", sheet = 2)

#make dates be dates
df$Date_sorted<-ymd(df$Date_sorted)

#make n species column
df$Grass<-as.numeric(df$Grass)
df$Legume<-as.numeric(df$Legume)
df$Brassicae<-as.numeric(df$Brassicae)
df$n_species<-rowSums(df[,5:7])

# filter out day were pos ctl didn't work
df<-filter(df, Date_sorted != "2023-05-25")
head(df)

# make treatment and day factors
df$Treatment   <- factor(df$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
df$Date_sorted   <- factor(df$Date_sorted)

################we need normalize by the day/rep ###############
#for soil samples -- they were run on the same day
## there were 4 samples from June 15 that were negative
## 3 soil samples and 1 brassicae. These values we didn't adjust left them as raw values
# think i should probably use a link logit model to this in the future, however back transforming the effect size can be tricky.
#linear model to get coefs to adjusted see binomial model script for binomcial model.
#m1<-lm(data=df, BONCAT_freq~Treatment + Date_sorted + Rep + Block)
#plot(m1)
#coef(m1)

# not adjusted plot:
df%>% 
  ggplot(aes(x=Date_sorted, y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 18, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))

# adjusted plot
df%>% 
  ggplot(aes(x=Date_sorted, y=BONCAT_freq_adj)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 18, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))

# avg the technical reps that are adj for day.

df1<-df%>% group_by(Species1, Species2, Species3, n_species, Nitrogen, Grass, Legume, Brassicae,   Rep, Block, Treatment, Label_short) %>%
  summarise(BONCAT_freq = mean(BONCAT_freq_adj), 
            n_events_cells= round(mean(n_events_cells), digits = 0),
            n_events_BONCAT= round(mean(success_adj), digits = 0) , 
            n = n())

df1$n_events_cells
###### plots ###########

# increasing species boxplot
# adjusted and technical reps averaged:
df1%>%  filter(Species1!="Soil") %>%
  ggplot(aes(x=as.factor(n_species), y=BONCAT_freq)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 18, )+
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
ylab("percent active")

df$Treatment   <- factor(df$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))


# plot for each treatment

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="activity.svg",width = 7, height=4)

  df1  %>%
  ggplot(aes(x=Treatment, y=BONCAT_freq, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0.5))+
  ylab("percent active")+
  facet_grid( ~n_species, scales = "free", space = "free")+
  ggtitle("number of species")+
  geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LB", "*", "")), size=10)+
  geom_text(aes(,y=18, label = ifelse(df1$Treatment=="LG", "*", "")), size=10)

dev.off()  
#



#note about outliers: the high value outlier in LG treatment is the avg of 2 technical reps. (sample #56)
# the high outlier in teh LGB treatment (sample # 88) does not have technical reps.
df1 %>% filter(Treatment=="LG")

## linear model for activity ###########
# in simple lm LB and LG are higher than soil.
m1<-lm(BONCAT_freq  ~Treatment + Block + Rep, data=df1)
summary(m1)

t<-df1%>% filter(Treatment!="Soil")
m1<-lm(BONCAT_freq  ~Treatment + Block, data=t)
summary(m1)
# LG is different than the base line. 

#mixed model
#library(lme4)
#lm<-lme4::lmer(data=df, BONCAT_freq~Treatment + (1|Date_sorted))
#library(lmerTest)
#lm<-lmer(data=df, BONCAT_freq~Treatment + (1|Date_sorted))
#s<-summary(lm)
#anova(lm)



#---------binomial models ---------
#binomail model on # failures # successes and proportion of each ##

#prop<-df1 %>%
#    mutate(n_failures =  n_events_cells-n_events_BONCAT)
#y<-cbind(prop$n_events_BONCAT, prop$n_failures)
# make vector of successes and failures
#m1<-glm(data= prop, y~Treatment +Block, family = binomial)
#m1
#summary(m1)
#plot(m1)
# I need a quasibinomial model because residual deviance is higher than the degrees of freedom = model is over disposed

## quasibinomial glm on success and failures##
#m1<-glm(data= prop, y~Treatment +Block , family = quasibinomial)
#m1
#summary(m1)
#anova(m1, test= "LRT")  
#anova(m1, test= "Chisq")

# post hoc tests
# subset to compare treatments
#prop$Treatment <- as.character(prop$Treatment)

# no soil
#prop1<-prop%>% filter(Treatment!="Soil")
#y<-cbind(prop1$n_events_BONCAT, prop1$n_failures)
#prop1$Treatment
#m1<-glm(y~prop1$Treatment,  quasibinomial)
#anova(m2, test= "LRT")
#summary(m1)

# L verse G
#prop1<-prop%>% filter(Treatment=="L"|Treatment=="G")
#y<-cbind(prop1$n_events_BONCAT, prop1$n_failures)
#prop1$Treatment
#m2<-glm(y~prop1$Treatment,  quasibinomial)
#anova(m2, test= "LRT")

######binomial model with percent data##############
# make vector of successes and failures
prop<-df1 %>%
  mutate(BONCAT_freq = round(BONCAT_freq, 0)) %>%
  mutate(n_failures =  100-BONCAT_freq)
y<-cbind(prop$BONCAT_freq, prop$n_failures)
y
#model
m1<-glm(data= prop, y~Treatment +Block, family = binomial)
m1
summary(m1)
#plot(m1)
#plots look okay

#block not significant so we can take block out
m1<-glm(data= prop, y~Treatment, family = binomial)
m1
summary(m1)
# all treatments are different that than the soil
#plot(m1)

# post hoc tests
# subset to compare treatments
#prop$Treatment <- as.character(prop$Treatment)

# no soil
prop1<-prop%>% filter(Treatment!="Soil")
prop1$Treatment   <- factor(prop1$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

prop1$Treatment
y<-cbind(prop1$BONCAT_freq, prop1$n_failures)
prop1$Treatment
y
m1<-glm(data=prop1, y~Treatment,  binomial)
#anova(m2, test= "LRT")
summary(m1)

hist(df1$BONCAT_freq, breaks=10)
plot(df1$BONCAT_freq ~ df1$n_species, las=1)


# post hoc tests
# subset to compare treatments
prop$Treatment <- as.character(prop$Treatment)
# no soil
prop1<-prop%>% filter(Treatment!="Soil")
y<-cbind(prop1$n_events_BONCAT, prop1$n_failures)
prop1$Treatment
m2<-glm(y~prop1$Treatment,  binomial)
summary(m2)
#anova(m2, test= "LRT")
#there are difference among treatments when soil is removed.


#########################effect size + predictions ##################
# grab effect size
effectsize<-s$coefficients
effectsize<-data.frame(effectsize)
colnames(effectsize) <- c("effect.size", "std.error", "df", "tvalue", "pvalue") 
head(effectsize)


#make summary table
avg <-df %>% group_by(Treatment) %>%
summarise(mean.activity = mean(BONCAT_freq), sd.activity = sd(BONCAT_freq))
table<-cbind(effectsize, avg)
head(table)

#expectations
# G + B
gb<- table[which(table$Treatment == "Soil"), 1]+ table[which(table$Treatment == "B"), 1] + table[which(table$Treatment == "G"), 1] # mean
gbe<-sqrt(table[which(table$Treatment == "Soil"), 2]+ table[which(table$Treatment == "B"), 2] + table[which(table$Treatment == "G"), 2])# standard dev 
# 15.0 + or  - 2.67 
#actually
#table[which(table$Treatment == "GB"), 1] 
# 6.7 + or - 5.3


# L + B
lb<- table[which(table$Treatment == "Soil"), 1]+ table[which(table$Treatment == "B"), 1] + table[which(table$Treatment == "L"), 1] #mean
lbe<-sqrt(table[which(table$Treatment == "Soil"), 2]+table[which(table$Treatment == "B"), 2] + table[which(table$Treatment == "G"), 2]) #stdev
# expectation:  12.2 +/- # 3.06
#table[which(table$Treatment == "LB"), 2] #mean
#table[which(table$Treatment == "LB"), 3]  #stdev
# actual : 10.5 +/- 7


# L + G
lg<-table[which(table$Treatment == "Soil"), 1]+  table[which(table$Treatment == "G"), 1] + table[which(table$Treatment == "L"), 1] #mean
lge<-sqrt(table[which(table$Treatment == "Soil"), 2]+table[which(table$Treatment == "L"), 2] + table[which(table$Treatment == "G"), 2]) #stdev
# expectation:  8.6 +/- # 2.5
#table[which(table$Treatment == "LG"), 2] #mean
#table[which(table$Treatment == "LG"), 3]  #stdev
# actual : 12.8 +/- 8.5

# L + G + B
lgb<- table[which(table$Treatment == "Soil"), 1]+ table[which(table$Treatment == "G"), 1] + table[which(table$Treatment == "L"), 1] + table[which(table$Treatment == "B"), 1] #mean
lgbe<-sqrt(table[which(table$Treatment == "Soil"), 2]+table[which(table$Treatment == "L"), 2] + table[which(table$Treatment == "G"), 2] +  table[which(table$Treatment == "B"), 2] )#stdev
# expectation:  16.6 +/- 3.5
#table[which(table$Treatment == "LGB"), 2] #mean
#table[which(table$Treatment == "LGB"), 3]  #stdev
# actual : 7.8 +/- 7.8


###predictors verse actual
head(table)
prediction<-c(NA,NA, NA,NA, gb, lb, lg, lgb)
prediction.error <- c(NA, NA, NA,NA, gbe, lbe, lge, lgbe)
table<-cbind(table, prediction)
table <- cbind(table, prediction.error)
head(table)
table$color <- c(0, 0,0, 0, 1, 1, 1 ,1)

table %>% ggplot(aes(x=Treatment, y=prediction, col="red")) +
            geom_boxplot()+
            geom_point(aes(x=Treatment, y=mean.activity, col="black")) 

####plot

  #ggplot(data=df, aes(Treatment, BONCAT_freq)) +
  #geom_jitter(width = .2, size=1 )+
  #geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  #theme_bw(base_size = 22, )+
  #theme(axis.text.x = element_text(angle=60, hjust=1)) +
  #geom_point(data= table, aes(Treatment, y= prediction, col="black")) +
  #geom_point(data= table, aes(Treatment, y= mean.activity, col="red")) 


mycols <- ("black", "red")

mycols<-c(rep("black",4), rep("red", 4))
  plot(
    x = table$Treatment,
    y = table$mean.activity,
    xlab = "Treatment",
    ylab = "Microbial Activity",
    pch = 20, # solid dots increase the readability of this data plot
    col = mycols,
    fill = mycols
  )
  


legend(
  x ="topleft",
  legend = paste("Color", levels(diamonds$color)), # for readability of legend
  col = diamond_color_colors,
  pch = 19, # same as pch=20, just smaller
  cex = .7 # scale the legend to look attractively sized
)


##############################################
# increasing species line w/o brasssisae 
df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Brassicae))) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  stat_summary(geom = "line", fun = mean)

# increasing species line w/o grass
df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Grass))) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  stat_summary(geom = "line", fun = mean)


# increasing species line w/o legume
df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Legume))) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  stat_summary(geom = "line", fun = mean)






df%>% filter(Treatment!="ctl") %>%
  ggplot(aes(x=n_species, y=Percent_BONCAT)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("N_species")

m1<-lm(Percent_BONCAT~n_species, data = df)
summary(m1)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/")
svg(file="brassica.svg", width=5, height=5 )
#windows(6,4)
df%>% filter(Treatment!="ctl", method=="filter", n_species=="1", Percent_BONCAT>0.1) %>%
    group_by(Brassica)%>%
    summarise(mean=mean(Percent_BONCAT), sd=sd(Percent_BONCAT)) %>%
ggplot(aes(x=Brassica, y=mean, fill=Brassica)) +
  geom_bar( position = "dodge", stat = "identity", alpha=.7)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9), col="grey28") +

  scale_fill_manual(values=c("#F4743B", "#70AE6E"))+
  theme_bw(base_size = 22 )+
  theme(legend.position="none")+
  xlab("Brassica")+
  ylab("% active microbes")

dev.off()

df1<-df%>% filter(Treatment!="ctl", method=="filter", n_species=="1", Percent_BONCAT>0.1)
m1<-lm(Percent_BONCAT~Brassica, data=df1)
summary(m1)


df%>% filter(Treatment!="ctl") %>%
  ggplot(aes(x=Grass, y=Percent_BONCAT)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("Grass")



df%>% filter(Treatment!="ctl") %>%
  ggplot(aes(x=Legume, y=Percent_BONCAT)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("")




#####################colored line plots#########################
# increasing species line w/o brasssisae 
df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Brassicae))) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  stat_summary(geom = "line", fun = mean)

# increasing species line w/o grass
df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Grass))) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  stat_summary(geom = "line", fun = mean)


# increasing species line w/o legume
df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Legume))) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  stat_summary(geom = "line", fun = mean)

# Binvary plots #######################################Binvary plots ###########################################
# not a clear pattern on any.
# Brassicae
df1%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=as.factor(Brassicae), y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))
  #scale_y_log10()+
  #xlab("method")

df1%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=as.factor(Grass), y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))
#scale_y_log10()+
#xlab("method")
 
df1%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=as.factor(Legume), y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))
#scale_y_log10()+
#xlab("method")

