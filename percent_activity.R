# Mixtures Boncat
# March 2023
# last edited: Jan 2025
#Jennifer Harris


rm(list=ls())
library(readxl)
library(tidyverse)
library(lubridate)



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


# increasing species boxplot
df%>% # filter(Species1!="Soil") %>%
  ggplot(aes(x=as.factor(n_species), y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  #geom_line()+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 18, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))

# increasing species line
#df%>% filter(Species1!="Soil") %>%
#  ggplot(aes(x=n_species, y=BONCAT_freq)) +
#  geom_jitter(width = .2, size=2 )+
#  theme_bw(base_size = 18, )+
#  theme(axis.text.x = element_text(angle=60, hjust=1))+
#  stat_summary(geom = "line", fun = mean)

#df%>% filter(Species1!="Soil") %>%
#  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Date))) +
#  geom_jitter(width = .2, size=2 )+
#  theme_bw(base_size = 18, )+
#  stat_summary(geom = "line", fun = mean)

#ggplot(df, aes(y = BONCAT_freq, x = n_species)) +
#  geom_point() +
#  theme_classic(base_size = 15) +
#  stat_smooth(method = "lm", formula = 'y ~ x', se=F,fullrange = T) +
#  facet_wrap(~Rep)

#ggplot(df, aes(y = BONCAT_freq, x = Treatment)) +
#  geom_point() +
#  theme_classic(base_size = 15) +
#  stat_smooth(method = "lm", formula = 'y ~ x', se=F,fullrange = T) +
#  facet_wrap(~Rep)

################we need normalize by the day/rep  #################3
df$Treatment   <- factor(df$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))

### simple linear model
# in simple lm LB and LG are higher than soil.
m1<-lm(BONCAT_freq  ~Treatment, data=df)
summary(m1)

library(lme4)
lm<-lme4::lmer(data=df, BONCAT_freq~Treatment + (1|Date_sorted))
lm
coef(lm)
summary(lm)
anova(lm)

library(lmerTest)
lm<-lmer(data=df, BONCAT_freq~Treatment + (1|Date_sorted))
s<-summary(lm)
anova(lm)
# grab effect size
effectsize<-s$coefficients
effectsize<-data.frame(effectsize)
colnames(effectsize) <- c("effect.size", "std.error", "df", "tvalue", "pvalue") 
head(effectsize)

###### treatment effect ###########

df %>% #filter(Species1!= "Soil")  %>%
  ggplot(aes(x=Treatment, y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))
 
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

