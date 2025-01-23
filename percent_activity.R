# Mixtures Boncat
# March 2023
# last edited: Jan 2025
#Jennifer Harris

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
df%>% #filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq)) +
  geom_jitter(width = .2, size=2 )+
  theme_bw(base_size = 18, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  stat_summary(geom = "line", fun = mean)

# increasing species line w/o brasssisae 
df%>% #filter(Species1!="Soil") %>%
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

df%>% filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq, col=as.factor(Rep))) +
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



############################################treatment ###########

df$Treatment   <- factor(df$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))

df %>% #filter(Species1!= "Soil")  %>%
  ggplot(aes(x=Treatment, y=BONCAT_freq)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))
 

### model
# in simple lm LB and LG are higher than soil.
m1<-lm(BONCAT_freq  ~Treatment, data=df)
summary(m1)

#make summary table
avgs<-df %>% group_by(Treatment) %>%
summarise(mean = mean(BONCAT_freq), sd = sd(BONCAT_freq))

#expectations
# G + B
avgs[which(avgs$Treatment == "B"), 2] + avgs[which(avgs$Treatment == "G"), 2] # mean
sqrt(avgs[which(avgs$Treatment == "B"), 3] + avgs[which(avgs$Treatment == "G"), 3])# standard dev 
# 12.4 + or  - 3.03 
#actually
avgs[which(avgs$Treatment == "GB"), 2] 
# 6.7 + or - 5.3


# L + B
avgs[which(avgs$Treatment == "B"), 2] + avgs[which(avgs$Treatment == "L"), 2] #mean
sqrt(avgs[which(avgs$Treatment == "B"), 3] + avgs[which(avgs$Treatment == "G"), 3]) #stdev
# expectation:  12.2 +/- # 3.06
avgs[which(avgs$Treatment == "LB"), 2] #mean
avgs[which(avgs$Treatment == "LB"), 3]  #stdev
# actual : 10.5 +/- 7


# L + G
avgs[which(avgs$Treatment == "G"), 2] + avgs[which(avgs$Treatment == "L"), 2] #mean
sqrt(avgs[which(avgs$Treatment == "L"), 3] + avgs[which(avgs$Treatment == "G"), 3]) #stdev
# expectation:  8.6 +/- # 2.5
avgs[which(avgs$Treatment == "LG"), 2] #mean
avgs[which(avgs$Treatment == "LG"), 3]  #stdev
# actual : 12.8 +/- 8.5

# L + G + B
avgs[which(avgs$Treatment == "G"), 2] + avgs[which(avgs$Treatment == "L"), 2] + avgs[which(avgs$Treatment == "B"), 2] #mean
sqrt(avgs[which(avgs$Treatment == "L"), 3] + avgs[which(avgs$Treatment == "G"), 3] +  avgs[which(avgs$Treatment == "B"), 3] )#stdev
# expectation:  16.6 +/- 3.5
avgs[which(avgs$Treatment == "LGB"), 2] #mean
avgs[which(avgs$Treatment == "LGB"), 3]  #stdev
# actual : 7.8 +/- 7.8


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


