# Mixtures Boncat
# March 2023
#Jennifer Harris

library(readxl)
library(tidyverse)



# import data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")

df<-read_excel("Flow_cyto_log.xlsx", sheet = 2)


# pellet verse filter

df%>% filter(Treatment!="ctl", method!="NA") %>%
  ggplot(aes(x=as.factor(method), y=Percent_BONCAT)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("method")

df%>% filter(Treatment!="ctl", method!="NA") %>%
  ggplot(aes(x=as.factor(method), y=Cells)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("method")


df%>% filter(Treatment!="ctl", method!="NA") %>%
  ggplot(aes(x=as.factor(method), y=Percent_cells)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("method")



df%>% filter(Treatment!="ctl" & method=="pellet") %>%
  ggplot(aes(x=as.factor(Nitrogen), y=Percent_BONCAT)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("Nitrogen")
# not a huge different between nitrogen and no nitrogen
df$Percent_BONCAT
m1<-lm(Percent_BONCAT~Nitrogen, data = df)
summary(m1)
# not significant in a t test



df$Treatment<- factor(df$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

df%>% filter(Treatment!="ctl", method== "filter"  ) %>%
  ggplot(aes(x=Treatment, y=Percent_BONCAT)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10()+
  xlab("Treatment")+
 
# some patterns by treatment
### model
df1<-df%>%filter(Treatment!="ctl", n_species=="1", method=="filter")
m1<-lm(Percent_BONCAT~Brassica, data=df1)
summary(m1)



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


