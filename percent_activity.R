# Mixtures Boncat
# created: March 2023
# last edited: April 25
# author: Jennifer Harris

#rstudioapi::restartSession(clean = TRUE)
rm(list=ls())

#load libraries 
library(readxl)
library(tidyverse)
library(lubridate)



# set colors
IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)



### clean data  #############

# import data
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
df1<-read_excel("Flow_cyto_master.xlsx", sheet = 2)

# avg the technical reps that are adj for day.
df1
df1<-df1%>% group_by( Rep, Group, Treatment,  Pot_ID) %>%
  summarise(BONCAT_freq = mean(BONCAT_freq_adj), 
            n_events_cells= round(mean(n_events_cells), digits = 0),
            n_events_BONCAT= round(mean(success_adj), digits = 0) , 
            n = n())
# import data
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
df<-read_excel("Flow_cyto_master.xlsx", sheet = 1)
colnames(df)
head(df)
biggie<-left_join(df1, df)
write.csv(biggie, "flow_cyto.csv")

############################################
# import data
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
fc <-read.csv("flow_cyto.csv")
head(fc)


#make dates be dates
fc$Date_Sorted<-ymd(fc$Date_Sorted)

# filter out day were pos ctl didn't work
fc<-filter(fc, Date_Sorted != "2023-05-25")
head(fc)

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)

# make long name var
fc$long_name <- fc$Treatment
fc$long_name<-gsub("L", "Legume", fc$long_name)
fc$long_name<-gsub("G", "Grass", fc$long_name)
fc$long_name<-gsub("B", "Brassica", fc$long_name)
fc$long_name<-gsub("LegumeBrassica", "Legume_Brassica", fc$long_name)
fc$long_name<-gsub("LegumeGrass", "Legume_Grass", fc$long_name)
fc$long_name<-gsub("GrassBrassica", "Grass_Brassica", fc$long_name)
fc$long_name 
fc$long_name<- factor(fc$long_name, levels = c("Soil", "Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))



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
#df%>% 
#  ggplot(aes(x=Date_sorted, y=BONCAT_freq)) +
#  geom_jitter(width = .2, size=1 )+
#  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
#  theme_bw(base_size = 18, )+
#  theme(axis.text.x = element_text(angle=60, hjust=1))

# adjusted plot
#df%>% 
#  ggplot(aes(x=Date_sorted, y=BONCAT_freq_adj)) +
#  geom_jitter(width = .2, size=1 )+
#  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
#  theme_bw(base_size = 18, )+
#  theme(axis.text.x = element_text(angle=60, hjust=1))


###### plot ###########

fc<-fc  %>%   filter(Treatment!="Soil")
lab = as.character(fc$Treatment)
lab<-gsub("LGB", "A", lab)
lab<-gsub("GB", "A", lab)
lab<-gsub("B", "A", lab)
lab<-gsub("LA", "B", lab)
lab<-gsub("LG", "B", lab)
lab

lab<-gsub("G", "A", lab)
lab<-gsub("L", "A", lab)

p1<-fc  %>%
  filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=BONCAT_freq, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0))+
  ylab("percent active")+
  xlab("")+
  #facet_grid( ~n_species, scales = "free", space = "free")+
  geom_text(y=15, label = lab, size=5)
p1

#binomial model with percent data##
# make vector of successes and failures
prop<-fc %>%
  mutate(BONCAT_freq = round(BONCAT_freq, 0)) %>%
  mutate(n_failures =  100-BONCAT_freq)
y<-cbind(prop$BONCAT_freq, prop$n_failures)
y

#model
m1<-glm(data= prop, y~Treatment +Block, family = binomial)
summary(m1)
#block not significant so we can take block out

# model
m1<-glm(data= prop, y~Treatment, family = binomial)
m1
summary(m1)



###################################number of cells ########################

# remove outlier
fc <- filter(fc, Trt_ID!="B+N6")
head(fc)

# plot


#G    B   GB   LB   LG  LGB    L 
#"a" "ab"  "a" "ab" "ab" "ab"  "b

fc<-fc  %>%   filter(Treatment!="Soil")
lab = as.character(fc$Treatment)
unique(lab)
lab<-gsub("LGB", "AX", lab)
lab<-gsub("GB", "A", lab)
lab<-gsub("LG", "AX", lab)
lab<-gsub("LB", "AX", lab)

unique(lab)
lab<-gsub("B", "AX", lab)
lab<-gsub("G", "A", lab)
lab<-gsub("L", "X", lab)

lab<-gsub("X", "B", lab)
unique(lab)


p2<-fc  %>%
  filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=active_cel_per_g, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0))+
  #ylab("percent active")+
  xlab("")+
  geom_text(y=3000, label = lab, size=5)

#facet_grid( ~n_species, scaactive_cel_per_g#facet_grid( ~n_species, scales = "free", space = "free")+
p2




#scale_shape_discrete() 
require(gridExtra)
grid.arrange(p1, p2, ncol=2)

 
# stats
library(multcompView)
a1<- aov(active_cel_per_g~ Treatment, data = fc)
summary(a1)

#tukey test
tukey_results <- TukeyHSD(a1)
print(tukey_results)
p_values <- tukey_results$Treatment[, 4]#Extract the p-values for the factor of interest
cld <- multcompLetters(p_values)# The 'multcompLetters()' function takes a named vector of p-values
print(cld)#Print the results






# corr plot activity with functions ####

# import data
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
fc <-read.csv("flow_cyto.csv")
head(fc)

#make dates be dates
fc$Date_Sorted<-ymd(fc$Date_Sorted)
# filter out day were pos ctl didn't work
#fc<-filter(fc, Date_Sorted != "2023-05-25") 
#head(fc)
# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)


###### clean data weed seed
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed <- read_csv("GH Mix Germination Data(Sheet1).csv")
head(weed)


f<-weed %>% filter(Species=="foxtail")
# tidy data  

colnames(f) <- c("Trt_ID" ,
                 "Treatment",
                 "Rep",
                 "Block",
                 "Species",
                 "foxtail_num_nongerm",
                 "foxtail_num_germ",
                 "foxtail_total_seeds",
                 "foxtail_prop_nongerm")
f<-f%>% select(-Species)
## pigweed


p<-weed %>% filter(Species=="pigweed")
# tidy data  

colnames(p) <- c("Trt_ID" ,
                 "Treatment",
                 "Rep",
                 "Block",
                 "Species",
                 "pigweed_num_nongerm",
                 "pigweed_num_germ",
                 "pigweed_total_seeds",
                 "pigweed_prop_nongerm")
p<-p%>% select(-Species)
weed<-left_join(p,f)
write.csv(weed, "weed_seed_decay.csv")


weed<-read.csv("weed_seed_decay.csv", row.names = 1)
  head(weed)
  
  
df<- left_join(weed, fc)
head(df)  

##

plot(df$active_cel_per_g, df$foxtail_prop_nongerm)
plot(df$BONCAT_freq, df$pigweed_prop_nongerm)

m1<-lm(df$pigweed_prop_nongerm~df$active_cel_per_g)
summary(m1)

m1<-lm(pigweed_prop_nongerm~, data=df)
summary(m1)



### nitrogen fixed with microbial activity



# biomass increase verse activity
