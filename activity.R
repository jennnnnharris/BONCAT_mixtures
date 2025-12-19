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
# 
# # import data
# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
# df1<-read_excel("Flow_cyto_master.xlsx", sheet = 2)
# 
# # avg the technical reps that are adj for day.
# df1
# df1<-df1%>% group_by( Rep, Group, Treatment,  Pot_ID) %>%
#   summarise(BONCAT_freq = mean(BONCAT_freq_adj), 
#             n_events_cells= round(mean(n_events_cells), digits = 0),
#             n_events_BONCAT= round(mean(success_adj), digits = 0) , 
#             n = n())
# # import data
# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
# df<-read_excel("Flow_cyto_master.xlsx", sheet = 1)
# biggie<-left_join(df1, df)
# write.csv(biggie, "flow_cyto.csv")

####### import data #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
fc <-read.csv("processed_flow_cyto.csv")
head(fc)

#make dates be dates
fc$Date_Sorted<-mdy(fc$Date_Sorted)

# remove outliers
# filter out day were pos ctl didn't work
#fc <- filter(fc, Trt_ID!="B+N6")
#fc<-filter(fc, Date_Sorted != "2023-05-25")
head(fc)

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)

# # make long name var
# fc$long_name <- fc$Treatment
# fc$long_name<-gsub("L", "Legume", fc$long_name)
# fc$long_name<-gsub("G", "Grass", fc$long_name)
# fc$long_name<-gsub("B", "Brassica", fc$long_name)
# fc$long_name<-gsub("LegumeBrassica", "Legume_Brassica", fc$long_name)
# fc$long_name<-gsub("LegumeGrass", "Legume_Grass", fc$long_name)
# fc$long_name<-gsub("GrassBrassica", "Grass_Brassica", fc$long_name)
# fc$long_name 
# fc$long_name<- factor(fc$long_name, levels = c("Soil", "Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))


######we need normalize by the day/rep ###
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
lab<-gsub("LGB", "X", lab)
lab<-gsub("LB", "X", lab)
lab<-gsub("GB", "A", lab)
lab<-gsub("LG", "C", lab)

lab<-gsub("G", "A", lab)
lab<-gsub("L", "A", lab)
lab<-gsub("B", "X", lab)

lab<-gsub("X", "B", lab)

# Treatment emmean     SE  df asymp.LCL asymp.UCL .group
# GB         -2.81 0.0148 Inf     -2.85     -2.77  a    
# G          -2.80 0.0173 Inf     -2.85     -2.76  a    
# L          -2.77 0.0106 Inf     -2.80     -2.74  a    
# LB         -2.71 0.0126 Inf     -2.74     -2.67   b   
# LGB        -2.70 0.0128 Inf     -2.73     -2.66   b   
# B          -2.69 0.0112 Inf     -2.72     -2.66   b   
# LG         -2.42 0.0109 Inf     -2.45     -2.39    c 

p1<-fc  %>%
  filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=boncat_freq, fill = Treatment)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0))+
  labs(title = "A",
       x="",
       y= "percent active cells")+
  geom_text(y=17.6, label = lab, size=5)
p1

#binomial model with percent data##
# make vector of successes and failures
#prop<-fc %>%
#  mutate(boncat_freq = round(boncat_freq, 0)) %>%
#  mutate(n_failures =  100-boncat_freq)

prop<-fc %>%
  mutate(success = n_events_BONCAT) %>%
  mutate(n_failures =  n_events_cells)
y<-cbind(prop$success, prop$n_failures)

#model
m1<-glm(data= prop, y~Treatment +Block, family = binomial)
summary(m1)
#block not significant so we can take block out
m1<-glm(data= prop, y~Treatment, family = binomial)
m1
summary(m1)

# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
emm_object <- emmeans(m1, ~ Treatment)
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




# p1<-fc  %>%
#   filter(Treatment!="Soil") %>%
#   ggplot(aes(x=as.factor(Legume), y=boncat_freq, fill = as.factor(Legume))) +
#   geom_jitter(width = .2, size=2 )+
#   geom_boxplot(alpha=.5, outlier.shape = NA)+
#   scale_fill_manual(values = IBM)+
#   theme_classic(base_size = 16)+
#   theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
#         plot.title = element_text(hjust = 0))+
#   ylab("percent active")
# 
# p1
# m1<-glm(data= prop, y~Legume, family = binomial)
# m1
# summary(m1)


###################################number of cells ########################

# plot
#G    B   GB   LB   LG  LGB    L 
#"a" "ab"  "a" "ab" "ab" "ab"  "b

fc<-fc  %>%   filter(Treatment!="Soil")
lab = as.character(fc$Treatment)
unique(lab)
lab<-gsub("LGB", "AX", lab)
lab<-gsub("GB", "AX", lab)
lab<-gsub("LG", "AX", lab)
lab<-gsub("LB", "AX", lab)
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
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0))+
  labs(title = "B",
       x="",
       y= "active cells/g rhizosphere")+
  geom_text(y=3150, label = lab, size=5)

p2



require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
svg("activity.svg", width=8, height=4)
grid.arrange(p1, p2, ncol=2)
dev.off()

 
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


p1<-fc  %>%
  filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Legume), y=active_cel_per_g, fill = as.factor(Legume))) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_fill_manual(values = IBM)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none",
        plot.title = element_text(hjust = 0))+
  ylab("percent active")

p1
m1<-lm(data= fc, active_cel_per_g~Legume)

m1
summary(m1)



# corr plot activity with functions ####
# weed seed
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed<-read.csv("weed_seed_decay.csv", row.names = 1)
head(weed)
  
  
df<- left_join(weed, fc)
head(df)  

## check distribution
hist(df$active_cel_per_g)
hist(log(df$active_cel_per_g)+1)
hist(df$boncat_freq)
hist(df$pigweed_prop_nongerm)
hist(log(df$foxtail_num_germ))


# plots
plot(log10(df$active_cel_per_g+1), log10(df$foxtail_num_nongerm+1))
plot(df$active_cel_per_g, df$pigweed_prop_nongerm)
plot(df$boncat_freq, df$pigweed_prop_nongerm)
plot(df$boncat_freq, df$foxtail_prop_nongerm)


# model
m1<-lm(log(df$active_cel_per_g+1) ~log(df$foxtail_num_nongerm+1))
summary(m1)

m1<-lm(pigweed_prop_nongerm~df$boncat_freq, data=df)
summary(m1)

m1<-lm(df$foxtail_prop_nongerm~df$active_cel_per_g)
summary(m1)
m1<-lm(foxtail_prop_nongerm~df$boncat_freq, data=df)
summary(m1)

# no clear signals

### nitrogen fixed with microbial activity
fc
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
nfix<-read.csv("Nfix.csv")

colnames(nfix) 
df<-left_join(nfix, fc)

# distributions
hist(df$active_cel_per_g)
hist(df$perc.Ndfa)
hist(df$n_fix_per_legume)


# plots
plot(df$active_cel_per_g, df$perc.Ndfa)
plot(df$boncat_freq, df$perc.Ndfa)
plot(df$active_cel_per_g, df$n_fix_per_legume)
plot(df$boncat_freq, df$n_fix_per_legume)


# frequency boncat
m1<-lm(df$perc.Ndfa~df$boncat_freq, data=df)
summary(m1)
#plot(m1)
m1<-lm(df$n_fix_per_legume~df$boncat_freq, data=df)
summary(m1)
#plot(m1)

# cells boncat
m1<-lm(df$n_fix_per_legume~df$active_cel_per_g)
summary(m1)
m1<-lm(df$perc.Ndfa~df$active_cel_per_g)
summary(m1)

# log transformed
m1<-lm(df$n_fix_per_legume~log(df$active_cel_per_g+1))
summary(m1)
m1<-lm(df$perc.Ndfa~log(df$active_cel_per_g+1))
summary(m1)




# biomass increase verse activity
fc
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("biomass_potlevel.csv")
colnames(biomass) 
df<-left_join(biomass, fc)


# distributions
hist(df$active_cel_per_g)
hist(df$Total.Root.g)
hist(df$Stem.Biomass.g)


# plots

plot(df$active_cel_per_g, df$Stem.Biomass.g)
plot(df$active_cel_per_g, df$Root.Biomass.g)
plot(df$boncat_freq, df$Total.Root.g)
plot(df$boncat_freq, df$Stem.Biomass.g)


# frequency boncat
m1<-lm(df$Stem.Biomass.g~df$boncat_freq, data=df)
summary(m1)
#plot(m1)
m1<-lm(df$Total.Root.g~df$boncat_freq, data=df)
summary(m1) ### marginal trend
#plot(m1)

# cells boncat
m1<-lm(df$Stem.Biomass.g~df$active_cel_per_g)
summary(m1)  # sig
m1<-lm(df$Total.Root.g~df$active_cel_per_g)
summary(m1)

# log transformed
m1<-lm(df$Stem.Biomass.g~log(df$active_cel_per_g+1))
summary(m1)  # sig
m1<-lm(df$Total.Root.g~log(df$active_cel_per_g+1))
summary(m1)


##### difference between predicted and not ##
fc
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("predicted.biomass.csv")
colnames(biomass) 
df<-left_join(biomass, fc)


# distributions
hist(df$active_cel_per_g)
hist(df$boncat_freq)
hist(df$root.difference)
hist(df$shoot.difference)


# plots

plot(df$active_cel_per_g, df$shoot.difference)
plot(df$active_cel_per_g, df$root.difference)
plot(df$boncat_freq, df$shoot.difference)
plot(df$boncat_freq, df$root.difference)


# frequency boncat
m1<-lm(df$shoot.difference~df$boncat_freq, data=df)
summary(m1)
#plot(m1)

m1<-lm(df$root.difference~df$boncat_freq, data=df)
summary(m1) ### marginal sig
#plot(m1)

# cells boncat
m1<-lm(df$shoot.difference~df$active_cel_per_g)
summary(m1)  
m1<-lm(df$root.difference~df$active_cel_per_g)
summary(m1)

# log transformed
m1<-lm(df$shoot.difference~log(df$active_cel_per_g+1))
summary(m1)  # sig
m1<-lm(df$root.difference~log(df$active_cel_per_g+1))
summary(m1)






