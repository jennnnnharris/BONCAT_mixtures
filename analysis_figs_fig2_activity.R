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
  "#a4bdfc",  # "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)




####### import data #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
#setwd("C:/Users/jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto")
fc <-read.csv("processed_flow_cyto.csv")
head(fc)

#make dates be dates
fc$Date_Sorted<-mdy(fc$Date_Sorted)

# remove outliers
# filter out day were pos ctl didn't work
fc <- filter(fc, Trt_ID!="B+N6")
fc<-filter(fc, Date_Sorted != "2023-05-25")
head(fc)

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)


# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_experiment_planning.xlsx")
head(block)  
fc<-left_join(fc, block)
fc
fc<-fc  %>%   filter(Treatment!="Soil")

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


###### percent active ###########

# fc<-fc  %>%   filter(Treatment!="Soil")
# lab = as.character(fc$Treatment)
# lab<-gsub("LGB", "AX", lab)
# lab<-gsub("LB", "AX", lab)
# lab<-gsub("GB", "A", lab)
# lab<-gsub("LG", "X", lab)
# 
# lab<-gsub("G", "AX", lab)
# lab<-gsub("L", "AX", lab)
# lab<-gsub("B", "AX", lab)
# 
# lab<-gsub("X", "B", lab)

# Treatment emmean    SE  df asymp.LCL asymp.UCL .group
# GB         -2.91 0.184 Inf     -3.40     -2.41  a    
# L          -2.88 0.182 Inf     -3.37     -2.39  ab   
# G          -2.69 0.168 Inf     -3.14     -2.24  ab   
# B          -2.67 0.166 Inf     -3.11     -2.22  ab   
# LGB        -2.54 0.157 Inf     -2.96     -2.12  ab   
# LB         -2.49 0.153 Inf     -2.90     -2.08  ab   
# LG         -2.22 0.137 Inf     -2.58     -1.85   b   

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
       y= "percent active cells")
p1




# stats % active ###############
library(lme4)
library(multcompView)
library(multcomp)
library(emmeans)
# filter 
fc<-fc  %>%   filter(Treatment!="Soil")
# quick check
hist(fc$boncat_freq)
fc$boncat_freq


################## linear model 
m1<- lm(boncat_freq ~ Treatment+block, data = fc)
summary(m1)

#block not sig
m1<- lm(boncat_freq ~ Treatment, data = fc)
summary(m1)

# plot model 
plot(m1) 

# anova
anov = aov(boncat_freq ~ Treatment+boncat_freq, data = fc)
summary(anov) # not sig


####################binomial model#
prop<-fc %>%
  mutate(success = n_events_BONCAT) %>%
  mutate(n_failures =  n_events_cells)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= prop, y~Treatment +block, family = binomial)
summary(m1)

#block not significant so we can take block out
model_binom<-glm(data= prop, y~Treatment+block, family = binomial)
model_binom
summary(model_binom)

library(car)

# Run Anova on your beta model
Anova(model_binom, type = "II")
# Get the EMMs for your treatment groups
emm_object <- emmeans(model_binom, ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
# This will perform the pairwise tests on the log-odds scale, 
# apply the sidek adjustment, and assign letters based on the results.

cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  # The Letters argument is optional, but common for CLDs
                  Letters = letters) 


print(cld_result)



#########beta reg ###########
#install.packages("betareg")
library(betareg)

fc$prop_boncat_freq<-fc$boncat_freq/100

# Formula: y ~ x
# Where y_trans is your (0, 1) proportion
model_beta <- betareg(prop_boncat_freq ~ Treatment+block, data = fc)
summary(model_beta)
# model comparison
AIC(model_binom)
AIC(model_beta)
AIC(m1)

# beta regression has the best fit. 
summary(model_beta)

# model output
# Install and load 'car' if you haven't yet
# install.packages("car")
library(car)

# Run Anova on your beta model
Anova(model_beta, type = "II")

# pairwise tests
library(emmeans)
library(multcomp)
# Get the EMMs for your treatment groups
emm_object <- emmeans(model_beta, ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
# This will perform the pairwise tests on the log-odds scale, 
# apply the sidek adjustment, and assign letters based on the results.

cld_result <- cld(emm_object, 
                  adjust = "tukey", 
                  alpha = 0.05,
                  # The Letters argument is optional, but common for CLDs
                  Letters = letters) 


print(cld_result)




###################################number of cells ########################

# plot
#G    B   GB   LB   LG  LGB    L 
#"a" "ab"  "a" "ab" "ab" "ab"  "b

fc<-fc  %>%   filter(Treatment!="Soil")
# lab = as.character(fc$Treatment)
# unique(lab)
# lab<-gsub("LGB", "AX", lab)
# lab<-gsub("GB", "AX", lab)
# lab<-gsub("LG", "AX", lab)
# lab<-gsub("LB", "AX", lab)
# lab<-gsub("B", "AX", lab)
# lab<-gsub("G", "A", lab)
# lab<-gsub("L", "X", lab)
# 
# lab<-gsub("X", "B", lab)
# unique(lab)


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
       y= "active cells/g rhizosphere")
  #geom_text(y=3150, label = lab, size=5)

p2



require(gridExtra)
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
#svg("activity.svg", width=8, height=4)
#grid.arrange(p1, p2, ncol=2)
#dev.off()

 
# stats



######STATS number of cells
library(lme4)
library(multcompView)
library(multcomp)
library(emmeans)

# overall model 
m1<- lm(active_cel_per_g ~ Treatment+block, data = fc)
summary(m1)

#block not sig
m1<- lm(active_cel_per_g ~ Treatment, data = fc)
summary(m1)

# plot model 
plot(m1) #homoscedasticity looks not great

# anova
anov = aov(active_cel_per_g ~ Treatment, data = fc)
summary(anov)

# tukey test
tukey <- TukeyHSD(anov)
print(tukey) # All difference except LG-LB


# emmeans 
emm_object <- emmeans(m1, specs = ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

cld <- multcompLetters4(m1, tukey)
print(cld)
#L   LB   LG  LGB   GB    B    G 
#"a" "ab" "ab" "ab" "ab" "ab"  "b" 



# log transformed Y
hist(fc$active_cel_per_g)
hist(log(fc$active_cel_per_g))
hist(log10(fc$active_cel_per_g))

fc$log_active_cel_per_g = log(fc$active_cel_per_g)


# overall model 
logm1<- lm(log_active_cel_per_g ~ Treatment+block, data = fc)
summary(m1)
#block weak effect - p= 0.05

#plot model
plot(m1) # this looks alot better

# anova
anov = aov(log_active_cel_per_g ~ Treatment+block, data = fc)
summary(anov)

# tukey test
#tukey <- TukeyHSD(anov)
#print(tukey) # All difference except LG-LB
#plot(m1) #homoscedasticity looks not great

# emmeans 
emm_object <- emmeans(m1, specs = ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

cld <- multcompLetters4(m1, tukey)
print(cld)
#L   LB   LG  LGB   GB    B    G 
#"a" "ab" "ab" "ab" "ab" "ab"  "b" 



# coded differently

# overall model 
m1<- lm(log_active_cel_per_g ~ Legume+Grass+block, data = fc)
summary(m1)


# monocultures
fc1<-fc %>% filter(n_species==1)
a1<- aov(active_cel_per_g~ Treatment, data = fc1)
summary(a1)
tukey_results <- TukeyHSD(a1)
print(tukey_results)
p_values <- tukey_results$Treatment[, 4]#Extract the p-values for the factor of interest
cld <- multcompLetters(p_values)# The 'multcompLetters()' function takes a named vector of p-values
print(cld)#Print the results



#legume effect
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
summary(m1)





###################predicting microbial activity with biomass ##############

#load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(lubridate)

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("percent.biomass.csv")
# only n+ 
biomass<-biomass %>% filter(N==1)
biomass




# LG
# filter fc for L treatment
fc %>% filter(Treatment=="L") 
fc %>% dplyr::select(boncat_freq)
df<-fc %>% filter(Treatment=="L") %>% dplyr::select(boncat_freq, active_cel_per_g )
df
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LG") %>% filter(Species=="legume") %>% arrange(Trt_ID)
sp1
vector<-sp1$percent/100
data_frame_multiplied1 <- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
data_frame_multiplied1
# filter for G treat
fc %>% filter(Treatment=="G") 
df<-fc %>% filter(Treatment=="G") %>% dplyr::select(boncat_freq, active_cel_per_g )


# get biomass info and multiply by G biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LG") %>% filter(Species=="grass") %>% arrange(Trt_ID)
sp1
vector<-sp1$percent/100
vector

data_frame_multiplied2<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )

print(data_frame_multiplied2)
# add together
LG<-data_frame_multiplied1 + data_frame_multiplied2
LG$Trt_ID<-c("predict_LG+N1", "predict_LG+N2", "predict_LG+N3", "predict_LG+N4", "predict_LG+N5", "predict_LG+N6" )
LG$Treatment<-c("LG.predict", "LG.predict", "LG.predict", "LG.predict", "LG.predict", "LG.predict" )
LG$Rep <- c(1, 2, 3, 4, 5, 6 )
LG


# LB
# filter fc for L treatment
# check
fc %>% filter(Treatment=="L") 
df<-fc %>% filter(Treatment=="L")%>% filter(Rep!="6" & Rep!="4") %>% dplyr::select(boncat_freq, active_cel_per_g ) 
df
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LB") %>% filter(Species=="legume") %>% arrange(Trt_ID) %>% filter( Rep!="6") 
sp1
vector<-sp1$percent/100
vector
data_frame_multiplied1 <- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
data_frame_multiplied1
# filter for B treat
fc %>% filter(Treatment=="B")
df<-fc %>% filter(Treatment=="B") %>%  dplyr::select(boncat_freq, active_cel_per_g )
df
# get biomass info and multiply by B biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LB") %>% filter(Species=="brassica") %>%  filter(Rep!="6") %>% arrange(Trt_ID)
sp1
vector<-sp1$percent/100
vector

data_frame_multiplied2<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )

print(data_frame_multiplied2)
# add together
LB<-data_frame_multiplied1 + data_frame_multiplied2
LB$Trt_ID<-c("predict_LB+N1", "predict_LB+N2", "predict_LB+N3", "predict_LB+N4" )
LB$Treatment<-c("LB.predict", "LB.predict", "LB.predict", "LB.predict")
LB$Rep <- c(1, 2, 3, 4 )
LB




# GB
# filter fc for G treatment
fc %>% filter(Treatment=="G") 
df<-fc %>% filter(Treatment=="G")  %>% filter(Rep!="6"  & Rep!="4" ) %>%dplyr::select(boncat_freq, active_cel_per_g ) 
df
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="GB") %>% filter(Species=="grass") %>% arrange(Trt_ID)  %>% filter(Rep!="6"  & Rep!="4" )
sp1
vector<-sp1$percent/100
vector
data_frame_multiplied1 <- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
data_frame_multiplied1
# filter for B treat
fc %>% filter(Treatment=="B")
df<-fc %>% filter(Treatment=="B") %>%  dplyr::select(boncat_freq, active_cel_per_g )
df
# get biomass info and multiply by B biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="GB") %>% filter(Species=="brassica") %>% arrange(Trt_ID)  %>% filter(Rep!="6" & Rep!="4")
sp1
vector<-sp1$percent/100
vector

data_frame_multiplied2<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )

print(data_frame_multiplied2)
# add together
GB<-data_frame_multiplied1 + data_frame_multiplied2
GB$Trt_ID<-c("predict_GB+N1", "predict_GB+N2", "predict_GB+N3", "predict_GB+N5" )
GB$Treatment<-c("GB.predict", "GB.predict", "GB.predict", "GB.predict" )
GB$Rep <- c(1, 2, 3, 5)
GB





library(dplyr)
# LGB
# filter fc for L treatment
fc %>% filter(Treatment=="L") 
df<-fc %>% filter(Treatment=="L") %>%  filter(Rep!="6"& Rep!="4" )  %>% dplyr::select(boncat_freq, active_cel_per_g )
df
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LGB") %>% filter(Species=="legume") %>% arrange(Trt_ID) %>% filter(Rep!="6"  & Rep!="4" ) 
sp1

vector<-sp1$percent/100
vector
data_frame_multiplied0 <- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
data_frame_multiplied0

# filter fc for G treatment
fc %>% filter(Treatment=="G") 
df<-fc %>% filter(Treatment=="G") %>%  filter(Rep!="6"& Rep!="4" ) %>%    dplyr::select(boncat_freq, active_cel_per_g ) 
df
# get biomass info and multiply by biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LGB") %>% filter(Species=="grass") %>% arrange(Trt_ID) %>%  filter(Rep!="6"& Rep!="4" )
sp1
vector<-sp1$percent/100
vector
data_frame_multiplied1 <- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
data_frame_multiplied1

# filter for B treat
fc %>% filter(Treatment=="B")
df<-fc %>% filter(Treatment=="B") %>%    dplyr::select(boncat_freq, active_cel_per_g )
df
# get biomass info and multiply by G biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LGB") %>% filter(Species=="brassica") %>% arrange(Trt_ID) %>%  filter(Rep!="6"& Rep!="4" )
sp1
vector<-sp1$percent/100
vector

data_frame_multiplied2<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )

print(data_frame_multiplied2)
# add together
LGB<-data_frame_multiplied0 + data_frame_multiplied1 + data_frame_multiplied2
LGB$Trt_ID<-c("predict_LGB+N1", "predict_LGB+N2", "predict_LGB+N3", "predict_LGB+N4")
LGB$Treatment<-c("LGB.predict", "LGB.predict", "LGB.predict", "LGB.predict")
LGB$Rep <- c(1, 2, 3, 4)
LGB


# combine
predict<-rbind(GB, LG, LB, LGB)
df<-full_join(fc, predict)



mycols <- c( #IBM colors
  "navy", # dark royal blue L
 "#648FFF", # french blue G
  "#785EF0", # light purple B
  "#DC267F", # magenta pink GB
  "grey",
  "#FE6100", # bright orange LB
  "grey",
  "#FFB000", # golden yellow LG
  "grey",
  "#865338", # medium mocha brown LGB
  "grey"
)


df<-df %>% filter(Treatment!="Soil")
#df<-df %>% filter(Treatment!="G" & Treatment!="B" & Treatment!="L")

unique(df$Treatment)

df$Treatment<-factor(df$Treatment, levels = c("L", "G", "B", "GB", "GB.predict", "LB", "LB.predict",  "LG", 
                                                "LG.predict", "LGB", "LGB.predict" ))

df
# plot
p1<-df  %>% 
  ggplot(aes(x=Treatment, y=active_cel_per_g, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none")+
  labs(title = "A",
       x="",
       y="active cells/g rhizosphere")
  #geom_text(y=.25, label =label , nudge_x = -.8, size=8)

p1


p2<-df  %>% 
  ggplot(aes(x=Treatment, y=boncat_freq, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none")+
labs(title = "B",
     x="",
     y=" active cells (%)")
#geom_text(y=.25, label =label , nudge_x = -.8, size=8)

p2
require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
svg("activity.predict.svg", width=10, height=3.5)
grid.arrange(p1, p2, ncol=2)
dev.off()

# stats #####
# anova active cells
#GB
df1<-df%>% filter(Treatment=="GB" | Treatment=="GB.predict")
m1<-lm(data= df1, active_cel_per_g~Treatment)
summary(m1)
m<-summary(m1)
m$coefficients

#LB
df1<-df%>% filter(Treatment=="LB" | Treatment=="LB.predict")
m1<-lm(data= df1, active_cel_per_g~Treatment)
summary(m1)
m<-summary(m1)
m$coefficients


#LG
df1<-df%>% filter(Treatment=="LG" | Treatment=="LG.predict")
m1<-lm(data= df1, active_cel_per_g~Treatment)
summary(m1)
summary(m1)
m<-summary(m1)
m$coefficients

#LGB
df1<-df%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
m1<-lm(data= df1, active_cel_per_g~Treatment)
summary(m1)




# anova binomial model percent 
#binomial model with percent data##
# make vector of successes and failures



#GB
df1<-df%>% filter(Treatment=="GB" | Treatment=="GB.predict")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment, family = binomial)
summary(m1)

#LB
df1<-df%>% filter(Treatment=="LB" | Treatment=="LB.predict")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment, family = binomial)
summary(m1)

#LG
df1<-df%>% filter(Treatment=="LG" | Treatment=="LG.predict")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment, family = binomial)
summary(m1)


#LGB
df1<-df%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment, family = binomial)
summary(m1)

######## load data############


################predicting microbial activity with simple proportions##########
##### import data 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
#setwd("C:/Users/jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto")
df <-read.csv("processed_flow_cyto.csv")
head(df)

#make dates be dates
df$Date_Sorted<-mdy(df$Date_Sorted)

# remove outliers
# filter out day were pos ctl didn't work
df <- filter(df, Trt_ID!="B+N6")
df<-filter(df, Date_Sorted != "2023-05-25")
head(df)

# make treatment and day factors
#df$Treatment   <- factor(df$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
df$Date_Sorted   <- factor(df$Date_Sorted)


# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_experiment_planning.xlsx")
head(block)  
df<-left_join(df, block)
df
df<- df  %>%   filter(Treatment!="Soil")
df$Treatment

get.predict<-function(df, sp1, sp2, sp3) {
  if(missing(sp3)){
    sp1.df<- filter(df, Treatment==sp1) %>% select(boncat_freq)
    sp1.predict<-sp1.df$boncat_freq/2
    sp2.df<-filter(df, Treatment==sp2)
    sp2.predict<-sp2.df$boncat_freq/2
    predict <- sp1.predict + sp2.predict
    print(predict)
    return(predict)
  }
  else
    
    print(paste("the 3rd species is",sp3))
  sp1.df<- filter(df, Treatment==sp1) %>% select(boncat_freq)
  sp1.predict<-sp1.df$boncat_freq/3
  sp2.df<-filter(df, Treatment==sp2)
  sp2.predict<-sp2.df$boncat_freq/3
  sp3.df<-filter(df, Treatment==sp3)
  sp3.predict<-sp3.df$boncat_freq/3
  predict <- sp1.predict + sp2.predict + sp3.predict
  print(predict)
  return(predict)
  
}




## we expect that a plant makes the same amount of biomass in monoculures vs mixtures
#example
#LB biomass = L monoculture/  2 + B monocultre/2 

#get shoot predictions

filter(df, df$Treatment=="G")
filter(df, df$Treatment=="B")

GB<-get.predict(df, "G", "B") 
LB<-get.predict(df, "L", "B") 
LG<-get.predict(df, "L", "G") 
LGB<-get.predict(df, "L", "G", "B")  
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
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB", "GB.predict", "LB", "LB.predict",  "LG", 
                                                "LG.predict", "LGB", "LGB.predict" ))

#### plot predi
