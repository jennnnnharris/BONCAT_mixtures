# Mixtures Boncat
# created: March 2023
# last edited: July 2026
# author: Jennifer Harris

rm(list=ls())
rstudioapi::restartSession(clean = TRUE)

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
fc <-read.csv("processed_flowcyto.csv")
head(fc)

#make dates be dates
fc$Date_Sorted<-mdy(fc$Date_Sorted)

# remove outliers
# filter out day were pos ctl didn't work
fc <- filter(fc, Trt_ID!="B+N6")
#fc<-filter(fc, Date_Sorted != "2023-05-25")
head(fc)

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)


# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_blockinfo.xlsx")
head(block)  
fc<-left_join(fc, block)
fc
fc<-fc  %>%   filter(Treatment!="Soil")


###### percent active ###########

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


####################binomial model#
prop<-fc %>%
  mutate(success = n_events_BONCAT) %>%
  mutate(n_failures =  n_events_cells)
y<-cbind(prop$success, prop$n_failures)
model_binom<-glm(data= prop, y~Treatment +block, family = binomial)
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


# legume effect
model_binom<-glm(data= prop, y~Legume +block, family = binomial)
summary(model_binom)
library(car)
Anova(model_binom, type = "II")





###################################number of cells ########################

fc<-fc  %>%   filter(Treatment!="Soil")
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
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
svg("activity.svg", width=8, height=4)
grid.arrange(p1, p2, ncol=2)
dev.off()


######stats number of cells############
library(lme4)
library(multcompView)
library(multcomp)
library(emmeans)
install.packages("MASS")
library(MASS)

# overall model 
m1<- lm(active_cel_per_g ~ Treatment+block, data = fc)
summary(m1)
plot(m1) #homoscedasticity looks not great
# anova
anov = aov(active_cel_per_g ~ Treatment+block, data = fc)
summary(anov)
# emmeans 
emm_object <- emmeans(m1, specs = ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
# get letters
cld <- multcompLetters4(m1, tukey)
print(cld)


# log transformed Y model
hist(fc$active_cel_per_g)
hist(log(fc$active_cel_per_g))
hist(log10(fc$active_cel_per_g))
fc$log_active_cel_per_g = log(fc$active_cel_per_g)
# overall model 
m1<- aov(log_active_cel_per_g ~ Treatment+block, data = fc)
summary(m1)
#plot model
plot(m1) # this looks alot better
# anova
anov = aov(log_active_cel_per_g ~ Treatment+block, data = fc)
summary(anov)
# emmeans 
emm_object <- emmeans(anov, specs = ~ Treatment)
# Perform all pairwise comparisons with Tukey adjustment
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)
cld <- multcompLetters4(m1, tukey)
print(cld)


# piosson 
m1  <- glm(active_cel_per_g ~ Treatment+block, data = fc, family = poisson(link = "log"))
plot(m1)
anova(m1)


# negative binomial
model_nb <- glm.nb(active_cel_per_g ~ Treatment+block, data = fc)
plot(model_nb)
anova(m1)

emm_object <- emmeans(model_nb, specs = ~ Treatment)
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)



#legume effect
model_nb <- glm.nb(active_cel_per_g ~ Legume+block, data = fc)
summary(model_nb)

#legume effect plot#
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




###################predicting microbial activity with biomass ##############

# restart R for package conflicts
rm(list=ls())

#load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(lubridate)

# load data
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
fc <-read.csv("processed_flowcyto.csv")
fc$Date_Sorted<-mdy(fc$Date_Sorted)

# filter out day were pos ctl didn't work
fc <- filter(fc, Trt_ID!="B+N6")
fc<-filter(fc, Date_Sorted != "2023-05-25")

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)

# add block info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
block <- read_excel("metadata_blockinfo.xlsx")
fc<-left_join(fc, block)
fc<-fc  %>%   filter(Treatment!="Soil")

# load biomass info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/biomass")
biomass<-read.csv("biomass_percent.csv")
# only n+ 
biomass<-biomass %>% filter(N==1)
biomass




#### predict values#########
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
LG$Treatment<-c("LG_expected", "LG_expected", "LG_expected", "LG_expected", "LG_expected", "LG_expected" )
LG$Rep <- c(1, 2, 3, 4, 5, 6 )
LG


# LB
# filter fc for L treatment
fc %>% filter(Treatment=="L") 
df<-fc %>% filter(Treatment=="L")%>% filter(Rep!="6" & Rep!="4") %>% dplyr::select(boncat_freq, active_cel_per_g ) 
df
#L has 1-5
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LB") %>% filter(Species=="legume") %>% arrange(Trt_ID) %>% filter( Rep!="6"  & Rep!="4") 
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
#B has 1-4
# get biomass info and multiply by B biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LB") %>% filter(Species=="brassica") %>%  filter(Rep!="6" & Rep!="4") %>% arrange(Trt_ID)
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
LB$Treatment<-c("LB_expected", "LB_expected", "LB_expected", "LB_expected")
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
GB$Treatment<-c("GB_expected", "GB_expected", "GB_expected", "GB_expected" )
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
LGB$Treatment<-c("LGB_expected", "LGB_expected", "LGB_expected", "LGB_expected")
LGB$Rep <- c(1, 2, 3, 4)
LGB


# combine
predict<-rbind(GB, LG, LB, LGB)
df<-full_join(fc, predict)
head(df)


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
"#865338"# medium mocha brown LGB
)


df<-df %>% filter(Treatment!="Soil")
#df<-df %>% filter(Treatment!="G" & Treatment!="B" & Treatment!="L")

unique(df$Treatment)

df$Treatment<-factor(df$Treatment, levels = c("L", "G", "B", "GB_expected", "GB","LB_expected", "LB",  "LG_expected",  "LG", 
                                              "LGB_expected" , "LGB"  ))

df
# plot
p1<-df  %>% 
  ggplot(aes(x=Treatment, y=active_cel_per_g, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 14)+
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
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none")+
labs(title = "B",
     x="",
     y=" active cells (%)")
#geom_text(y=.25, label =label , nudge_x = -.8, size=8)

p2
require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_04active")
svg("activity.predict.svg", width=10, height=5)
grid.arrange(p1, p2, ncol=2)
dev.off()

# stats #####
# # linear model  active cells
# #GB
# df1<-df%>% filter(Treatment=="GB" | Treatment=="GB_expected")
# m1<-lm(data= df1, active_cel_per_g~Treatment)
# summary(m1)
# m<-summary(m1)
# m$coefficients
# 
# #LB
# df1<-df%>% filter(Treatment=="LB" | Treatment=="LB_expected")
# m1<-lm(data= df1, active_cel_per_g~Treatment)
# summary(m1)
# m<-summary(m1)
# m$coefficients
# 
# 
# #LG
# df1<-df%>% filter(Treatment=="LG" | Treatment=="LG_expected")
# m1<-lm(data= df1, active_cel_per_g~Treatment)
# summary(m1)
# m<-summary(m1)
# m$coefficients
# 
# #LGB
# df1<-df%>% filter(Treatment=="LGB" | Treatment=="LGB_expected")
# m1<-lm(data= df1, active_cel_per_g~Treatment)
# summary(m1)


# beta regression for active cells #
df$Treatment<-as.character(df$Treatment)
library(MASS)


#GB
df
df1<-df%>% filter(Treatment=="GB" | Treatment=="GB_expected")
df1
m1<- glm.nb(active_cel_per_g ~ Treatment+Rep, data = df1)
anova(m1)
summary(m1)

#LB
df1<-df%>% filter(Treatment=="LB" | Treatment=="LB_expected")
m1<- glm.nb(active_cel_per_g ~ Treatment+Rep, data = df1)
anova(m1)

#LG
df1<-df%>% filter(Treatment=="LG" | Treatment=="LG_expected")
m1<- glm.nb(active_cel_per_g ~ Treatment+Rep, data = df1)
anova(m1)


#LGB
df1<-df%>% filter(Treatment=="LGB" | Treatment=="LGB_expected")
m1<- glm.nb(active_cel_per_g ~ Treatment+Rep, data = df1)
anova(m1)


# extract table
library(dplyr)
library(purrr)
library(MASS)
library(broom)
# Vector of treatment groups to test
treatments <- c("GB", "LB", "LG", "LGB")

# 1. Extract Model Coefficients (Estimate, Std. Error, z-value, p-value)
coefficients_table <- map_dfr(treatments, function(trt) {
  df_sub <- df %>% 
    filter(Treatment %in% c(trt, paste0(trt, "_expected")))
  
  model <- glm.nb(active_cel_per_g ~ Treatment + Rep, data = df_sub)
  
  tidy(model) %>%
    mutate(Comparison = trt, .before = 1)
})

# 2. Extract ANOVA / Deviance Test Results (Likelihood Ratio / Deviance Table)
anova_table <- map_dfr(treatments, function(trt) {
  df_sub <- df %>% 
    filter(Treatment %in% c(trt, paste0(trt, "_expected")))
  
  model <- glm.nb(active_cel_per_g ~ Treatment + Rep, data = df_sub)
  
  anova(model) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    mutate(Comparison = trt, .before = 1)
})

# View results
coefficients_table %>% filter(term!="Rep" & term!="(Intercept)")
print(coefficients_table)
print(anova_table)


#binomial model with percent data##
#GB
df1<-df%>% filter(Treatment=="GB" | Treatment=="GB_expected")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment+Rep, family = binomial)
summary(m1)

#LB
df1<-df%>% filter(Treatment=="LB" | Treatment=="LB_expected")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment+Rep, family = binomial)
summary(m1)

#LG
df1<-df%>% filter(Treatment=="LG" | Treatment=="LG_expected")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment+Rep, family = binomial)
summary(m1)


#LGB
df1<-df%>% filter(Treatment=="LGB" | Treatment=="LGB_expected")
prop<-df1 %>%
  mutate(success = round(boncat_freq, 0)) %>%
  mutate(n_failures =  100-success)
y<-cbind(prop$success, prop$n_failures)
m1<-glm(data= df1, y~Treatment+Rep, family = binomial)
summary(m1)




