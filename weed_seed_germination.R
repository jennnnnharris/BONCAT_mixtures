# weed seed germination

######## weed seed decay #######
#load libraries
library(readxl)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)

IBM <- c( #IBM colors
  "navy", # dark royal blue L
  "#648FFF", # french blue G
  "#785EF0", # light purple B
  "#DC267F", # magenta pink GB
  "#FE6100", # bright orange LB
  "#FFB000", # golden yellow LG
  "#865338" # medium mocha brown LGB
)

# load data 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df <- read.csv("weed_seed_decay.csv", row.names = 1 )
head(df)

# set factor
df$Treatment   <- factor(df$Treatment, levels= c("S", "L", "G", "B", "GB", "LB", "LG", "LGB"))
unique(df$Treatment)


# foxtail
#L G B LB GB LG LGB
# a a a ab ab ab b
df1<-df  %>%   filter(Treatment!="S") 
df1$Treatment   <- factor(df1$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

lab = as.character(df1$Treatment)
unique(lab)
lab<-gsub("LGB", "X", lab)
lab<-gsub("GB", "AX", lab)
lab<-gsub("LG", "AX", lab)
lab<-gsub("LB", "AX", lab)
lab<-gsub("B", "A", lab)
lab<-gsub("G", "A", lab)
lab<-gsub("L", "A", lab)

lab<-gsub("X", "B", lab)
unique(lab)
length(lab)


p1<- df1 %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_nongerm, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "none")+
  scale_fill_manual(values = IBM)+
  geom_text(y=.21, label = lab, size=3)+
  labs(title = "A",
       x="",
       y="non germinating foxtail (%)")
p1





# pigweed
# L   G   B   GB    LB  LG  LGB
# bc abc  d   ab    a   abc cd

lab = as.character(df1$Treatment)
unique(lab)
lab<-gsub("LGB", "CD", lab)
lab<-gsub("GB", "AX", lab)
lab<-gsub("LG", "AXC", lab)
lab<-gsub("LB", "A", lab)
lab<-gsub("B", "D", lab)
lab<-gsub("G", "ABC", lab)
lab<-gsub("L", "XC", lab)
lab<-gsub("X", "B", lab)

p2<-df1 %>%
  ggplot(aes(x=Treatment, y=pigweed_prop_nongerm, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "none")+
  scale_fill_manual(values = IBM)+
  geom_text(y=.73, label = lab, size=3)+
  labs(title = "B",
       x="",
       y="non germinating pigweed (%)")
p2


require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
svg("weeddecay.svg", width = 5, height = 2.5)
grid.arrange(p1, p2, ncol=2)
dev.off()




# model - binomial model with percent data##
# make vector of successes and failures
y<-cbind(df$num_nongerm, df$num_germ)
m1<-glm(data= df, y~(treatment+block+species)^2, family = binomial)
summary(m1)
# significant effect of block and species

# model without soil
df <- df %>% filter(treatment!="S")
y<-cbind(df$num_nongerm, df$num_germ)
m1<-glm(data= df, y~(treatment+block+species)^2, family = binomial)
summary(m1)
# significant effect of block and species and species treatment interaction


# model foxtail
df1 <- df %>% filter(species=="foxtail")
y<-cbind(df1$num_nongerm, df1$num_germ)
m1<-glm(data= df1, y~treatment+block, family = binomial)
summary(m1)
# significant effect of block and species and species treatment interaction
# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
emm_object <- emmeans(m1, specs = ~ treatment)
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



# model pigweed
df2 <- df %>% filter(species=="pigweed")
y<-cbind(df2$num_nongerm, df2$num_germ)
m1<-glm(data= df2, y~treatment+block, family = binomial)
summary(m1)

# significant effect of block and species and species treatment interaction
# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
emm_object <- emmeans(m1, specs = ~ treatment)
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

# L   G   B   GB    LB  LG  LGB
# bc abc  d   ab    a   abc cd

########predicted

######predictions weed seed########
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
# group by
df<- df %>% group_by(Rep, Nitrogen)
df

# pigdf
LG<-get.predict(df, "L", "G", "pigweed_prop_nongerm")
LB<-get.predict(df, "L", "B", "pigweed_prop_nongerm")
GB<-get.predict(df, "G", "B", "pigweed_prop_nongerm")
LGB<-get.predict(df, "L", "G", "pigweed_prop_nongerm", "B")
LG

predict<-rbind(GB, LB, LG, LGB)
predict$pigweed_prop_nongerm<-round(as.numeric(predict$pigweed_prop_nongerm), 3)
predict$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict


# foxtail
LG<-get.predict(df, "L", "G", "foxtail_prop_nongerm")
LB<-get.predict(df, "L", "B", "foxtail_prop_nongerm")
GB<-get.predict(df, "G", "B", "foxtail_prop_nongerm")
LGB<-get.predict(df, "L", "G", "foxtail_prop_nongerm", "B")
LG

predict1<-rbind(GB, LB, LG, LGB)
predict1$foxtail_prop_nongerm<-round(as.numeric(predict1$foxtail_prop_nongerm), 3)
predict1$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict1
predict<-full_join(predict, predict1)

## add to df
#df1<-df %>% filter(n_species!="1")

df1<-full_join(df, predict) 
df1<-df1 %>% ungroup()
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B", "GB", "GB.predict", "LB", "LB.predict",  "LG", 
                                                "LG.predict", "LGB", "LGB.predict" ))

#### plot predictions from monocultures for weed seed

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



# make labels 
label <- df1$Treatment
label
label <- gsub("LGB.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label

# plot
p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none")+
  labs(title = "C",
       x="",
       y="non germinating foxtail (%)")+
  geom_text(y=.25, label =label , nudge_x = -.8, size=8)


p1

# make labels 
label <- df1$Treatment
label <- gsub("LB.predict", "*" ,label )
label <- gsub("GB.predict", "*" ,label )
label<- gsub("L", "" ,label )
label<- gsub("G", "" ,label )
label <- gsub("B", "" ,label )
label <- gsub(".predict", "" ,label )
label


# plot
p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=pigweed_prop_nongerm, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  geom_text(y=.8, label =label , nudge_x = -.8, size=8)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",)+
  labs(title = "D",
       x="",
       y="non germinating pigweed (%)")

p2

# put the two plots together
require(gridExtra)
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_plant_physio")
#svg("weeddecay.predict.svg", width = 4.5, height = 6)
grid.arrange(p1, p2, ncol=1)
#dev.off()

##########stats weed prop non germinated #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
# group by
df<- df %>% group_by(Rep, Nitrogen)
df

# pigweed non germ 
LG<-get.predict(df, "L", "G", "pigweed_num_nongerm")
LB<-get.predict(df, "L", "B", "pigweed_num_nongerm")
GB<-get.predict(df, "G", "B", "pigweed_num_nongerm")
LGB<-get.predict(df, "L", "G", "pigweed_num_nongerm", "B")
predict<-rbind(GB, LB, LG, LGB)
predict$pigweed_num_nongerm<-round(as.numeric(predict$pigweed_num_nongerm), 0)
predict$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict

# foxtail
LG<-get.predict(df, "L", "G", "foxtail_num_nongerm")
LB<-get.predict(df, "L", "B", "foxtail_num_nongerm")
GB<-get.predict(df, "G", "B", "foxtail_num_nongerm")
LGB<-get.predict(df, "L", "G", "foxtail_num_nongerm", "B")
predict1<-rbind(GB, LB, LG, LGB)
predict1$foxtail_num_nongerm<-round(as.numeric(predict1$foxtail_num_nongerm), 0)
predict1$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict1
predict<-full_join(predict, predict1)


# pigweed 2
LG<-get.predict(df, "L", "G", "pigweed_num_germ")
LB<-get.predict(df, "L", "B", "pigweed_num_germ")
GB<-get.predict(df, "G", "B", "pigweed_num_germ")
LGB<-get.predict(df, "L", "G", "pigweed_num_germ", "B")
predict2<-rbind(GB, LB, LG, LGB)
head(predict2)
predict2$pigweed_num_germ<-round(as.numeric(predict2$pigweed_num_germ), 0)
predict2$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict<-full_join(predict, predict2)
predict


# foxtail
LG<-get.predict(df, "L", "G", "foxtail_num_germ")
LB<-get.predict(df, "L", "B", "foxtail_num_germ")
GB<-get.predict(df, "G", "B", "foxtail_num_germ")
LGB<-get.predict(df, "L", "G", "foxtail_num_germ", "B")
predict3<-rbind(GB, LB, LG, LGB)
predict3$foxtail_num_germ<-round(as.numeric(predict3$foxtail_num_germ), 0)
predict3$Treatment<-c(rep("GB.predict", 6), rep("LB.predict", 6), rep("LG.predict", 6), rep("LGB.predict", 6) )
predict3
predict<-full_join(predict, predict3)
head(predict)

## add to df
df1<-df %>% filter(n_species!="1")

df1<-full_join(df1, predict) 
df1<-df1 %>% ungroup()
as.factor(df1$Treatment)
head(df1)

# anova binomial model pigweed
#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
y<-cbind(df2$pigweed_num_nongerm, df2$pigweed_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

# anova binomial model foxtail
#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.predict")
y<-cbind(df2$foxtail_num_nongerm, df2$foxtail_num_germ)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)


# there are difference between the predicted verse measured for weed seed decay for some treatments

