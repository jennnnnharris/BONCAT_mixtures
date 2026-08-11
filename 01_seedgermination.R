# weed seed germination

#load libraries
library(tidyverse)


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
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/seedgermination")
df <- read.csv("weedseedgermination.csv", row.names = 1 )
head(df)

# set factor
df$Treatment   <- factor(df$Treatment, levels= c("S", "L", "G", "B", "GB", "LB", "LG", "LGB"))
unique(df$Treatment)


df<-df  %>%   filter(Treatment!="S") 
df$Treatment   <- factor(df$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

df<-df %>%
  mutate(df, foxtail_prop_germ = foxtail_num_germ/foxtail_total_seeds,
         pigweed_prop_germ= pigweed_num_germ/pigweed_total_seeds)


p1<- df %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_germ, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "none")+
  scale_fill_manual(values = mycols)+
  labs(title = "A",
       x="",
       y="germinating foxtail (%)")
p1



p2<-df %>%
  ggplot(aes(x=Treatment, y=pigweed_prop_germ, fill=Treatment)) + 
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(size=.5)+
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position = "none")+
  scale_fill_manual(values = mycols)+
  labs(title = "B",
       x="",
       y="germinating pigweed (%)")
p2


require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement")
svg("weedgermination.svg", width = 5, height = 2.5)
grid.arrange(p1, p2, ncol=2)
dev.off()



### stats #########
# model foxtail
y<-cbind(df$foxtail_num_germ, df$foxtail_num_nongerm)
m1<-glm(data= df, y~Treatment+Block, family = binomial)
summary(m1)

# Run Anova on your binom model
library(car)
Anova(m1, type = "II")

# significant effect of block and species and species treatment interaction
# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
emm_object <- emmeans(m1, specs = ~ Treatment)
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
# model foxtail
y<-cbind(df$pigweed_num_germ, df$pigweed_num_nongerm)
m1<-glm(data= df, y~Treatment+Block, family = binomial)
summary(m1)


library(car)
Anova(m1, type = "II")
# significant effect of block and species and species treatment interaction
# pairwise tests
library(emmeans)
# Get the EMMs for your treatment groups
emm_object <- emmeans(m1, specs = ~ Treatment)
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
factor(df$Treatment)

########expected ##########
get.expect<-function(df, sp1, sp2, trait, sp3) {
  if(missing(sp3)){
    print(trait)
    sp1.df<- filter(df, Treatment==sp1) %>% select(all_of(trait))
    print(sp1.df)
    sp1.expect<-sp1.df/2
    print(sp1.expect)
    sp2.df<-filter(df, Treatment==sp2) %>% select(all_of(trait))
    sp2.expect<-sp2.df /2
    print(sp2.expect)
    expect <- sp1.expect + sp2.expect
    print(expect)
    #return(expect)
  }
  else {
    print(paste("the 3rd species is",sp3))
    print(trait)
    sp1.df<- filter(df, Treatment==sp1) %>% select(all_of(trait))
    sp1.expect<-sp1.df/3
    #print(sp1.df)
    #print(sp1.expect)
    sp2.df<-filter(df, Treatment==sp2) %>% select(all_of(trait))
    sp2.expect<-sp2.df/3
    print(sp2.df)
    print(sp2.expect)
    sp3.df<-filter(df, Treatment==sp3) %>% select(all_of(trait))
    sp3.expect<-sp3.df/3
    #print(sp3.df)
    #print(sp3.expect)
    expect <- sp1.expect + sp2.expect + sp3.expect
    print(expect)
    return(expect)
  }
  
}

######expections weed seed########
library(tidyverse)

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
df<-df %>%
  mutate(df, foxtail_prop_germ = foxtail_num_germ/foxtail_total_seeds,
         pigweed_prop_germ= pigweed_num_germ/pigweed_total_seeds)
# group by
df<- df %>% group_by(Rep, Nitrogen)


# pigdf
LG<-get.expect(df, "L", "G", "pigweed_prop_germ")
LB<-get.expect(df, "L", "B", "pigweed_prop_germ")
GB<-get.expect(df, "G", "B", "pigweed_prop_germ")
LGB<-get.expect(df, "L", "G", "pigweed_prop_germ", "B")
LG

expect<-rbind(GB, LB, LG, LGB)
expect$pigweed_prop_germ<-round(as.numeric(expect$pigweed_prop_germ), 3)
expect$Treatment<-c(rep("GB.expect", 6), rep("LB.expect", 6), rep("LG.expect", 6), rep("LGB.expect", 6) )
expect


# foxtail
LG<-get.expect(df, "L", "G", "foxtail_prop_germ")
LB<-get.expect(df, "L", "B", "foxtail_prop_germ")
GB<-get.expect(df, "G", "B", "foxtail_prop_germ")
LGB<-get.expect(df, "L", "G", "foxtail_prop_germ", "B")
LG

expect1<-rbind(GB, LB, LG, LGB)
expect1$foxtail_prop_germ<-round(as.numeric(expect1$foxtail_prop_germ), 3)
expect1$Treatment<-c(rep("GB.expect", 6), rep("LB.expect", 6), rep("LG.expect", 6), rep("LGB.expect", 6) )
expect1
expect<-full_join(expect, expect1)

## add to df
#df1<-df %>% filter(n_species!="1")

df1<-full_join(df, expect) 
df1<-df1 %>% ungroup()
df1$Treatment<-factor(df1$Treatment, levels = c("L", "G", "B",  "GB.expect", "GB","LB.expect","LB", "LG.expect",   "LG", 
                                                "LGB.expect","LGB" ))
head(df1)
df1
  #### plot expections from monocultures for weed seed

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
  "#865338" # medium mocha brown LGB
 
)




# plot
p1<-df1  %>% 
  ggplot(aes(x=Treatment, y=foxtail_prop_germ, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position="none")+
  labs(title = "",
       x="",
       y="germinating foxtail (%)")
 # geom_text(y=1, label =label , nudge_x = +.6, size=8)


p1




# plot
p2<-df1  %>% 
  ggplot(aes(x=Treatment, y=pigweed_prop_germ, fill = Treatment)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  theme_classic(base_size = 12)+
  #geom_text(y=1, label =label , nudge_x = +.6, size=8)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position="none",)+
  labs(title = "",
       x="",
       y="germinating pigweed (%)")+
  ylim(c(.2,1))

p2

# put the two plots together
require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement/weed_germination")
svg("weed.expect.svg", width = 4.5, height = 6)
grid.arrange(p1, p2, ncol=1)
dev.off()

##########stats weed prop non germinated #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
df<-read.csv("weed_seed_decay.csv", row.names = 1)
head(df)
df<-df %>% filter(Treatment!="S")
# group by
df<- df %>% group_by(Rep, Nitrogen)
df

# pigweed non germ 
LG<-get.expect(df, "L", "G", "pigweed_num_nongerm")
LB<-get.expect(df, "L", "B", "pigweed_num_nongerm")
GB<-get.expect(df, "G", "B", "pigweed_num_nongerm")
LGB<-get.expect(df, "L", "G", "pigweed_num_nongerm", "B")
expect<-rbind(GB, LB, LG, LGB)
expect$pigweed_num_nongerm<-round(as.numeric(expect$pigweed_num_nongerm), 0)
expect$Treatment<-c(rep("GB.expect", 6), rep("LB.expect", 6), rep("LG.expect", 6), rep("LGB.expect", 6) )
expect

# foxtail
LG<-get.expect(df, "L", "G", "foxtail_num_nongerm")
LB<-get.expect(df, "L", "B", "foxtail_num_nongerm")
GB<-get.expect(df, "G", "B", "foxtail_num_nongerm")
LGB<-get.expect(df, "L", "G", "foxtail_num_nongerm", "B")
expect1<-rbind(GB, LB, LG, LGB)
expect1$foxtail_num_nongerm<-round(as.numeric(expect1$foxtail_num_nongerm), 0)
expect1$Treatment<-c(rep("GB.expect", 6), rep("LB.expect", 6), rep("LG.expect", 6), rep("LGB.expect", 6) )
expect1
expect<-full_join(expect, expect1)


# pigweed 2
LG<-get.expect(df, "L", "G", "pigweed_num_germ")
LB<-get.expect(df, "L", "B", "pigweed_num_germ")
GB<-get.expect(df, "G", "B", "pigweed_num_germ")
LGB<-get.expect(df, "L", "G", "pigweed_num_germ", "B")
expect2<-rbind(GB, LB, LG, LGB)
head(expect2)
expect2$pigweed_num_germ<-round(as.numeric(expect2$pigweed_num_germ), 0)
expect2$Treatment<-c(rep("GB.expect", 6), rep("LB.expect", 6), rep("LG.expect", 6), rep("LGB.expect", 6) )
expect<-full_join(expect, expect2)
expect


# foxtail
LG<-get.expect(df, "L", "G", "foxtail_num_germ")
LB<-get.expect(df, "L", "B", "foxtail_num_germ")
GB<-get.expect(df, "G", "B", "foxtail_num_germ")
LGB<-get.expect(df, "L", "G", "foxtail_num_germ", "B")
expect3<-rbind(GB, LB, LG, LGB)
expect3$foxtail_num_germ<-round(as.numeric(expect3$foxtail_num_germ), 0)
expect3$Treatment<-c(rep("GB.expect", 6), rep("LB.expect", 6), rep("LG.expect", 6), rep("LGB.expect", 6) )
expect3
expect<-full_join(expect, expect3)
head(expect)

## add to df
df1<-df %>% filter(n_species!="1")

df1<-full_join(df1, expect) 
df1<-df1 %>% ungroup()
as.factor(df1$Treatment)
head(df1)

# anova binomial model pigweed
#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.expect")
y<-cbind(df2$pigweed_num_germ, df2$pigweed_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.expect")
y<-cbind(df2$pigweed_num_germ, df2$pigweed_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.expect")
y<-cbind(df2$pigweed_num_germ, df2$pigweed_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.expect")
y<-cbind(df2$pigweed_num_germ, df2$pigweed_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

# anova binomial model foxtail
#GB
df2<-df1%>% filter(Treatment=="GB" | Treatment=="GB.expect")
y<-cbind(df2$foxtail_num_germ, df2$foxtail_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LB
df2<-df1%>% filter(Treatment=="LB" | Treatment=="LB.expect")
y<-cbind(df2$foxtail_num_germ, df2$foxtail_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LG
df2<-df1%>% filter(Treatment=="LG" | Treatment=="LG.expect")
y<-cbind(df2$foxtail_num_germ, df2$foxtail_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)

#LGB
df2<-df1%>% filter(Treatment=="LGB" | Treatment=="LGB.expect")
y<-cbind(df2$foxtail_num_germ, df2$foxtail_num_nongerm)
m1<-glm(data= df2, y~Treatment, family = binomial)
summary(m1)


# there are difference between the expected verse measured for weed seed decay for some treatments

