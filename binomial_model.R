## my data
rm(list=ls())
library(readxl)
library(tidyverse)
library(lubridate)
library(boot)




set.seed(1)
n <- 50
cov <- 10 #"size
x <- c(rep(0,n/2), rep(1, n/2))
p <- 0.4 + 0.2*x
y <- rbinom(n, cov, p)
y

?rbinom()

model0 <- glm(cbind(y, cov-y) ~ x, family="binomial")
summary(model0)

model_intercept <- 1/(1+(1/(exp(model0[[1]][1]))))
model_intercept # probabilty of success for x = 0
boot::inv.logit((model0$coefficients[1]))

model_intercept # probabilty of success for x = 0

# (Intercept) 
#       0.424
model_x <- 1/(1+(1/(exp(model0[[1]][2]))))

model_x
#         x 
# 0.6634292  # probalbity of succuss for x = 1


# amount of increase for X =1 compared to x =0

inc<-model_x- model_intercept # proablity increase

cov * inc # increase in events



## my data ############
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
# think i should probably use a link logit model to this in the future, however back transforming the effect size can be tricky.

df<-df%>%
  mutate(n_failures =  viablecells_count-BONCAT_count)

y<-cbind(df$BONCAT_count, df$n_failures)

m1<-glm(data= df, y~ Treatment + Date_sorted + Block, family = binomial)
summary(m1)

m1$coefficients[1]
m1$coefficients[9:20]

intercept<-boot::inv.logit(m1$coefficients[1])
intercept
m1$coefficients[9:20]
n<-boot::inv.logit(m1$coefficients[9:20])
n
del<-n - .5
del
1000 * .03
# mutliplintercept# mutliple del times the number of total events. 
# that would give the increase. 

df$viablecells_count * .2
