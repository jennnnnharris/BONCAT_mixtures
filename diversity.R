# diversity
# 16S analysis
# BONCAT mixtures
# Jennifer Harris
# Jan 7 2025
# Last Updated: May 19 2025

# R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

# r version 
### 1. Initial Setup ###

### Clear workspace ###

rm(list=ls())


################## Load required libraries ############

#basic
library(tidyverse)
library(vegan)
library(readxl)
library(lubridate)




#Phyloseq and mbiome #if trouble loading phyloseq, see: http://joey711.github.io/phyloseq/install
library(phyloseq)

#install.packages("Rtools")

#if(!"devtools" %in% installed.packages()){
#  install.packages("devtools")
#}
#devtools::install_github("gmteunisse/fantaxtic")

#install.packages("fantaxtic")
#library(fantaxtic)
#library(microbiome)
#library(DT)

# venn diagrams
#library(ggvenn)
#library(grid)

#Ancom
#library(ANCOMBC)
#load("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_mixtures/16S.RData")


#colors
mycols7<-c( "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF", "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
mycols8<-c("grey", "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF",  "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
mycols3<- c(  "#f4f1bb", "#ed6a5a","#9bc1bc")
mycols3<- c( "#006d77",  "#f4d35e", "#e94f37")
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")

####
#("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_mixtures/16S.RData")
#####Import data#####
## Set the working directory ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")

taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read_excel("metadata.xlsx", sheet = 1)

## Transpose ASVS table ##
asvs <- t(asvs)
asvs[1:5,1:5]#taxa are columns

## check the reads per sample ##
#rowSums(asvs)[order(rowSums(asvs))]

#T_DNA_23_S153 has really few reads so I am omitting it.
asvs<-asvs[which(row.names(asvs)!= "T_DNA_23_S153"),]


#get min number of reads in a sample
min.s<-min(rowSums(asvs))


## order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
row.names(metadat) <- metadat$SampleID
metadat

# length(intersect(colnames(asvs.raw) , metadat$SampleID)) # Apply setdiff function to see what's missing from the tree
# mynames<- setdiff(metadat$SampleID , colnames(asvs.raw) )
# length(mynames)
# mynames

#make taxon matrix row names OTUs
row.names(taxon) <- taxon$Feature.ID

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps
# 329K taxa 

####### rarefaction curve ######

#use not rarefied asvs. 
asvs[1:5,1:5]
S <- specnumber(asvs) # observed number of species
S
(raremax <- min(rowSums(asvs)))
Srare <- rarefy(asvs, raremax)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")

svg("observed.vs.rarefied.svg", )
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
dev.off()

svg("rarefaction.svg", height=8, width=8 )
rarecurve(asvs, step = 10000, sample = raremax, col = "blue", cex = 0.6)
dev.off()

i<-asvs[which(metadat$Fraction=="Active"),]
svg("r.plot2.svg", )
rarecurve(i, step = 10000, sample = raremax, col = "blue", cex = 0.6)
dev.off()

i<-asvs[which(metadat$Fraction=="Total"),]
i[1:5,1:15]
svg("r.plot3.svg", )
rarecurve(i, step = 10000, sample = raremax, col = "blue", cex = 0.6)
dev.off()

##can you rarefy by number of cells. 


#####remove plant contamination  ########
# Select unassigned Asvs that the only in the the roots and nodules
ps<-subset_taxa(ps, Class!="c__Chloroplast" )
ps<-subset_taxa(ps, Class!=" c__Chloroplast" )
ps<-subset_taxa(ps, Family!= " f__Mitochondria" )

ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 314K taxa and 84 samples when mitochondria removed. 



#######filtering total##################
#total
total<-subset_samples(ps, Fraction=="Total")
total<-prune_taxa(taxa_sums(total) > 0, total)
total

# remove true singletons 
total<-prune_taxa(taxa_sums(total) > 1, total)

#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(total)))/nsamples(total)
keep<-row.names(t(otu_table(total))[ mean.reads > 5, ])
total<-prune_taxa(keep, total)

# remove asvs that are in less than 5 samples
total <-ps_prune(total, min.samples = 3)
total
#1696 taxa

## plot
#plot(sort(taxa_sums(total), TRUE), type="h", ylim=c(0, 8000))

#######filtering active + inactive ##################
fc<-subset_samples(ps, Fraction!="Total" & Fraction!="CTL")
fc<-prune_taxa(taxa_sums(fc) > 0, fc)
fc

# remove true singletons # 183000 ASVS 
fc<-prune_taxa(taxa_sums(fc) > 1, fc)
fc

#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(fc)))/nsamples(fc)
keep<-row.names(t(otu_table(fc))[ mean.reads > 5, ])
fc<-prune_taxa(keep, fc)
fc

# remove asvs that are in less than 5 samples
fc <-ps_prune(fc, min.samples = 3)
fc

#1864 taxa
plot(sort(taxa_sums(fc), TRUE), type="h", ylim=c(0, 8000))



####DIVERSITY plots  ####

#overall 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1"))

# Data wrangling fo rdiversity of active microbes in each fraction
metadat.t <- sample_data(ps)
rich<-cbind(rich, sample_data(ps))
rich<-as.data.frame(rich)

# summary
rich %>% group_by(Fraction) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))

#colors
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
rich$Fraction   <- factor(rich$Fraction, levels= c("Total", "Active", "Inactive"))


# Does diversity increase with number of species?
p1 <- rich%>% filter(Treatment!="NA") %>%
  ggplot(aes(x=Fraction, y=Shannon,  fill=Fraction))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols3) +
  scale_fill_manual(values = mycols3)+
  geom_jitter(width = .1, size=1 )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Shannon Diversity ")+
  xlab("")

p1

p2<- rich%>% filter(Treatment!="NA") %>%
  ggplot(aes(x=Fraction, y=Observed,  fill=Fraction))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols3) +
  scale_fill_manual(values = mycols3)+
  geom_jitter(width = .1, size=1 )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Observed Diversity ")+
  xlab("")

p2

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="diversity.fraction.svg",width = 5, height=4)
p2
dev.off()


svg(file="activity.species.svg",width = 3, height=3)
#require(gridExtra)
#windows(8,4)
df  %>%  filter(Species1!="Soil") %>%
  ggplot(aes(x=n_species, y=BONCAT_freq )) +
  geom_jitter(width = .2, size=1 )+
  geom_smooth(method = lm, color= blues[4])+
  theme_classic(base_size = 14)+
  theme( legend.position="none",
         plot.title = element_text(hjust = 0.5))+
  ylab("percent active")+
  xlab("n species")

dev.off()


#total
rich<-estimate_richness(total, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))


# Data wrangling fo rdiversity of active microbes in each fraction
metadat.t <- sample_data(total)
rich<-cbind(rich, sample_data(total))
rich<-as.data.frame(rich)
rich %>%
  arrange( -Observed)

# summary
rich %>% group_by(Fraction) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))

#colors
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
#mycols8<- c( "grey", "#1F78B4",  "#eb05db", "#33A02C", "#FF7F00","#1a635a", "#6A3D9A", "yellow")


# Does diversity increase with number of species?
p1<-rich%>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  geom_jitter(width = .1, size=1 )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Shannon Diversity ")+
  xlab("")+
  labs(title = "Total")+
  facet_grid( ~n_species, scales = "free", space = "free")
p1

p2<-rich%>%
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  geom_jitter(width = .1, size=1 )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("N Asvs")+
  xlab("")+
  facet_grid( ~n_species, scales = "free", space = "free")

p2

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="diversity.total.svg",width = 6, height=7)

require(gridExtra)
#windows(6,8)
grid.arrange(p1, p2, ncol=1)
dev.off()

######active #####
rich<-estimate_richness(fc, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

# Data wrangling for diversity of active microbes in each fraction
metadat.t <- sample_data(fc)
rich<-as.data.frame(cbind(rich, sample_data(fc)))
rich %>%  arrange( -Observed)

# summary
rich %>% group_by(Fraction) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))

#colors
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))

rich$Trt_fraction <- factor(rich$Trt_fraction, levels =c("Soil_Active",   "Soil_Inactive" ,"B_Active",      "B_Inactive",    "G_Active" ,     "G_Inactive",    "GB_Active",     "GB_Inactive",  
                                                         "L_Active",      "L_Inactive",    "LB_Active" ,    "LB_Inactive" ,  "LG_Active",     "LG_Inactive",  
                                                         "LGB_Active"  ,  "LGB_Inactive"  ))


# Does diversity increase with number of species?

p1<-  rich%>%  filter(Fraction=="Active") %>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Shannon Diversity")+
  xlab("")+
  labs(title = "Active")+
  facet_grid( ~n_species, scales = "free", space = "free")+
  scale_shape_discrete() 

p1

p2<-rich%>%  filter(Fraction=="Active") %>%
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("N Asvs")+
  xlab("")+
  facet_grid( ~n_species, scales = "free", space = "free")+
  scale_shape_discrete() 

p2
#inactive
# Does diversity increase with number of species?
p3<-rich%>% filter(Fraction=="Inactive") %>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Shannon Diversity ")+
  xlab("")+
  labs(title = "Inactive")+
  facet_grid( ~n_species, scales = "free", space = "free")+
  scale_shape_discrete() 
p3

p4<-rich%>%  filter(Fraction=="Inactive") %>%
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  scale_color_manual(values=mycols8) +
  scale_fill_manual(values = mycols8)+
  geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  theme_minimal(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("N Asvs")+
  xlab("")+
  facet_grid( ~n_species, scales = "free", space = "free")+
  scale_shape_discrete() 


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="diversity.fc.svg",width = 10, height=6)

require(gridExtra)
#windows(6,8)
grid.arrange(p1, p3, p2, p4, ncol=2)
dev.off()





##DIVERSITY STATS######
####shannon###
# overall anova

#total
rich<-estimate_richness(total, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

# Data wrangling fo rdiversity of active microbes in each fraction
metadat.t <- sample_data(total)
rich<-cbind(rich, sample_data(total))
rich<-as.data.frame(rich)

# summary
rich %>% group_by(Treatment) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))

#Shannon
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
m1<-lm(Shannon ~ Treatment*N,  data = rich)
summary(m1)
# different than soil 
#remove soil
rich<-rich %>% filter(Treatment!="Soil")
rich$Treatment   <- factor(rich$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# run model
m1<-lm(Shannon ~ Treatment*N,  data = rich)
summary(m1)
plot(m1)
# LG is higher 
# marginal effect of nitrogen

#observed
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
m1<-lm(Observed ~ Treatment*N,  data = rich)
summary(m1)
# higher than soil
#remove soil
rich<-rich %>% filter(Treatment!="Soil")
rich$Treatment   <- factor(rich$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# run model
m1<-lm(Observed ~ Treatment*N,  data = rich)
summary(m1)
plot(m1)
# LG is higher 
# marginal effect of nitrogen

# inactive
rich<-estimate_richness(fc, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

# Data wrangling fo rdiversity of active microbes in each fraction
metadat.t <- sample_data(fc)
rich<-cbind(rich, sample_data(fc))
rich<-as.data.frame(rich)

# summary
rich %>% group_by(Treatment) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(Fraction) %>% summarise(mean(Observed), sd(Observed))


#Shannon
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
m1<-lm(Shannon ~ Treatment*Fraction,  data = rich)
summary(m1)
# different than soil 
#remove soil
rich<-rich %>% filter(Treatment!="Soil" & Fraction=="Inactive")
rich$Treatment   <- factor(rich$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# run model
m1<-lm(Shannon ~ Treatment,  data = rich)
summary(m1)
plot(m1)


#observed
#remove soil
rich<-rich %>% filter(Treatment!="Soil" & Fraction=="Inactive")
rich$Treatment   <- factor(rich$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
# run model
m1<-lm(Observed ~ Treatment,  data = rich)
summary(m1)
plot(m1)

