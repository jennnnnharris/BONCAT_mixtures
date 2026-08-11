#analysis_figs_16s_totaldna_diversity

# 16S analysis
# total DNA
# diversity
# BONCAT mixtures
# Jennifer Harris
# Jan 7 2025
# Last Updated: March 2026

#R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

########Initial Setup ##################

### Clear workspace ###
rm(list=ls())

# install phyloseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")


# Load required libraries #
library(tidyverse)
library(vegan)
library(readxl)
library(lubridate)
library(phyloseq)
library(multcompView)
library(BiodiversityR)

# colors

mycols <- c(
  "#2107EA", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)
mycols <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)


bw <- c( #
  "grey", # grey
  "black" # teal
)


#import data#
# Set the working directory 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")

taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat2.csv", header = T, row.names = 1)
#note metadat2 has block info

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

# order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
row.names(metadat) <- metadat$SampleID
metadat$N   <- factor(metadat$N)
metadat$Legume   <- factor(metadat$Legume)
metadat$Brassicae   <- factor(metadat$Brassicae)
metadat$Grass   <- factor(metadat$Grass)
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes absent", "Legumes present"))
head(metadat)


head(asvs)
#T_DNA_23_S153 has really few reads so I am omitting it.
asvs<-asvs[which(row.names(asvs)!= "T_DNA_23_S153"),]

# get avg number of reads in seq run 
mean(rowSums(asvs))


# select only total dna
asvs <-asvs[which(metadat$Fraction=="Total" ),] 
metadat <- metadat %>% filter(Fraction=="Total" )


#get min number of reads in a sample
min.s<-min(rowSums(asvs))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs<-rrarefy(asvs, min.s)

# make composition var
metadat$composition <- metadat$Treatment
metadat$composition<-gsub("L", "Legume", metadat$composition)
metadat$composition<-gsub("G", "Grass", metadat$composition)
metadat$composition<-gsub("B", "Brassica", metadat$composition)
metadat$composition<-gsub("LegumeBrassica", "Legume_Brassica", metadat$composition)
metadat$composition<-gsub("LegumeGrass", "Legume_Grass", metadat$composition)
metadat$composition<-gsub("GrassBrassica", "Grass_Brassica", metadat$composition)
metadat$composition 
metadat$composition<- factor(metadat$composition, levels = c("Soil", "Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))


#make taxon matrix row names OTUs
row.names(taxon) <- taxon$asv

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps
# 329K taxa when rarefied 


#####remove plant contamination  ###
# Select unassigned Asvs that the only in the the roots and nodules
ps<-subset_taxa(ps, Class!="c__Chloroplast" )
ps<-subset_taxa(ps, Class!=" c__Chloroplast" )
ps<-subset_taxa(ps, Family!= " f__Mitochondria" )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 314K taxa and 84 samples when mitochondria removed. 

#total
ps<-subset_samples(ps, Fraction=="Total")
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 220 k asvs



####### DIVERSITY ####

#overall 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1"))

# Data wrangling for diversity of active microbes in each fraction
metadat.t <- sample_data(ps)
rich<-cbind(rich, sample_data(ps))
rich<-as.data.frame(rich)

# summary
rich %>% group_by(Fraction) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))

#
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))

# caculate pilou's evenness where formula J= H/ln(S), 
rich$evenness = rich$Shannon/log(rich$Observed)



##plots for nitrogen #
p1<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=Shannon))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  theme_classic(base_size = 12)+
  scale_shape_manual(values = c(17, 16)) +
  labs(title = "A",
       x="",
       y= "Shannon diversity")+
  theme(legend.position = "none")

p1
p2<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=Observed))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  theme_classic(base_size = 12)+
  scale_shape_manual(values = c(17, 16)) +
  labs(title = "B",
       x="",
       y= "ASV richness") +
  theme(legend.position = "none")

p2
p3<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=evenness))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.2)+
  theme_classic(base_size = 12)+
  scale_shape_manual(values = c(17, 16)) +
  labs(title = "C",
       x="",
       y= "Pilou's Eveness")#+
# theme(legend.position = "none")

p3
require(gridExtra)
#windows(8, 3.5)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement")
svg("diversity.nitrogenkey.svg", width=7, height = 3)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()


## plots for treatment ##
# shannon diversity
rich<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil")



p1<- rich%>%  
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.1)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "D",
       x="",
       y= "Shannon Diversity")+
  scale_shape_manual(values = c(17, 16))+
  facet_grid(~n_species, space="free", scales="free" )
p1



# OBSERVED ASVS

p2<- rich%>%  
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.1)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "E",
       x="",
       y= "Observed ASVs")+
  scale_shape_manual(values = c(17, 16))+
  facet_grid(~n_species, space="free", scales="free" )
p2




# EVENESSS

p3<- rich%>%  
  ggplot(aes(x=Treatment, y=evenness,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1, width=.1)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "F",
       x="",
       y= "Eveness")+
  scale_shape_manual(values = c(17, 16))+
  facet_grid(~n_species, space="free", scales="free" )
p3





setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Supplement")
svg(file="diversity.treatment.svg",width = 8.5, height=3.5)
require(gridExtra)
#windows(10,4)
grid.arrange(p1, p2,p3, ncol=3)
dev.off()



####STATS#####
library(emmeans)
library(multcompView)
# remove soil
rich <- rich %>% filter(Fraction=="Total") %>% filter(Treatment!="Soil")
rich$Treatment   <- factor(rich$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
rich$N <- factor(rich$N, levels= c("0", "1"))


#shannon 
# overall model 
m1<-lm(Shannon ~ N*Treatment,  data = rich)
summary(m1)
# ANOVA
anova1<- aov(Shannon ~ Treatment*N, data = rich)
summary(anova1) # block droped because not sig

## ANOVA andtukey post hoc on simple model
m1<-lm(Shannon ~ Treatment,  data = rich)
anova1<- aov(Shannon ~ Treatment, data = rich)
summary(anova1)

# tukey pos hoc
tukey <- TukeyHSD(anova1)
print(tukey) 

# emmeans post hoc# Perform all pairwise comparisons with Tukey adjustment
emm_object <- emmeans(m1, specs = ~ Treatment)
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

# get letters for plotting
cld <- multcompLetters4(anova1, tukey)
print(cld)



##observed
# overall model 
m1<-lm(Observed ~ N*Treatment,  data = rich)
summary(m1)
# ANOVA
anova1<- aov(Observed ~ Treatment*N, data = rich)
summary(anova1) # block droped because not sig

## ANOVA andtukey post hoc on simple model
m1<-lm(Observed ~ Treatment,  data = rich)
anova1<- aov(Observed ~ Treatment, data = rich)
summary(anova1)

# tukey pos hoc
tukey <- TukeyHSD(anova1)
print(tukey) 

# emmeans post hoc# Perform all pairwise comparisons with Tukey adjustment
emm_object <- emmeans(m1, specs = ~ Treatment)
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

# get letters for plotting
cld <- multcompLetters4(anova1, tukey)
print(cld)


# Evenness
# overall model 
m1<-lm(evenness ~ N*Treatment,  data = rich)
summary(m1)
# ANOVA
anova1<- aov(evenness ~ Treatment*N, data = rich)
summary(anova1) # block droped because not sig

## ANOVA andtukey post hoc on simple model
m1<-lm(evenness ~ Treatment,  data = rich)
anova1<- aov(evenness ~ Treatment, data = rich)
summary(anova1)

# tukey pos hoc
tukey <- TukeyHSD(anova1)
print(tukey) 

# emmeans post hoc# Perform all pairwise comparisons with Tukey adjustment
emm_object <- emmeans(m1, specs = ~ Treatment)
pairwise_comparison <- pairs(emm_object, adjust = "tukey") 
summary(pairwise_comparison)

# get letters for plotting
cld <- multcompLetters4(anova1, tukey)
print(cld)

