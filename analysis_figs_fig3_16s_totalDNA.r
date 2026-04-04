# 16S analysis
# total DNA
# beta diversity
# BONCAT mixtures
# Jennifer Harris
# Jan 7 2025
# Last Updated: May 2025

#R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

######## 1. Initial Setup ##################

### Clear workspace ###
#rstudioapi::restartSession(clean = TRUE)
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
IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)

"#8cadff"
bw <- c( #IBM colors
  "grey", # grey
  "black" # teal
)

# set shapes
myshapes <- c(1, 12,15 ,21, 22, 23 , 24)
myshapes2 <- c(21 , 12, 24,1, 15 , 22, 23 )



  
#import data#
# Set the working directory 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")

taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat2.csv", header = T, row.names = 1)
#note metadat2 had block info

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



####### 2. DIVERSITY ####

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



##plots##

# nitrogen effect
p1<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=Shannon))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  labs(title = "A",
       x="",
       y= "Total DNA Shannon diversity")
  #annotate("text", x=1.5, y=8, label="*")
p1
p2<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=Observed))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  labs(title = "B",
       x="",
       y= "Total DNA ASV richness")

p2
p3<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=evenness))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  labs(title = "C",
       x="",
       y= "Total DNA Pilou's Eveness")
  #annotate("text", x=1.5, y=5500, label="*")
p3
require(gridExtra)
windows(8, 3.5)
grid.arrange(p1, p2, p3, ncol=3)


# shannon diversity
rich<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil")
N1<-rich %>% filter(N==1)
N0 <- rich %>% filter(N==0)

# nitrogen +  
lab<-as.character(N1$Treatment)
lab<- gsub("LGB", "A", lab)
lab <- gsub("LG", "A", lab)
lab <- gsub("LB", "A", lab)
lab <- gsub("GB", "A", lab)
lab <- gsub("L", "A", lab) 
lab <- gsub("G", "A", lab) 
lab <- gsub("B", "A", lab) 
lab


N1$Nitrogen_label<-as.factor(N1$Nitrogen_label)
p1<-N1%>%  
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  geom_text(y=7.89, label = lab, size=5)+
  labs(title = "E",
       subtitle = "Nitrogen +",
       x="",
       y= "Total DNA Shannon Diversity")+
  scale_shape_manual(values = c(17, 16))+
  ylim(6.6, 8)
p1
# nitrogen - shannon 
#  L  G    B   GB   LB   LG  LGB    
# "c" "ab"  "a" "ab" "ab" "ab" "bc" 
lab<-as.character(N0$Treatment)
lab<- gsub("LGB", "XC", lab)
lab <- gsub("LG", "AX", lab)
lab <- gsub("LB", "AX", lab)
lab <- gsub("GB", "AX", lab)
lab
lab <- gsub("L", "C", lab) 
lab <- gsub("G", "AX", lab) 
lab <- gsub("X", "B", lab) 
lab


N0$Nitrogen_label<-as.factor(N0$Nitrogen_label)
p2<-N0%>%  
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  geom_text(y=7.89, label = lab, size=5)+
  labs(title = "F",
       subtitle = "Nitrogen -",
       x="",
       y= "Total DNA Shannon Diversity")+
  scale_shape_manual(values = c(17, 16))+
  ylim(6.6, 8)

p2
grid.arrange(p1, p2, ncol=2)


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
svg(file="diversity.svg",width = 8, height=4)
require(gridExtra)
#windows(10,4)
grid.arrange(p1, p2, ncol=2)
dev.off()





# OBSERVED ASVS
rich<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil")
N1<-rich %>% filter(N==1)
N0 <- rich %>% filter(N==0)

N1$Nitrogen_label<-as.factor(N1$Nitrogen_label)
p1<-N1%>%  
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  #geom_text(y=7.89, label = lab, size=5)+
  labs(title = "A",
       subtitle = "Nitrogen +",
       x="",
       y= "Total DNA ASV richness")+
  scale_shape_manual(values = c(17, 16))+
  ylim(2500, 5500)
p1

# nitrogen minus OBSERVED ASV
# L   G    B   GB   LB   LG  LGB   
#"b" "ab"  "a"  "a" "ab"  "a" "ab"  
lab<-as.character(N0$Treatment)
lab<- gsub("LGB", "AX", lab)
lab <- gsub("LG", "A", lab)
lab <- gsub("LB", "AX", lab)
lab <- gsub("GB", "A", lab)
lab <- gsub("L", "X", lab) 
lab <- gsub("G", "AX", lab) 
lab <- gsub("B", "A", lab) 
lab <- gsub("X", "B", lab) 

N0$Nitrogen_label<-as.factor(N0$Nitrogen_label)
p2<-N0%>%  
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  geom_text(y=5300, label = lab, size=4)+
  labs(title = "B",
       subtitle = "Nitrogen -",
       x="",
       y= "Total DNA  ASV Richness")+
  scale_shape_manual(values = c(17, 16))+
  ylim(2500, 5500)
p2

# EVENESSS
rich$Nitrogen_label<-as.factor(rich$Nitrogen_label)
# L  G    B   GB   LB   LG  LGB    
#"b" "a"  "a" "ab" "ab" "ab" "ab"  

lab<-as.character(rich$Treatment)
lab<- gsub("LGB", "AX", lab)
lab <- gsub("LG", "AX", lab)
lab <- gsub("LB", "AX", lab)
lab <- gsub("GB", "AX", lab)
lab <- gsub("L", "X", lab) 
lab <- gsub("G", "A", lab) 
lab <- gsub("B", "A", lab) 
lab <- gsub("X", "B", lab) 
lab

p3<-rich%>%  
  ggplot(aes(x=Treatment, y=evenness,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_text(y=.945, label = lab, size=4)+
  labs(title = "C",
       subtitle = "Evenness",
       x="",
       y= " Pilou's Evenness")+
  scale_shape_manual(values = c(17, 16))+
  ylim(.85, .95)

  p3

windows(8,3.5)
grid.arrange(p1, p2, p3, ncol=3)




#STATS#
# remove soil
rich <- rich %>% filter(Fraction=="Total") %>% filter(Treatment!="Soil")
rich$Treatment   <- factor(rich$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
rich$N <- factor(rich$N, levels= c("0", "1"))
#rich$Treatment <- as.character(rich$Treatment)

# Shannon
m1<-lm(Shannon ~ N,  data = rich)
summary(m1)

anova1<- aov(Shannon ~ Treatment*N, data = rich)
summary(anova1) # block droped because not sig

# shannon N+ anonva
N1<-rich %>% filter(N==1)
N0 <- rich %>% filter(N==0)
anova1<- aov(Shannon ~ Treatment, data = N1)
summary(anova1)
library(multcompView)
tukey.a1 <- TukeyHSD(anova1)
print(tukey.a1) 
#plot(anova1) #homoscedasticity looks fine
#Extract the p-values for the factor of interest
p_values <- tukey.a1$Treatment[, 4]
# Generate the grouping letters using multcompLetters()
cld <- multcompLetters(p_values)
print(cld)

# shannon N- anova
anova1<- aov(Shannon ~ Treatment, data = N0)
summary(anova1)
library(multcompView)
tukey.a1 <- TukeyHSD(anova1)
print(tukey.a1) 
#plot(anova1) #homoscedasticity looks fine
#Extract the p-values for the factor of interest
p_values <- tukey.a1$Treatment[, 4]
# Generate the grouping letters using multcompLetters()
cld <- multcompLetters(p_values)
print(cld)

##observed
#anova
anova1<- aov(Observed ~ Treatment*N, data = rich)
summary(anova1)

# break up N+ and N-
N1<-rich %>% filter(N==1)
N0 <- rich %>% filter(N==0)
anova1<- aov(Observed ~ Treatment, data = N1)
summary(anova1)
library(multcompView)
tukey.a1 <- TukeyHSD(anova1)
print(tukey.a1) 
#plot(anova1) #homoscedasticity looks fine
#Extract the p-values for the factor of interest
p_values <- tukey.a1$Treatment[, 4]
# Generate the grouping letters using multcompLetters()
cld <- multcompLetters(p_values)
print(cld)

## N-
anova1<- aov(Observed ~ Treatment, data = N0)
summary(anova1)
library(multcompView)
tukey.a1 <- TukeyHSD(anova1)
print(tukey.a1) 
#plot(anova1) #homoscedasticity looks fine
#Extract the p-values for the factor of interest
p_values <- tukey.a1$Treatment[, 4]
# Generate the grouping letters using multcompLetters()
cld <- multcompLetters(p_values)
print(cld)


# Evenness
anova1<- aov(evenness ~ Treatment*N, data = rich)
summary(anova1)
# no interaction with nitrogen, so I'm just gonna bin +N =N together

#tukey
anova1<- aov(evenness ~ Treatment, data = rich)
summary(anova1)
library(multcompView)
tukey.a1 <- TukeyHSD(anova1)
print(tukey.a1) 
#plot(anova1) #homoscedasticity looks fine
#Extract the p-values for the factor of interest
p_values <- tukey.a1$Treatment[, 4]
# Generate the grouping letters using multcompLetters()
cld <- multcompLetters(p_values)
print(cld)



######## 3. Filtering rare taxa ##################

# remove true singletons 
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
ps
# 183095 

#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps


# remove asvs that are in less than 5 samples
#ps <-ps_prune(ps, min.samples = 3) # no features to group!
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
#1778 asvs

## plot
#plot(sort(taxa_sums(ps), TRUE), type="h", ylim=c(0, 8000))


######## 4. CAP -Treatment- ##################
# Constrained ordination
# Perform vegdist analysis of BC distances #

# subset data
ps1 <-subset_samples(ps, Treatment !="Soil" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Treatment*N+block,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

anova.cca(cap_result, by="terms")

### 3. grab info for the plot
smry <- summary(cap_result)
smry
sc_si <- scores(cap_result, display="sites", choices=c(1,2), scaling=1)
sc_si

# Extract the model's adjusted R2
RsquareAdj(cap_result)$adj.r.squared

# percent varience of total varience on RDA 1 and RDA 2
perc <- round(100*(summary(cap_result)$cont$importance[2, 1:2]), 2)
perc


### 4. plot 
#
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal") 
svg("cap.total1.svg", width = 3, height = 7)
#windows(4,7)
par(mfrow=c(2,1))
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         main="", cex.lab =.8,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
 # par(adj = 0)
#title(main= "A")
par(adj=.5)
points(sc_si, 
       col= mycols[metadat2$Treatment],
       pch= c(16,17)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=.8,
       bg=mycols[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 30,
            cex=.8)

# legend("bottomleft", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
#        fill= IBM,
#        cex=.8,
# #        title = "",
# #        bty = "n")
# legend("bottom", legend=c("Nitrogen +", "Nitrogen -"  ),
#       pch=c(16,17 ),
#        cex=.8,
#        title = "",     bty = "n")
# 


###### CAP nitrogen +  and - comparision

ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         cex.lab = .8,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
#par(adj = 0)
#title(main= "B")
#par(adj=.5)
points(sc_si, 
       col= bw[as.factor(metadat2$Nitrogen_label)],
       pch= c(21,24)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=.8,
       bg=bw[as.factor(metadat2$Nitrogen_label)])
ordiellipse(sc_si, metadat2$Nitrogen_label,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            col= bw,
            alpha = 30,
            cex=.8)
# legend("bottomleft", legend=c("Nitrogen -", "Nitrogen +"  ),
#        pch=c(16,17 ),
#        cex=.8,
#        title = "",     bty = "n")

dev.off()

####### 5.  N+ only and N - only  #####
# Constrained ordination
# Perform vegdist analysis of BC distances #

# subset data
ps1 <-subset_samples(ps, Treatment !="Soil" & N=="1" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil"  & N=="1")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
head(metadat2)
# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Treatment+block,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

# # 3. Permutation test for significance of constraints
anova_cap <- anova(cap_result, permutations = 999, by = "term")
anova_cap

# Perform all Pairwise Comparisons
#The function will iterate through all pairs of the 'Habitat' factor

pairwise_results <- multiconstrained(
  formula = dist_matrix ~ Treatment+block,  # Same formula as the main CAP model
  data = metadat2,
  constrained = capscale,       # Specify the constrained ordination method
  permutations = 999            # Number of permutations for the test
)
# View the raw pairwise results
print(pairwise_results)
# Extract the raw p-values from the results
raw_pvalues <- pairwise_results[, "Pr(>F)"]
# Apply the Holm (Holm-Bonferroni) Adjustment
adjusted_pvalues <-p.adjust(raw_pvalues, method = "fdr")
# make table
table <- data.frame(
  Group1 = str_split_i(rownames(pairwise_results), "vs. ", 1),
  Group2 = str_split_i(rownames(pairwise_results), "vs. ", 2),
  Pseudo_F = pairwise_results[, "F"],
  Raw_P = raw_pvalues,
  p_value = round(adjusted_pvalues,4)
)
print(table)
factor(table$Group1)
table$Group1<-factor(table$Group1,levels= c("L ", "G ",  "B ", "GB " , "LB ", "LG " ) )
table$Group1


factor(table$Group2)
table$Group2[order(table$Group2)]
table$Group2
table$Group2<-factor(table$Group2,levels= c( "G", "B", "GB", "LB", "LG", "LGB") )



#Combine the results for final interpretation
# 4) Plot heatmap
p_recol_covercrop <- ggplot(table, aes(x = Group2, y = Group1, fill = p_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = p_value), size = 2) +
  scale_fill_gradient(
    low = "#648FFF",
    high = "white",
    na.value = "grey90",
    limits = c(0, 0.1),
    name = "P-value"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
pdf("Pval_withN.pdf",  width=3.5, height=2.5 )
p_recol_covercrop
dev.off()



### 4. grab info for the plot
smry <- summary(cap_result)
smry
sc_si <- scores(cap_result, display="sites", choices=c(1,2), scaling=1)
sc_si

# Extract the model's adjusted R2
RsquareAdj(cap_result)$adj.r.squared

# percent varience of total varience on RDA 1 and RDA 2
perc <- round(100*(summary(cap_result)$cont$importance[2, 1:2]), 2)
perc


###  5. plot 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
svg("cap.total2.svg", width = 3, height = 7)
#windows(4,7)
par(mfrow=c(2,1))
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         main="Nitrogen + ", cex.lab = .8,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj = 0)
#title(main= "C")
par(adj=.5)
points(sc_si, 
       col= mycols[metadat2$Treatment],
       pch= c(21)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=.8,
       bg=mycols[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 40,
            cex=.8)


####### N- only  
# Constrained ordination
# Perform vegdist analysis of BC distances #

# subset data
ps1 <-subset_samples(ps, Treatment !="Soil" & N=="0" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil"  & N=="0")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Treatment+block,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

# 3. Permutation test for significance of constraints
anova_cap <- anova(cap_result, permutations = 999, by = "term")
anova_cap


### 4. grab info for the plot
smry <- summary(cap_result)
smry
sc_si <- scores(cap_result, display="sites", choices=c(1,2), scaling=1)
sc_si
# Extract the model's adjusted R2
RsquareAdj(cap_result)$adj.r.squared
# percent varience of total varience on RDA 1 and RDA 2
perc <- round(100*(summary(cap_result)$cont$importance[2, 1:2]), 2)
perc

### 5. plot 
#par(cex= 1) # make all fonts in graphs little bigger
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         main="Nitrogen -", cex.lab = .8,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj = 0)
#title(main= "D")
par(adj=.5)
points(sc_si, 
       col= mycols[metadat2$Treatment],
       pch= c(24)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=.8,
       bg=mycols[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 50,
            cex=.8)


dev.off()

#
# Perform all Pairwise Comparisons
#The function will iterate through all pairs of the 'Habitat' factor

pairwise_results <- multiconstrained(
  formula = dist_matrix ~ Treatment+block,  # Same formula as the main CAP model
  data = metadat2,
  constrained = capscale,       # Specify the constrained ordination method
  permutations = 999            # Number of permutations for the test
)
# View the raw pairwise results
print(pairwise_results)
# Extract the raw p-values from the results
raw_pvalues <- pairwise_results[, "Pr(>F)"]
# Apply the Holm (Holm-Bonferroni) Adjustment
adjusted_pvalues <-p.adjust(raw_pvalues, method = "fdr")
# make table
table <- data.frame(
  Group1 = str_split_i(rownames(pairwise_results), "vs. ", 1),
  Group2 = str_split_i(rownames(pairwise_results), "vs. ", 2),
  Pseudo_F = pairwise_results[, "F"],
  Raw_P = raw_pvalues,
  p_value = round(adjusted_pvalues,4)
)
print(table)
factor(table$Group1)
table$Group1<-factor(table$Group1,levels= c("L ", "G ",  "B ", "GB " , "LB ", "LG " ) )
table$Group1


factor(table$Group2)
table$Group2[order(table$Group2)]
table$Group2
table$Group2<-factor(table$Group2,levels= c( "G", "B", "GB", "LB", "LG", "LGB") )



#Combine the results for final interpretation


# 4) Plot heatmap
p_recol_covercrop <- ggplot(table, aes(x = Group2, y = Group1, fill = p_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = p_value), size = 2) +
  scale_fill_gradient(
    low = "#648FFF",
    high = "white",
    na.value = "grey90",
    limits = c(0, 0.1),
    name = "P-value"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
pdf("Pval_withoutN.pdf",  width=3.5, height=2.5)
p_recol_covercrop
dev.off()


#species scores #####
plot(p1.cap, display = c("sites", "species"))
plot(p1.cap, display = c("sp", "wa"))
plot(p1.cap, display = "species", type = "text")
p1.cap$inertia
species_scores <- as.data.frame(scores(
  x = p1.cap,
  display = "species" # or "sp"
))
species_scores[order(species_scores$CAP1),]
head(species_scores)
species_scores$ASV<-row.names(species_scores)

species_scores %>% filter(CAP1< -0.2) %>%
ggplot( mapping = aes(CAP1, CAP2))+
  geom_point()+
  geom_label(aes(label=ASV))

plot(species_scores$CAP1, species_scores$CAP2)
###
species_scores %>% filter(CAP1< -0.4)

#### nitrogen ##
p1 <-plot_ordination(ps1,  p1.cap, color= "N")+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=N), linetype=1)+
  theme(text=element_text(size=15),
        legend.position="left")+
  scale_color_manual(values =  mycols2, name="Nitrogen")+
  scale_shape_discrete(name= "Nitrogen")+
  coord_fixed(1)+
  ggtitle("model: composition * N")
p1

# cap plot total factor
p1 <-plot_ordination(ps1,  p1.cap, color="composition")+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  theme(text=element_text(size=15),
        legend.position="none")+
  scale_color_manual(values =  IBM, name="composition")+
  scale_shape_discrete(name= "Nitrogen")+
  coord_fixed(1/2)+
  ggtitle("model: composition * N")+
  facet_grid(~mixture)
p1

# cap plot total factor
p1 <-plot_ordination(ps1,  p1.cap, color="composition")+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  theme(text=element_text(size=15),
        legend.position="none")+
  scale_color_manual(values =  IBM, name="composition")+
  scale_shape_discrete(name= "Nitrogen")+
  coord_fixed(1/2)+
  ggtitle("model: composition * N")+
  facet_grid(~Legume_label)
p1


#####################CAP pairwise tests ##########

# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Treatment*N*Block,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA




# 3. Permutation test for significance of constraints
anova_cap <- anova(cap_result, permutations = 999, by = "term")
anova_cap
# 6. Perform all Pairwise Comparisons
# The function will iterate through all pairs of the 'Habitat' factor

pairwise_results <- multiconstrained(
  formula = dist_matrix ~ Treatment,  # Same formula as the main CAP model
  data = metadat2,
  constrained = capscale,       # Specify the constrained ordination method
  permutations = 999            # Number of permutations for the test
)

# 7. View the raw pairwise results
print(pairwise_results)

# 8. Extract the raw p-values from the results
raw_pvalues <- pairwise_results[, "Pr(>F)"]

# 9. Apply the Holm (Holm-Bonferroni) Adjustment
adjusted_pvalues_holm <- p.adjust(raw_pvalues, method = "bonferroni")
?p.adjust
# 10. Combine the results for final interpretation
final_table <- data.frame(
  Pair = rownames(pairwise_results),
  Pseudo_F = pairwise_results[, "F"],
  Raw_P = raw_pvalues,
  Holm_Adj_P = adjusted_pvalues_holm
)

# 11. Print the final results table
print(final_table)



###### pairwise adonsis ##########

# A. Install the necessary package (if you haven't already)
# install.packages("devtools") # If you don't have devtools
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


# B. Load the package and run the pairwise test
library(pairwiseAdonis)

# Assuming:
# - 'comm_data' is your community matrix (samples as rows, species as columns)
# - 'env_data' is your environmental/metadata data frame
# - 'Treatment' is the column with your 7 group levels in 'env_data'

  
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# Calculate the Bray-Curtis distance matrix
# Perform vegdist analysis of BC distances #
ps1 <-subset_samples(ps, Treatment !="Soil" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
bray_dist <- vegdist(otu_table(ps1), method = "bray")
# Run the pairwise PERMANOVA
pairwise_results <- pairwise.adonis2(bray_dist ~ Treatment, # Use the distance matrix directly
  data = metadat2,
  permutations = 999,
  method = "bray",
  p.adjust.m = "bonferroni")

pairwise_results
pairwise_results$LG_vs_LB





######## 4. PCOA ##################

#PCOA 
  # Perform vegdist analysis of BC distances #
  ps1 <-subset_samples(ps, Treatment !="Soil" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  # subset metadata
  metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")
  #set factors 
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
  # Calculate Bray-Curtis distance between samples
  otus.bray<-vegdist(otu_table(ps1), method = "bray")
  # Perform PCoA analysis of BC distances #
  otus.pcoa <- cmdscale(otus.bray, k=(44-1), eig=TRUE)
  # Store coordinates for first two axes in new variable #
  otus.p <- otus.pcoa$points[,1:2]
  # Calculate % variance explained by each axis #
  otus.eig<-otus.pcoa$eig
  perc.exp<-otus.eig/(sum(otus.eig))*100
  variance_explained1<-perc.exp[1]
  variance_explained2<-perc.exp[2]
  
#plot

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
#svg("baseRplot.svg", width=10, height = 6)
#windows(8,4)
par(mfrow=c(1,2))
  ordiplot(otus.pcoa,choices=c(1,2), type="none",xlab=paste("PCoA1 (",round(variance_explained1,2),"% variance explained)"),
           ylab=paste("PCoA2 (",round(variance_explained2,2),"% variance explained)"))+
  par(adj = 0)
  title(main= "A")
  par(adj=.5)
  points(otus.p, 
         col= mycols[metadat2$Treatment],
         pch= c(22,24)[as.factor(metadat2$N)],
         lwd=1,cex=1.5,
         bg=mycols[metadat2$Treatment])
  ordiellipse(otus.pcoa, metadat2$Treatment,  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              col= IBM,
              alpha = 30)
  
  legend("topleft", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
         fill= IBM,
         cex=1,
         title = "Composition",
         bty = "n")
  #legend("bottomleft", legend=c("Nitrogen -", "Nitrogen +"  ),
  #       pch=c(22,24 ),
  ##       cex=1,
  #       title = "",     bty = "n")
  

# plot
ordiplot(otus.pcoa,choices=c(1,2), type="none",xlab=paste("PCoA1 (",round(variance_explained1,2),"% variance explained)"),
           ylab=paste("PCoA2 (",round(variance_explained2,2),"% variance explained)"))+
    par(adj = 0)
  title(main= "B")
  par(adj=.5)
points(otus.p, 
       col= mycols2[as.factor(metadat2$N)],
       pch= c(22,24)[as.factor(metadat2$N)],
       lwd=1,cex=1.5,
       bg=mycols2[as.factor(metadat2$N)],)
ordiellipse(otus.pcoa, as.factor(metadat2$N),  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            col= mycols2,
            alpha = 30)  +
  legend("topleft", legend=c("Nitrogen -", "Nitrogen +"  ),
         pch=c(22,24 ),
         cex=1,
         title = "",     bty = "n")
  
  
  dev.off()
  
  
  
  
  
  
  
  
# plot with ggplot
otus.p <-  as.data.frame(otus.p)
colnames(otus.p)<- c("PCoA1", "PCoA2")
df<-cbind(otus.p, metadat2)

# plot
p1<-ggplot(df, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_point(aes(shape = as.factor(N) ), size=3)+
   labs(title = "A",
     x = paste0("PCoA 1 (", round(variance_explained1,2), "%)"),
    y = paste0("PCoA 2 (", round(variance_explained2,2), "%)")) +
  theme_bw() +
  #coord_fixed(1)+  # Ensure the axes are scaled equally (important for ordination plots)
  scale_color_manual(values =  mycols, name="composition")+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  theme(text=element_text(size=12), strip.text.x=element_text(size=12),
        legend.position="none")


# Display the plot
print(p1) 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
svg("pcoa.svg", width=4, height = 4)
p1
dev.off()

   
  
######PCOA nitrogen + no nitrogen##

mycols2 <- c(
  "grey",   # grey
  "black" # teal
)
# plot with ggplot
p2<-ggplot(df, aes(x = PCoA1, y = PCoA2, color = as.factor(N))) +
  geom_point(aes(shape = as.factor(N) ), size=3)+
  labs(title = "B",
       x = paste0("PCoA 1 (", round(variance_explained1,2), "%)"),
       y = paste0("PCoA 2 (", round(variance_explained2,2), "%)")) +
  theme_bw() +
  
  scale_color_manual(values =  mycols2, name="Nitrogen")+
  stat_ellipse(aes(group=as.factor(N)), linetype=1)+
  theme(text=element_text(size=12), strip.text.x=element_text(size=12),
        legend.position="none")

p2
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
svg("Npcoa.svg", width=4, height = 4)
p2
dev.off()


####### PCOA PLOTS facet wrap#####
p3<- df %>%
ggplot( aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_point(aes(shape = as.factor(N) ), size=3)+
  labs(title = "D",
       x = paste0("PCoA 1 (", round(variance_explained1,2), "%)"),
       y = paste0("PCoA 2 (", round(variance_explained2,2), "%)")) +
  theme_bw() +
 # coord_fixed(ratio=1/2)+  # Ensure the axes are scaled equally (important for ordination plots)
  scale_color_manual(values =  mycols, name="composition")+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  theme(text=element_text(size=12), strip.text.x=element_text(size=12),
        legend.position="none")+
  facet_wrap(~Legume_label, ncol=2)
p3

p4<- df %>%
  ggplot( aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_point(aes(shape = as.factor(N) ), size=3)+
  labs(title = "C",
       x = paste0("PCoA 1 (", round(variance_explained1,2), "%)"),
       y = paste0("PCoA 2 (", round(variance_explained2,2), "%)")) +
  theme_bw() +
  # coord_fixed(ratio=1/2)+  # Ensure the axes are scaled equally (important for ordination plots)
  scale_color_manual(values =  mycols, name="composition")+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  theme(text=element_text(size=12), strip.text.x=element_text(size=12),
        legend.position="none")+
  facet_wrap(~mixture, ncol=2)
p4

require(gridExtra)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
svg("facetplot.svg", width=7, height = 7)
grid.arrange(p4, p3, ncol=1)
dev.off()


######PCOA STATS: BETA DISPERSION#
#full dna  between trts
ps2<-subset_samples(ps, Fraction =="Total")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
ps2

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")

# subset metadata
metadat2<-filter(metadat,Fraction =="Total")

#calculate beta dispersion
dispersion <- betadisper(otus.bray, group=metadat2$Treatment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)
#Number of permutations: 999
#
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
##Groups     7 0.040222 0.0057459 8.0956    999  0.001 ***
#  Residuals 78 0.055362 0.0007098  

#Total but not including soil
ps2<-subset_samples(ps, Fraction =="Total" & Treatment!="Soil")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
ps2

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")

# subset metadata
metadat2<-filter(metadat,Fraction =="Total" & Treatment!="Soil")

#calculate beta dispersion
dispersion <- betadisper(otus.bray, group=metadat2$Treatment)
permutest(dispersion, by = "terms")
plot(dispersion, hull=FALSE, ellipse=TRUE)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)   
#Groups     6 0.018732 0.00312194 4.3247    999  0.004 **
# Residuals 76 0.054863 0.00072189   

#full dna  between Fraction
ps2<-subset_samples(ps, Fraction!="CTL")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)

otus.bray<-vegdist(otu_table(ps2), method = "bray")
metadat2<-filter(metadat,Fraction!="CTL")

#calculate beta dispersion
dispersion <- betadisper(otus.bray, group=metadat2$Fraction)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)


##### no total dna -- between Fraction
ps2<-subset_samples(ps, Fraction!="CTL" & Fraction!="Total")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
ps2

otus.bray<-vegdist(otu_table(ps2), method = "bray")
metadat2<-filter(metadat,Fraction!="CTL"& Fraction!="Total")

#calculate beta dispersion
dispersion <- betadisper(otus.bray, group=metadat2$Fraction)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)

######PERMANOVA ######

asvs.clean<-otu_table(ps1)
metadat2<-as.data.frame(as.matrix(sample_data(ps1)))
head(metadat2)
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
head(asvs.bray)
asvs.perm<- adonis2(asvs.clean ~ (N+Treatment)^2, data = metadat2, permutations = 999, method="bray", by="terms")
asvs.perm
summary(asvs.perm)

#Df SumOfSqs      R2      F Pr(>F)    
#Grass             1   0.6725 0.02898 2.5422  0.001 ***
#  Brassicae         1   0.3851 0.01660 1.4558  0.004 ** 
#  Legume            1   1.1009 0.04745 4.1618  0.001 ***
#  Grass:Brassicae   1   0.3533 0.01522 1.3354  0.015 *  
#  Grass:Legume      1   0.3167 0.01365 1.1970  0.040 *  
#  Brassicae:Legume  1   0.2711 0.01168 1.0248  0.325    
#Residual         76  20.1049 0.86642                  
#Total            82  23.2045 1.00000     

adonis2(formula = asvs.clean ~ (Grass*Brassicae*Legume + Grass:Brassicae:Legume), data = metadat2, permutations = 999, method = "bray")


#adonis2(formula = asvs.clean ~ (Grass * Brassicae * Legume + Grass:Brassicae:Legume), data = metadat2, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Grass             1   0.6725 0.02898 2.5422  0.001 ***
#  Brassicae         1   0.3851 0.01660 1.4558  0.006 ** 
#  Legume            1   1.1009 0.04745 4.1618  0.001 ***
#  Grass:Brassicae   1   0.3533 0.01522 1.3354  0.013 *  
#  Grass:Legume      1   0.3167 0.01365 1.1970  0.056 .  
# Brassicae:Legume  1   0.2711 0.01168 1.0248  0.309    
#Residual         76  20.1049 0.86642                  
#Total            82  23.2045 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(formula = asvs.clean ~ (N*Treatment + Grass:Brassicae:Legume), data = metadat2, permutations = 999, method = "bray")

######PERMANOVA by Treatment##

t1<-subset_samples(total, n_species=="1" & Legume=="0")
asvs.clean<-otu_table(t1)
metadat2<-as.data.frame(as.matrix(sample_data(t1)))
head(metadat2)
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
head(asvs.bray)
asvs.perm<- adonis2(asvs.clean ~ Treatment*N, data = metadat2, permutations = 999, method="bray")
asvs.perm
# G and B are different

# 2 species
t1<-subset_samples(total, n_species=="2")
asvs.clean<-otu_table(t1)
metadat2<-as.data.frame(as.matrix(sample_data(t1)))
metadat2
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
head(asvs.bray)
asvs.perm <- adonis2(asvs.clean ~ Treatment*N, data = metadat2, permutations = 999, method="bray")
asvs.perm
# LG and GB are different

# 2 species
t1<-subset_samples(total, n_species!="1", Treatment!="GB"  & Treatment!="LB")
asvs.clean<-otu_table(t1)
metadat2<-as.data.frame(as.matrix(sample_data(t1)))
metadat2
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
head(asvs.bray)
asvs.perm <- adonis2(asvs.clean ~ Treatment*N, data = metadat2, permutations = 999, method="bray")
asvs.perm
# LG and GB are different



#L treatments
# 2 species
t1<-subset_samples(total, Legume=="1", Treatment!="L")
asvs.clean<-otu_table(t1)
metadat2<-as.data.frame(as.matrix(sample_data(t1)))
metadat2
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
head(asvs.bray)
asvs.perm <- adonis2(asvs.clean ~ Treatment*N, data = metadat2, permutations = 999, method="bray")
asvs.perm
# LG and GB are different



######PERMANOVA by Fraction###
ps2<-subset_samples(ps, Fraction!="CTL")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
#
#get asvs table
asvs.clean<-otu_table(ps2)
metadat2<-filter(metadat,Fraction!="CTL")
#
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
#
asvs.perm<- adonis2(asvs.clean ~ Fraction, data = metadat2, permutations = 999, method="bray")
asvs.perm

##active verse inactive
ps2<-subset_samples(ps, Fraction!="CTL" & Fraction!="Total")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
#
#get asvs table
asvs.clean<-otu_table(ps2)
metadat2<-filter(metadat,Fraction!="CTL" & Fraction!="Total")
#
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
#
asvs.perm<- adonis2(asvs.clean ~ Fraction, data = metadat2, permutations = 999, method="bray")
asvs.perm
#adonis2(formula = asvs.clean ~ Fraction, data = metadat2, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Fraction  1   1.0113 0.03579 2.7098  0.001 ***


######## 6. PERCENT abundance figure ##################
#get data
#otu_table(total)[1:5, 1:5]

# Transform ASV counts to proportions / relative abundance

total <-subset_samples(total, Treatment !="Soil" )
total<-prune_taxa(taxa_sums(total) > 0, total)
total

total.prop <- transform_sample_counts(subset_samples(total),function(ASV) ASV/sum(ASV))
otu_table(total.prop)[1:5, 1:5]

# remove taxa groups finer than Phyla
total.prop <- phyloseq::tax_glom(total.prop, "Phyla")
# check number of Phylas and names
length(unique(as.data.frame(tax_table(total.prop))$Phyla))
#
#head(as.data.frame(tax_table(total.prop)))
# 23 Phyla

# Melt data frame with phyloseq function for plotting
df <- psmelt(total.prop)
head(df)
#
#aggregate by Phyla
df <- df %>% group_by(across(c(-OTU, -Feature.ID))) %>%
summarise(Abundance= sum(Abundance))
#
#rename low abundance taxa as other
df.mean<-df %>%
group_by(Phyla) %>%
summarise(mean=mean(Abundance)) %>%
arrange(.,mean)
df.mean

remove<-filter(df.mean, mean<.0005 )$Phyla
#for anything in this list, remove it and call it other

for (i in remove) {
    print(i)
    df$Phyla<-sub(i, "other", df$Phyla)
  
    }

#df$Phy <- "other"
# check Phyla names
unique(df$Phyla)
# check number of Phylas
length(unique(df$Phyla))
#18
# remove the o__ in front of Phyla
df$Phyla<-sub(" p__", "", df$Phyla)
unique(df$Phyla)
# Reorder levels to put other at the end - otherwise taxa are in alphabetical Phyla
#install.packages("forcats")
library(forcats)
df$Phyla<-fct_relevel(as.factor(df$Phyla), "other", after = Inf)
# make really low abundance taxa other
unique(df$Phyla)
colnames(df)

df$Treatment <- factor(df$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))



#set working directory for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")

# Visualize

svg("percent.abund.total.svg", width = 12, height = 7)
windows(12,7)
ggplot(df, aes(x=Sample, y=Abundance, fill=Phyla))+
  geom_bar(stat = "identity")+
  facet_wrap(~Treatment, scales="free_x", nrow = 1)+
  #           labeller=as_labeller(c(BC='Bush clover', YSC='Yellow sweet clover',
  #                                  FP='Field pea', CC='Crimson clover',
  #                                  WC='White clover', CP='Cowpea')))+
  #labs(y="Relative abundance", title="Nodule")+
  theme_classic()+
  theme(legend.position='right', axis.title.x=element_blank(),
        axis.text.x=element_blank(), text=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=16), legend.text=element_text
        (size=14))+
  ggtitle("Total DNA")+ 
  scale_fill_manual(values=c("#855C75", "#dba7c7",
                             "#996e29", "#D9AF6B",
                             "#AF6458", "#eba69b" , 
                             "#736F4C", "#bdb993", 
                             "#455669", "#8095ab",
                             "#625377", "#baafc9",
                             "#68855C", "#c7dbbf", 
                             #"#9C9C5E", "#e3e3bc",
                             "#467378", "#b9c8c9",
                             "#8C785D", "#ded6cc",
                             "#7C7C7C", "#C3C3C3"))
dev.off()



####old ##
taxon<- as.data.frame(tax_table(ps))
df<-as.data.frame(otu_table(ps))

# make it percent
df<-(df/rowSums(df))*100
df<-as.data.frame(t(df))
df.taxa<-cbind(df, taxon)


# summarize by phyla
df1<-aggregate(cbind( BCAT_11_S31 ,BCAT_12_S41 ,  BCAT_13_S51 , BCAT_14_S61 ,  BCAT_15_S71 , BCAT_23_S2 , 
                      BCAT_24_S12 , BCAT_25_S22 ,  BCAT_26_S32 ,  BCAT_28_S52 , BCAT_30_S62 ,  BCAT_38_S72 ,
                      BCAT_39_S3  , BCAT_40_S13 ,  BCAT_41_S23 ,  BCAT_43_S33 ,  BCAT_44_S43 ,  BCAT_45_S53 , 
                      BCAT_53_S63 , BCAT_54_S73 ,  BCAT_55_S4  ,  BCAT_56_S14  ,  BCAT_57_S24  ,  BCAT_59_S44 ,  
                      BCAT_60_S54 , BCAT_68_S64 ,  BCAT_69_S74 ,   BCAT_70_S5   ,   BCAT_71_S15  ,  BCAT_72_S25 ,  
                      BCAT_73_S35 , BCAT_8_S1   ,  BCAT_83_S45 ,   BCAT_84_S55  ,  BCAT_86_S65  ,  BCAT_87_S75 ,  
                      BCAT_88_S6  , BCAT_9_S11  , fc_ct_iso_S83,   fc_ct_sh_S84 ,
                      i_10_S36 ,      i_11_S46 ,  
                      i_12_S56 ,    i_14_S76   ,    i_15_S7   ,    i_23_S17   ,    i_25_S37 ,    i_26_S47  ,   
                      i_27_S57 ,    i_30_S77   ,    i_38_S8   ,    i_39_S18   ,   i_40_S28  ,    i_41_S38  , 
                      i_43_S48 ,    i_44_S58   ,    i_45_S68  ,    i_53_S78   ,    i_54_S9  ,    i_55_S19  ,   
                      i_56_S29 ,    i_57_S39  ,    i_58_S49  ,    i_59_S59   ,    i_60_S69 ,    i_68_S79  ,    
                      i_69_S10 ,    i_70_S20  ,    i_71_S30  ,    i_72_S40   ,    i_73_S50 ,    i_8_S16   ,    
                      i_83_S60 ,    i_84_S70  ,    i_86_S80  ,   i_88_S81    ,   i_9_S26   ,   pcr_ctl25_S173 ,
                      pcr_ctl35_S82,  T_DNA_1_S85,
                      T_DNA_10_S97  , T_DNA_11_S108 , T_DNA_12_S119 , T_DNA_13_S130 ,
                      T_DNA_14_S141 , T_DNA_15_S152 , T_DNA_16_S163 , T_DNA_17_S87 ,  T_DNA_18_S98  ,  T_DNA_19_S109 ,
                      T_DNA_2_S96  ,  T_DNA_20_S120 , T_DNA_21_S131 , T_DNA_22_S142 , T_DNA_24_S164 , T_DNA_25_S88  ,
                      T_DNA_26_S99  , T_DNA_27_S110  ,T_DNA_28_S121 , T_DNA_29_S132 , T_DNA_3_S107  , T_DNA_31_S154 ,
                      T_DNA_32_S165 , T_DNA_33_S89  , T_DNA_34_S100 , T_DNA_35_S111 , T_DNA_36_S122 , T_DNA_37_S133 ,
                      T_DNA_38_S144 , T_DNA_39_S155 , T_DNA_4_S118  , T_DNA_40_S166 , T_DNA_41_S90  , T_DNA_42_S101 ,
                      T_DNA_43_S112 , T_DNA_44_S123 , T_DNA_45_S134 , T_DNA_46_S145 , T_DNA_47_S156 , T_DNA_48_S167 ,
                      T_DNA_49_S91  , T_DNA_5_S129  , T_DNA_50_S102 , T_DNA_51_S113 , T_DNA_52_S124 , T_DNA_53_S135 ,
                      T_DNA_54_S146 , T_DNA_55_S157 , T_DNA_56_S168 , T_DNA_57_S92  , T_DNA_58_S103 , T_DNA_59_S114 ,
                      T_DNA_6_S140  , T_DNA_60_S125 , T_DNA_61_S136 , T_DNA_62_S147 , T_DNA_63_S158 , T_DNA_64_S169 ,
                      T_DNA_65_S93  , T_DNA_66_S104 , T_DNA_67_S115 , T_DNA_68_S126 , T_DNA_69_S137 , T_DNA_7_S151  ,
                      T_DNA_70_S148 , T_DNA_71_S159 , T_DNA_72_S170 , T_DNA_73_S94  , T_DNA_74_S105 , T_DNA_75_S116 ,
                      T_DNA_76_S127 , T_DNA_77_S138 , T_DNA_78_S149 , T_DNA_79_S160 , T_DNA_8_S162  , T_DNA_80_S171 ,
                      T_DNA_81_S95  , T_DNA_82_S106 , T_DNA_83_S117 , T_DNA_84_S128 , T_DNA_85_S139 , T_DNA_86_S150 ,
                      T_DNA_87_S161 , T_DNA_88_S172 , T_DNA_9_S86 
                      
)~ Phyla, data = df.taxa, FUN = sum, na.rm = TRUE)

#head(df1)
# summ row 1 and 2 b\c they are both unassigned taxa
n<-dim(df1)[2]
#n
row1<-df1[1,2:n]+ df1[2,2:n] 
# call empty phyla unassigned
row1<-c("Unassigned", row1)
# put in df
#row1

df1[1,] <- as.vector(row1)
df1<-df1[c(-2),] #remove p___ row because it is already counted

# summarize by treatment
row.names(df1)<-df1$Phyla
df1<-df1[,-1]                        # first row is colnames names -- re name columnes. 
#head(df1)
df1<-t(df1) #transform
#head(df1)
df1<-as.data.frame(df1)
head(df1)
#df1$Trt_fraction <- metadat$Trt_fraction

df1<-metadat %>%
select(Trt_fraction, n_species, Treatment, Fraction) %>%
  cbind(., df1)


head(df1)

df1<-df1%>%
group_by(Trt_fraction, Fraction, n_species, Treatment) %>%
summarise_all(mean)

# gather by sample

df<- df1%>%
 gather("Phyla", value, 5:45 ) %>%
 filter(. , value>0)

df


# make really low abundance taxa other
df$Phyla[df$value<1] <- "other"
df0<-aggregate(cbind(value) ~ Trt_fraction+Phyla, data = df, FUN = sum, na.rm =TRUE)
head(df0)
#df0<-df0[order(df0$Trt_fraction),]
#head(df0)

hist(df0$value)
unique(df0$Phyla)
RColorBrewer::brewer.pal(12, "Paired")

mycols18<-c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
"#6A3D9A", "#FFFF99", "#B15928", "#5c3218")

mycols18<- c( "#1F78B4","#A6CEE3","#E31A1C",  "#FB9A99", "#eb05db",  "#ffccef", "#33A02C","#B2DF8A",  "#FF7F00",  "#FDBF6F", "#6A3D9A" , "#CAB2D6",
            "grey")

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/")
svg(file="percent_barplot.svg",width = 12, height=10)
windows(12,12)
df0%>% 
  filter(Trt_fraction!="CTL_CTL") %>%
  ggplot(aes(fill=Phyla, y=value, x=Trt_fraction)) + 
  geom_bar(position="fill", stat= "identity")+
 scale_fill_manual(values=mycols18) +
 # scale_fill_viridis_d() +
  #ggtitle("Top phyla") +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()





####PERCENT abundance figure: ####
# grab data
taxon<- as.data.frame(tax_table(ps))
df<-as.data.frame(otu_table(ps))

# make it percent
df<-(df/rowSums(df))*100
df<-as.data.frame(t(df))
df.taxa<-cbind(df, taxon)

# rename columns if you want, not required. 
colnames(df)
length(n)
length(colnames(df))
n<-c("BCAT_11_S31" ,"BCAT_12_S41" ,  "BCAT_13_S51" , "BCAT_14_S61" ,  "BCAT_15_S71" , "BCAT_23_S2" , 
  "BCAT_24_S12" , "BCAT_25_S22" ,  "BCAT_26_S32" ,  "BCAT_28_S52" , "BCAT_30_S62" ,  "BCAT_38_S72" ,
  "BCAT_39_S3"  , "BCAT_40_S13" ,  "BCAT_41_S23" ,  "BCAT_43_S33" ,  "BCAT_44_S43" ,  "BCAT_45_S53" , 
  "BCAT_53_S63" , "BCAT_54_S73" ,  "BCAT_55_S4"  ,  "BCAT_56_S14"  ,  "BCAT_57_S24"  ,  "BCAT_59_S44" ,  
  "BCAT_60_S54" , "BCAT_68_S64" ,  "BCAT_69_S74" ,   "BCAT_70_S5"   ,   "BCAT_71_S15"  ,  "BCAT_72_S25" ,  
  "BCAT_73_S35" , "BCAT_8_S1"   ,  "BCAT_83_S45" ,   "BCAT_84_S55"  ,  "BCAT_86_S65"  ,  "BCAT_87_S75" ,  
  "BCAT_88_S6"  , "BCAT_9_S11"  , "fc_ct_iso_S83",   "fc_ct_sh_S84" ,
  "i_10_S36" ,      "i_11_S46" ,  
  "i_12_S56" ,    "i_14_S76"   ,    "i_15_S7"   ,    "i_23_S17"   ,    "i_25_S37" ,    "i_26_S47"  ,   
  "i_27_S57" ,    "i_30_S77"   ,    "i_38_S8"   ,    "i_39_S18"   ,   "i_40_S28"  ,    "i_41_S38"  , 
  "i_43_S48" ,    "i_44_S58"   ,    "i_45_S68"  ,    "i_53_S78"   ,    "i_54_S9"  ,    "i_55_S19"  ,   
  "i_56_S29" ,    "i_57_S39"  ,    "i_58_S49"  ,    "i_59_S59"   ,    "i_60_S69" ,    "i_68_S79"  ,    
  "i_69_S10" ,    "i_70_S20"  ,    "i_71_S30"  ,    "i_72_S40"   ,    "i_73_S50" ,    "i_8_S16"   ,    
  "i_83_S60" ,    "i_84_S70"  ,    "i_86_S80"  ,   "i_88_S81"    ,   "i_9_S26"   ,   "pcr_ctl25_S173" ,
  "pcr_ctl35_S82",  "T_DNA_1_S85",
  "T_DNA_10_S97"  , "T_DNA_11_S108" , "T_DNA_12_S119" , "T_DNA_13_S130" ,
  "T_DNA_14_S141" , "T_DNA_15_S152" , "T_DNA_16_S163" , "T_DNA_17_S87" ,  "T_DNA_18_S98"  ,  "T_DNA_19_S109" ,
  "T_DNA_2_S96"  ,  "T_DNA_20_S120" , "T_DNA_21_S131" , "T_DNA_22_S142" , "T_DNA_24_S164" , "T_DNA_25_S88"  ,
  "T_DNA_26_S99"  , "T_DNA_27_S110"  ,"T_DNA_28_S121" , "T_DNA_29_S132" , "T_DNA_3_S107"  , "T_DNA_31_S154" ,
  "T_DNA_32_S165" , "T_DNA_33_S89"  , "T_DNA_34_S100" , "T_DNA_35_S111" , "T_DNA_36_S122" , "T_DNA_37_S133" ,
  "T_DNA_38_S144" , "T_DNA_39_S155" , "T_DNA_4_S118"  , "T_DNA_40_S166" , "T_DNA_41_S90"  , "T_DNA_42_S101" ,
  "T_DNA_43_S112" , "T_DNA_44_S123" , "T_DNA_45_S134" , "T_DNA_46_S145" , "T_DNA_47_S156" , "T_DNA_48_S167" ,
  "T_DNA_49_S91"  , "T_DNA_5_S129"  , "T_DNA_50_S102" , "T_DNA_51_S113" , "T_DNA_52_S124" , "T_DNA_53_S135" ,
  "T_DNA_54_S146" , "T_DNA_55_S157" , "T_DNA_56_S168" , "T_DNA_57_S92"  , "T_DNA_58_S103" , "T_DNA_59_S114" ,
  "T_DNA_6_S140"  , "T_DNA_60_S125" , "T_DNA_61_S136" , "T_DNA_62_S147" , "T_DNA_63_S158" , "T_DNA_64_S169" ,
  "T_DNA_65_S93"  , "T_DNA_66_S104" , "T_DNA_67_S115" , "T_DNA_68_S126" , "T_DNA_69_S137" , "T_DNA_7_S151"  ,
  "T_DNA_70_S148" , "T_DNA_71_S159" , "T_DNA_72_S170" , "T_DNA_73_S94"  , "T_DNA_74_S105" , "T_DNA_75_S116" ,
  "T_DNA_76_S127" , "T_DNA_77_S138" , "T_DNA_78_S149" , "T_DNA_79_S160" , "T_DNA_8_S162"  , "T_DNA_80_S171" ,
  "T_DNA_81_S95"  , "T_DNA_82_S106" , "T_DNA_83_S117" , "T_DNA_84_S128" , "T_DNA_85_S139" , "T_DNA_86_S150" ,
  "T_DNA_87_S161" , "T_DNA_88_S172" , "T_DNA_9_S86"    )
 
#length(n)
#colnames(df)<-n
#rownames(df)<-NULL
# make rownames null
# summarize by phyla
df1<-aggregate(cbind( BCAT_11_S31 ,BCAT_12_S41 ,  BCAT_13_S51 , BCAT_14_S61 ,  BCAT_15_S71 , BCAT_23_S2 , 
                      BCAT_24_S12 , BCAT_25_S22 ,  BCAT_26_S32 ,  BCAT_28_S52 , BCAT_30_S62 ,  BCAT_38_S72 ,
                      BCAT_39_S3  , BCAT_40_S13 ,  BCAT_41_S23 ,  BCAT_43_S33 ,  BCAT_44_S43 ,  BCAT_45_S53 , 
                      BCAT_53_S63 , BCAT_54_S73 ,  BCAT_55_S4  ,  BCAT_56_S14  ,  BCAT_57_S24  ,  BCAT_59_S44 ,  
                      BCAT_60_S54 , BCAT_68_S64 ,  BCAT_69_S74 ,   BCAT_70_S5   ,   BCAT_71_S15  ,  BCAT_72_S25 ,  
                      BCAT_73_S35 , BCAT_8_S1   ,  BCAT_83_S45 ,   BCAT_84_S55  ,  BCAT_86_S65  ,  BCAT_87_S75 ,  
                      BCAT_88_S6  , BCAT_9_S11  , fc_ct_iso_S83,   fc_ct_sh_S84 ,
                      i_10_S36 ,      i_11_S46 ,  
                      i_12_S56 ,    i_14_S76   ,    i_15_S7   ,    i_23_S17   ,    i_25_S37 ,    i_26_S47  ,   
                      i_27_S57 ,    i_30_S77   ,    i_38_S8   ,    i_39_S18   ,   i_40_S28  ,    i_41_S38  , 
                      i_43_S48 ,    i_44_S58   ,    i_45_S68  ,    i_53_S78   ,    i_54_S9  ,    i_55_S19  ,   
                      i_56_S29 ,    i_57_S39  ,    i_58_S49  ,    i_59_S59   ,    i_60_S69 ,    i_68_S79  ,    
                      i_69_S10 ,    i_70_S20  ,    i_71_S30  ,    i_72_S40   ,    i_73_S50 ,    i_8_S16   ,    
                      i_83_S60 ,    i_84_S70  ,    i_86_S80  ,   i_88_S81    ,   i_9_S26   ,   pcr_ctl25_S173 ,
                      pcr_ctl35_S82,  T_DNA_1_S85,
                      T_DNA_10_S97  , T_DNA_11_S108 , T_DNA_12_S119 , T_DNA_13_S130 ,
                      T_DNA_14_S141 , T_DNA_15_S152 , T_DNA_16_S163 , T_DNA_17_S87 ,  T_DNA_18_S98  ,  T_DNA_19_S109 ,
                      T_DNA_2_S96  ,  T_DNA_20_S120 , T_DNA_21_S131 , T_DNA_22_S142 , T_DNA_24_S164 , T_DNA_25_S88  ,
                      T_DNA_26_S99  , T_DNA_27_S110  ,T_DNA_28_S121 , T_DNA_29_S132 , T_DNA_3_S107  , T_DNA_31_S154 ,
                      T_DNA_32_S165 , T_DNA_33_S89  , T_DNA_34_S100 , T_DNA_35_S111 , T_DNA_36_S122 , T_DNA_37_S133 ,
                      T_DNA_38_S144 , T_DNA_39_S155 , T_DNA_4_S118  , T_DNA_40_S166 , T_DNA_41_S90  , T_DNA_42_S101 ,
                      T_DNA_43_S112 , T_DNA_44_S123 , T_DNA_45_S134 , T_DNA_46_S145 , T_DNA_47_S156 , T_DNA_48_S167 ,
                      T_DNA_49_S91  , T_DNA_5_S129  , T_DNA_50_S102 , T_DNA_51_S113 , T_DNA_52_S124 , T_DNA_53_S135 ,
                      T_DNA_54_S146 , T_DNA_55_S157 , T_DNA_56_S168 , T_DNA_57_S92  , T_DNA_58_S103 , T_DNA_59_S114 ,
                      T_DNA_6_S140  , T_DNA_60_S125 , T_DNA_61_S136 , T_DNA_62_S147 , T_DNA_63_S158 , T_DNA_64_S169 ,
                      T_DNA_65_S93  , T_DNA_66_S104 , T_DNA_67_S115 , T_DNA_68_S126 , T_DNA_69_S137 , T_DNA_7_S151  ,
                      T_DNA_70_S148 , T_DNA_71_S159 , T_DNA_72_S170 , T_DNA_73_S94  , T_DNA_74_S105 , T_DNA_75_S116 ,
                      T_DNA_76_S127 , T_DNA_77_S138 , T_DNA_78_S149 , T_DNA_79_S160 , T_DNA_8_S162  , T_DNA_80_S171 ,
                      T_DNA_81_S95  , T_DNA_82_S106 , T_DNA_83_S117 , T_DNA_84_S128 , T_DNA_85_S139 , T_DNA_86_S150 ,
                      T_DNA_87_S161 , T_DNA_88_S172 , T_DNA_9_S86 
  
               )~ Phyla, data = df.taxa, FUN = sum, na.rm = TRUE)

head(df1)
# summ row 1 and 2 b\c they are both unassigned taxa
n<-dim(df1)[2]
n
row1<-df1[1,2:n]+ df1[2,2:n] 
# call empty phyla unassigned
row1<-c("Unassigned", row1)
# put in df
row1
row1<-as.vector(row1)
df1[1,] <- row1
df1<-df1[c(-2),] #remove p___ row because it is already counted

# gather by sample
df1<-  gather(df1, "sample", value, 2:n )
head(df1)
#remove zeros
df1<-df1[df1$value!=0,]
head(df1)


# make really low abundance taxa other
df1$Phyla[df1$value<3] <- "other"
df0<-aggregate(cbind(value) ~ sample+Phyla, data = df1, FUN = sum, na.rm =TRUE)
head(df0)
df0<-df0[order(df0$sample),]
head(df0)

hist(df0$value)
RColorBrewer::brewer.pal(26, "Spectral")
mycols18<- c( "#1F78B4","#A6CEE3","#E31A1C",  "#FB9A99", "#33A02C","#B2DF8A",  "#FF7F00",  "#FDBF6F", "#6A3D9A" , "#CAB2D6",
               "#B15928", "#FFFF99",  "#eb05db","#edceeb","#1a635a","#9ad6ce" , "#969696", "#232423")

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/")
svg(file="percent_barplot.svg",width = 12, height=10)
windows(12,12)
df1%>% 
  ggplot(aes(fill=Phyla, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values=mycols18) +
  #scale_fill_viridis(discrete = TRUE) +
  #ggtitle("Top phyla") +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

dev.off()


######pc1 cor ####

# filter data 
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps
# filter for active fraction 
ps1 <-subset_samples(ps, Fraction=="Total" & Treatment!="Soil" & Treatment!="CTL")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1 #1778 taxa
# filter metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil" & Treatment!="CTL")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))


# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# 1. Perform PCA
# scale. = TRUE is critical to ensure all variables are treated equally
pca_result <- prcomp(dist_matrix, center = TRUE, scale. = TRUE)

# 2. Extract the first principal component (PC1)
# This is your new one-dimensional variable
pc1_variable <- pca_result$x[, 1]

# get variables of interest - weed seed decay, biomass, nfix, activity
# clean up metadat2
metadat2<-metadat2 %>% select(SampleID, Trt_ID, Pot_ID, Treatment, Rep, Block)

# biomass
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("biomass_potlevel.csv", row.names = 1)
biomass<-biomass %>% select(Trt_ID, Total.Root.g, Stem.Biomass.g ) %>%
  mutate(Root.Biomass=Total.Root.g, Shoot.Biomass= Stem.Biomass.g) %>%
  select(Trt_ID, Root.Biomass, Shoot.Biomass) 
df<-left_join(metadat2, biomass)

# fc 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/")
fc <-read.csv("processed_flow_cyto.csv")
head(fc)
fc<-fc %>% select(Trt_ID, boncat_freq, active_cel_per_g )
df<-left_join(df, fc)

# n fix 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
nfix<-read.csv("Nfix.csv")
head(nfix)
nfix<-nfix %>% 
  filter(duplicate=="N") %>%select(Trt_ID, perc.Ndfa, n_fix_per_legume )
df<-left_join(df, nfix)

# weeds
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed<-read.csv("weed_seed_decay.csv", row.names = 1)
head(weed)
weed<-weed %>% select(Trt_ID, foxtail_prop_nongerm, pigweed_prop_nongerm, foxtail_num_nongerm )
df<-left_join(df, weed)

# z transform
df$z_Root.Biomass<-as.vector(scale(df$Root.Biomass))
df$z_Shoot.Biomass<-as.vector(scale(df$Shoot.Biomass))
df$z_foxtail_prop_nongerm<-as.vector(scale(log(df$foxtail_prop_nongerm+1)))
df$z_foxtail_num_germ<-as.vector(scale(log(df$foxtail_num_nongerm+1)))
df$z_pigweed_prop_nongerm<-as.vector(scale(df$pigweed_prop_nongerm))
df$z_boncat_freq<-as.vector(scale(df$boncat_freq))
df$z_active_cel_per_g<-as.vector(scale(log(df$active_cel_per_g+1)))
df$z_perc.Ndfa<-as.vector(scale(log(df$perc.Ndfa+1)))
df$z_n_fix_per_legume<-as.vector(scale(df$n_fix_per_legume))

#hist
hist(df$z_Root.Biomass) # okay 
hist(df$z_Shoot.Biomass) # okay 
hist(df$z_pigweed_prop_nongerm) # okay
hist(df$z_foxtail_num_germ) # okay
hist(df$z_foxtail_prop_nongerm) # bad
hist(df$z_boncat_freq) # okay
hist(df$z_active_cel_per_g) # okay after log transform
hist(df$z_perc.Ndfa) # kinda skewed but better after log transform
hist(df$z_n_fix_per_legume)

# 4. Add it to your data frame
df$PC1 <- pc1_variable


# corrplots\

par(mfrow=c(2,2))
plot(df$PC1, df$z_Root.Biomass, main= "Total, p=.001, Rsq =.18")
plot(df$PC1, df$z_Shoot.Biomass, main= "Total, p=.001, Rsq =.59")

#par(mfrow=c(1,2))
plot(df$PC1, df$z_boncat_freq, main="total, not sig")
plot(df$PC1, df$z_active_cel_per_g, main= "Total, p=.05, Rsq=.14")

plot(df$PC1, df$z_foxtail_num_germ)
plot(df$PC1, df$z_pigweed_prop_nongerm, )

plot(df$PC1, df$z_perc.Ndfa)
plot(df$PC1, df$z_n_fix_per_legume)

# lm 
m1<-lm(df$PC1~df$z_Root.Biomass)
summary(m1) # sig
m1<-lm(df$PC1~df$z_Shoot.Biomass)
summary(m1) # sig


m1<-lm(df$PC1~df$z_boncat_freq)
summary(m1)
m1<-lm(df$PC1~df$z_active_cel_per_g)
summary(m1) # sig


m1<-lm(df$PC1~df$z_foxtail_prop_nongerm)
summary(m1)
m1<-lm(df$PC1~df$z_pigweed_prop_nongerm)
summary(m1) 

m1<-lm(df$PC1~df$z_perc.Ndfa)
summary(m1) 

m1<-lm(df$PC1~df$z_n_fix_per_legume)
summary(m1) 
plot(m1)


#
###### extract key taxa ####

#  Constrained ordination
ps<-subset_samples(ps,  Treatment!="Soil" & Treatment!="CTL")

# subset data
ps <-subset_samples(ps, N=="1" )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil"  & N=="1" & Treatment!="CTL")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps), method = "bray")

commdata = as.data.frame(otu_table(ps))
# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Treatment,
                       data = metadat2,
                       comm = commdata,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

# look for key taxa
plot(cap_result, display = c("sites", "species"))
plot(cap_result, display = "species")
species_scores <- as.data.frame(scores(
  x = cap_result,
  display = "species" # or "sp"
))

head(species_scores)
species_scores$asv<-row.names(species_scores)

# Calculate vector length (distance from origin)
species_scores <- species_scores %>%
  mutate(dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist))

# Select the top 10 species
top_spp <- head(species_scores, 15)
top_spp

ggplot() +
  # Draw a circle/origin cross for reference
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
  
  # Add the vectors (arrows)
  geom_segment(data = top_spp, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkred") +
  
  # Add labels with some padding
  geom_text(data = top_spp, 
            aes(x = CAP1, y = CAP2, label = asv), 
            color = "black", fontface = "italic", vjust = -0.5) +
  
  theme_bw() +
  labs(title = "Top 10 Species Contributing to CAP Variation",
       x = "CAP1", y = "CAP2")


### genus overall ########
# Create a vector of the genera you want and Subset the phyloseq object
asvkp<-unique(top_spp$asv)
tax<-as.data.frame((tax_table(ps)))
tax <- tax[tax$asv %in% asvkp, ]
tax<-tax %>% group_by(Phyla, Order, Family, Genus, asv) %>% summarise()
tax
target<-unique(tax$Genus)
ps1 <- subset_taxa(ps, Genus %in% target)


# make df relative abundance
df<-as.data.frame((otu_table(ps1)))
df<-df/rowSums(df)
df<- as.data.frame(t(df))
df$asv<-row.names(df)
df

# add taxa info
tax<-as.data.frame((tax_table(ps1)))
df<-left_join(df, tax)
df

#summarize by genus
df<-df %>% group_by(Genus, Phyla) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
df<-df %>% filter(Genus!="" & Genus!=" g__" & Genus!= " g_")
df

# save genus info
Genus <- df$Genus
Phyla <- df$Phyla
df$Genus = NULL
df$Phyla = NULL
df
# add treatment info
df<-as.data.frame(t(df))
colnames(df) <- Genus
df

metadat2<-metadat2 %>% select(Treatment)
df<-cbind(df, metadat2)
df

# pivot
df<-df %>%
  pivot_longer(cols = starts_with(" g_"), 
               names_to = "Genus",
               values_to = "percent_abundance")
df

# add phyla info
tax<-data.frame( Genus = Genus,
            Phyla = Phyla)
df<-left_join(df, tax)
df
##
df$Genus<-gsub(" g__", "", df$Genus)
df$Genus<-gsub("_501058", "", df$Genus)
df$Genus<-gsub("_487784", "", df$Genus)
df$Genus<-gsub("_C", "", df$Genus)


# for plotting?
#df$Genus<-gsub(" g__SZUA-359", " g__Saccharimonadales", df$Genus)
#df$Genus<-gsub(" g__UBA1020", " g__Saccharimonadales", df$Genus)
# plot
mycols<-c("blue", "purple","#B2DF8A", "#FF7F00")
df$percent_abundance <- as.numeric(df$percent_abundance) 
df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")

svg("total.genus1.svg", height = 4, width=8)
df %>%
  ggplot(aes(x=Treatment, y=percent_abundance, fill = Phyla))+
  geom_boxplot(outliers=FALSE)+
  geom_jitter()+
  scale_fill_manual(values= mycols)+
  theme_classic(base_size = 12) +
  facet_grid(~Genus, scales="free", space="free")+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5) #,
        #legend.position = "none"
        )+
 labs(y="Relative Abundance")
#coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
dev.off()

######### test###############
n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Caulobacter")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

df1<-df %>% filter(Genus=="Nitrosocosmicus") #
kruskal.test(percent_abundance ~ Treatment, data = df1) 
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value

df1<-df %>% filter(Genus=="Novosphingobium") # 
kruskal.test(percent_abundance ~ Treatment, data = df1) # 
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value

df1<-df %>% filter(Genus=="Polaromonas")
kruskal.test(percent_abundance ~ Treatment, data = df1) # 
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n4<-m1$p.value

df1<-df %>% filter(Genus=="Rhizobium")
kruskal.test(percent_abundance ~ Treatment, data = df1) # sig
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n5<-m1$p.value

df1<-df %>% filter(Genus=="Saccharimonadales")
kruskal.test(percent_abundance ~ Treatment, data = df1) # sig
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n6<-m1$p.value


# raw results
raw_pvalues<-c( n1, n2, n3, n4, n5, n6)
taxa <-  n
#  # 9. Apply the Holm (Holm-Bonferroni) Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues,
  fdr_adj_p_nonsci = format(adjusted_pvalues, scientific = FALSE)
  
)
print(final_table)



### pairwise testing

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Caulobacter")  %>% filter(Treatment=="L"| Treatment=="B")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Caulobacter")  %>% filter(Treatment=="L"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value 

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Caulobacter")  %>% filter(Treatment=="B"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value 

# raw results
raw_pvalues<-c( n1, n2, n3)
taxa <-  "Caulobacter"
pair <- c("L-B", "L-G", "B-G")
#  # 9. Apply the Holm (Holm-Bonferroni) Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  pair = pair,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues
  #fdr_adj_p_nonsci = format(adjusted_pvalues, scientific = FALSE)
  
)
print(final_table)




### 
n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Nitrosocosmicus")  %>% filter(Treatment=="L"| Treatment=="B")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Nitrosocosmicus")  %>% filter(Treatment=="L"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value 

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Nitrosocosmicus")  %>% filter(Treatment=="B"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value 

# raw results
raw_pvalues<-c( n1, n2, n3)
taxa <-  "Nitrosocosmicus"
pair <- c("L-B", "L-G", "B-G")
# p Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  pair = pair,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues,
  fdr_adj_p_nonsci = format(adjusted_pvalues, scientific = FALSE)
  
)
print(final_table)

### rhizobium 
n<-unique(df$Genus)

df1<-df %>% filter(Genus=="Rhizobium")  %>% filter(Treatment=="L"| Treatment=="B")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Rhizobium")  %>% filter(Treatment=="L"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value 

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Rhizobium")  %>% filter(Treatment=="B"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value 

# raw results
raw_pvalues<-c( n1, n2, n3)
taxa <-  "Rhizobium"
pair <- c("L-B", "L-G", "B-G")
# p Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  pair = pair,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues)
  

print(final_table)


## "Polaromonas"
n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Polaromonas")  %>% filter(Treatment=="L"| Treatment=="B")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Polaromonas")  %>% filter(Treatment=="L"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value 

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Polaromonas")  %>% filter(Treatment=="B"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value 

# raw results
raw_pvalues<-c( n1, n2, n3)
taxa <-  "Polaromonas"
pair <- c("L-B", "L-G", "B-G")
#  # 9. Apply the Holm (Holm-Bonferroni) Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  pair = pair,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues)
  
print(final_table)


## "Novosphingobium"
n<-unique(df$Genus)
n
df1<-df %>% filter(Genus=="Novosphingobium")  %>% filter(Treatment=="L"| Treatment=="B")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Novosphingobium")  %>% filter(Treatment=="L"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value 

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Novosphingobium")  %>% filter(Treatment=="B"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value 

# raw results
raw_pvalues<-c( n1, n2, n3)
taxa <-  "Novosphingobium"
pair <- c("L-B", "L-G", "B-G")
#  p Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  pair = pair,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues)

print(final_table)


## Saccharimonadales
n<-unique(df$Genus)
n
df1<-df %>% filter(Genus=="Saccharimonadales")  %>% filter(Treatment=="L"| Treatment=="B")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n1<-m1$p.value

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Saccharimonadales")  %>% filter(Treatment=="L"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n2<-m1$p.value 

n<-unique(df$Genus)
df1<-df %>% filter(Genus=="Saccharimonadales")  %>% filter(Treatment=="B"| Treatment=="G")
m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
n3<-m1$p.value 

# raw results
raw_pvalues<-c( n1, n2, n3)
taxa <-  "Saccharimonadales"
pair <- c("L-B", "L-G", "B-G")
#  p Adjustment

adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
adjusted_pvalues
final_table <- data.frame(
  taxa = taxa,
  pair = pair,
  Raw_P = raw_pvalues,
  fdr_adj_p = adjusted_pvalues)

print(final_table)



########deseq ########
library(phyloseq)
library(DESeq2)
library(ggplot2)

# asv level ##
asvkp<-unique(top_spp$asv)


Workshop_OTU <-otu_table(ps)+1
Workshop_metadat <- sample_data(ps)
Workshop_taxo <- tax_table(ps) 
ps.plusone <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps.plusone



# Focus only on specific genera
ps.subset <- subset_taxa(ps.plusone, asv %in% asvkp)


#  DESeq2 
ds <- phyloseq_to_deseq2(ps.subset, ~ factor(Legume))
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check results
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  


# DEseq2
ds <- phyloseq_to_deseq2(ps.subset, ~ Grass)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  


# DEseq
ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  


### aggregate 2 genus level ####

library(phyloseq)
library(DESeq2)
library(ggplot2)

# 1. Convert phyloseq object to DESeq2 format
# We add 1 to all counts to handle zeros (Laplace smoothing)
# import it phyloseq
target<- c( " g__GWC2-73-18" ,              " g__SZUA-359"    ,       " g__UBA1020"    ,       
                 " g__Caulobacter_487784", " g__Rhizobium_C_501058" ," g__Nitrosocosmicus"   )
ps.subset <- subset_taxa(ps1, Genus %in% target)


# aggregate by genus
df<-as.data.frame(t(otu_table(ps.subset)))
tax<- as.data.frame(tax_table(ps.subset))
tax <- tax %>% select(-Confidence)
df<-cbind(df,tax)
df

# aggregate df
df<-df %>% group_by(Genus) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
df<-as.data.frame(df)
row.names(df)<-df$Genus
df<-df %>% select(-Genus)
df<-df+1
df
# aggregate tax
tax
tax<-tax %>% select(-asv, -Species) %>%
  group_by( Phyla,Class,  Genus) %>%
  summarise()
tax<-as.data.frame(tax)
tax
df
row.names(tax) <- tax$Genus
tax<-as.data.frame(tax[order(tax$Genus),])
as.matrix(tax)

Workshop_OTU <- otu_table(df, taxa_are_rows = TRUE)
Workshop_metadat <- sample_data(ps1)
Workshop_taxo <- tax_table(as.matrix(tax))
ps.subset <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps.subset


#  DESeq2 
ds <- phyloseq_to_deseq2(ps.subset, ~ factor(Legume))
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab

# DEseq2
ds <- phyloseq_to_deseq2(ps.subset, ~ Grass)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  


# DEseq
ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  


# DEseq
ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae*Grass) #~ Site + Group	Controls for "Site" differences while testing for "Group" effects.
diagdds <- DESeq(ds, test="Wald", fitType="local")
resultsNames(diagdds)
resB <- results(diagdds, name="Brassicae_1_vs_0")
res_interact <- results(diagdds, name="Brassicae1.Grass1")
resG<- results(diagdds, name="Grass_1_vs_0")
# check
res_interact
alpha = 0.05
sigtab = res_interact[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  
# check
resB 
sigtab = resB[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  
# check
resG
sigtab = resG[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  


# DEseq
ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae*Legume) #~ Site + Group	Controls for "Site" differences while testing for "Group" effects.
diagdds <- DESeq(ds, test="Wald", fitType="local")
resultsNames(diagdds)
res_interact <- results(diagdds, name="Brassicae1.Legume1")
res_interact
# check
alpha = 0.05
sigtab = res_interact[which(res_interact$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab 



# DEseq
ds <- phyloseq_to_deseq2(ps.subset, ~ Grass*Legume) #~ Site + Group	Controls for "Site" differences while testing for "Group" effects.
diagdds <- DESeq(ds, test="Wald", fitType="local")
resultsNames(diagdds)
res_interact <- results(diagdds, name="Grass1.Legume1")
res_interact
resG <- results(diagdds, name="Grass_1_vs_0")
resL <- results(diagdds, name="Legume_1_vs_0")
# check
alpha = 0.05
sigtab = res_interact[which(res_interact$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab 
# check
resG
sigtab = resG[which(resG$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  




# DEseq
ps.subset1 <- subset_samples(ps.subset, Grass=="1")
ds <- phyloseq_to_deseq2(ps.subset1, ~ Brassicae)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  

# DEseq
ps.subset1 <- subset_samples(ps.subset, Grass=="0")
ds <- phyloseq_to_deseq2(ps.subset1, ~ Brassicae)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  






# DEseq
ps.subset1 <- subset_samples(ps.subset, Brassicae=="1")
ds <- phyloseq_to_deseq2(ps.subset1, ~ Grass)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  

# DEseq
ps.subset1 <- subset_samples(ps.subset, Brassicae=="0")
ds <- phyloseq_to_deseq2(ps.subset1, ~ Grass)
diagdds <- DESeq(ds, test="Wald", fitType="local")
res <- results(diagdds)
# check
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
sigtab  
