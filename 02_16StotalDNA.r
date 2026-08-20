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
#library(readxl)
library(lubridate)
library(phyloseq)
library(multcompView)
library(BiodiversityR)

# colors

mycols <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)

#"#8cadff"
bw <- c( #IBM colors
  "grey", # grey
  "black" # teal
)

# set shapes
myshapes <- c(1, 12,15 ,21, 22, 23 , 24)
myshapes2 <- c(21 , 12, 24,1, 15 , 22, 23 )


#import data#
# Set the working directory 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")

taxon <- read.csv("16S_taxonomy.csv", header=T)
asvs <- read.table("16S_feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("16S_metadata.csv", header = T, row.names = 1)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

# order metadata
metadat
metadat<-as.data.frame(metadat[order(row.names(metadat)),])
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat$N   <- factor(metadat$N)
metadat$Legume   <- factor(metadat$Legume)
metadat$Brassicae   <- factor(metadat$Brassicae)
metadat$Grass   <- factor(metadat$Grass)
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes absent", "Legumes present"))
head(metadat)
asvs[1:5,1:5]

#T_DNA_23_S153 has really few reads so I am omitting it.
asvs<-asvs[which(row.names(asvs)!= "T_DNA_23_S153"),]

# get avg number of reads in seq run 
mean(rowSums(asvs))


# select only total dna
asvs <-asvs[which(metadat$Fraction=="Total" ),] 
metadat <- metadat %>% filter(Fraction=="Total" )


### check rarefaction curves to see if sequencing depth is sufficient ###
#get min number of reads in a sample
#min.s<-min(rowSums(asvs))
# observe number of species
#S <- specnumber(asvs)

# check library saturation with rare curve
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/rarefaction")
#svg("total_rarecurve.svg",  width=10, height=10)
#rarecurve(asvs, step = 1000, col = "blue", xlab = "Sample Size", ylab = "Species Richness")
#abline(v = min.s, lty = 2)
#dev.off()

### plot rarefied species to observed species
S <- specnumber(asvs) # observed number of species
min.s<-min(rowSums(asvs))
Srare <- rarefy(asvs, min.s) # rarefied number of species
#svg("total_rarefiedvs_observed.svg",  width=10, height=10)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
#dev.off()


### Rarefy to obtain even numbers of reads by sample ###
min.s<-min(rowSums(asvs))
set.seed(336)
asvs<-rrarefy(asvs, min.s)

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
# 220K taxa and 86 samples when mitochondria and chloroplasts removed. 

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
       y= "Shannon diversity")
  #annotate("text", x=1.5, y=8, label="*")
p1
p2<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=Observed))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  labs(title = "B",
       x="",
       y= "ASV richness")

p2
p3<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(Nitrogen_label), y=evenness))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  labs(title = "C",
       x="",
       y= "Pilou's Eveness")
p3
require(gridExtra)
windows(8, 3.5)
grid.arrange(p1, p2, p3, ncol=3)


# diversity by treatment #

rich<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") #filter
rich$Nitrogen_label<-as.factor(rich$Nitrogen_label)

# shannon diversity

p1<-rich%>%  
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "D",
       x="",
       y= "ASV richness")+
  scale_shape_manual(values = c(17, 16))

p1

#  OBSERVED ASV
p2<-rich%>%  
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "E",
       x="",
       y= "ASV Richness")+
  scale_shape_manual(values = c(17, 16))+
  ylim(2500, 5500)
p2


p3<-rich%>%  
  ggplot(aes(x=Treatment, y=evenness,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  geom_jitter(aes(shape=Nitrogen_label), size=1.5, width=.1)+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust = 0),legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "F",
       x="",
       y= "Pilou's Evenness")+
  scale_shape_manual(values = c(17, 16))
p3


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



####### 3. Filtering rare taxa ##################

# remove true singletons 
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
ps
# 187k 

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


####### 4. CAP model L * G *B Nitrogen + #######

# subset data
# nitrogen + only
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
cap_result <- capscale(dist_matrix ~ Legume*Grass*Brassicae+block,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

anova.cca(cap_result, by="terms")

####### 4. CAP model L * G *B Nitrogen - #######
# subset data
# nitrogen + only
ps1 <-subset_samples(ps, Treatment !="Soil" & N=="0" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil"  & N=="0")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
head(metadat2)
# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Legume*Grass*Brassicae+block,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

anova.cca(cap_result, by="terms")


####### 5. CAP model Treatment* Nitrogen ##################
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
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal") 
#svg("cap.total.supplement.svg", width = 7.5, height = 4)

par(mfrow=c(1,2))
par(adj=0)
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         main="A", cex.lab =1,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj=.5)
points(sc_si, 
       col= mycols[metadat2$Treatment],
       pch= c(16,17)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=1,
       bg=mycols[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=F, 
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 30,
            cex=1)

###### CAP nitrogen +  and - comparision
par(adj=0)
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         cex.lab = 1, main = "B",
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj=.5)
points(sc_si, 
       col= bw[as.factor(metadat2$Nitrogen_label)],
       pch= c(21,24)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=1,
       bg=bw[as.factor(metadat2$Nitrogen_label)])
ordiellipse(sc_si, metadat2$Nitrogen_label,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            col= bw,
            alpha = 30,
            cex=1)
legend("bottomleft", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
       fill= mycols,
       cex=.8)

legend("bottom", legend=c("Nitrogen +", "Nitrogen -"  ),
      pch=c(16,17 ),
       cex=.8)


dev.off()

####### 6. N+ only and N - only CAP plots #####
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
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
#svg("cap.total2.svg", width = 8, height = 5)

par(mfrow=c(1,2))
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
          cex.lab = 1.2,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj = 0)
#title(main= "C")
par(adj=.5)
points(sc_si, 
       col= mycols[metadat2$Treatment],
       pch= c(21)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=1,
       bg=mycols[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=F, 
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 50,
            cex=1)


####### N- only  

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
         cex.lab = 1.2,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj = 0)
#title(main= "D")
par(adj=.5)
points(sc_si, 
       col= mycols[metadat2$Treatment],
       pch= c(24)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=1,
       bg=mycols[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=F, 
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 50,
            cex=1)


dev.off()

####### 7. pairwise tables ###########

# nitrogen +
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
#View the raw pairwise results
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
  geom_text(aes(label = p_value), size = 2.5) +
  scale_fill_gradient(
    low = "#648FFF",
    high = "white",
    na.value = "grey90",
    limits = c(0, 0.1),
    name = "P-value"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_recol_covercrop

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
pdf("Pval_withN.pdf",  width=4, height=3 )
p_recol_covercrop
dev.off()


####### N- only  


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
  geom_text(aes(label = p_value), size = 2.5) +
  scale_fill_gradient(
    low = "#648FFF",
    high = "white",
    na.value = "grey90",
    limits = c(0, 0.1),
    name = "P-value"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_recol_covercrop

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPtotal")
pdf("Pval_withoutN.pdf", width=4, height=3)
p_recol_covercrop
dev.off()




