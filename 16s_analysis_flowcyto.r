# 16S analysis
# total
# beta diversity
# BONCAT mixtures
# Jennifer Harris
# Jan 7 2025
# Last Updated: May 2025

#R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

### Initial Setup ###

### Clear workspace ###

#rstudioapi::restartSession(clean = TRUE)
rm(list=ls())


## Load required libraries ##

#basic
library(tidyverse)
library(vegan)
library(readxl)
library(lubridate)
library(phyloseq)
library(multcompView)
library(BiodiversityR)
#ANCOM
#BiocManager::install("ANCOMBC")
library(ANCOMBC)
#BiocManager::install("microbiome")
#install.packages("microbiome")
library(microbiome)

# set colors
mycols <- c( #IBM colors
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

bw <- c( #IBM colors
  "grey", # grey
  "black" # black
)
# set shapes
myshapes <- c(1, 12,15 ,21, 22, 23 , 24)
myshapes2 <- c(21 , 12, 24,1, 15 , 22, 23 )

#####Import data#####
## Set the working directory ###
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")


taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat.csv", header = T)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

## order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
row.names(metadat) <- metadat$SampleID
metadat

####filter for just flow cyto samples #
asvs<-asvs[which(metadat$Fraction!="Total"),]
dim(asvs)

metadat<-metadat[which(metadat$Fraction!="Total"),]
dim(metadat)
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes present", "Legumes absent"))



#make taxon matrix row names OTUs
#taxon[1:5,1:5]
row.names(taxon) <- taxon$asv

#get min number of reads in a sample
min.s<-min(rowSums(asvs))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs, min.s)

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs.r), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 110K taxa when rarefied 


###remove plant contamination  ##
# Select unassigned Asvs that the only in the the roots and nodules
ps<-subset_taxa(ps, Class!="c__Chloroplast" )
ps<-subset_taxa(ps, Class!=" c__Chloroplast" )
ps<-subset_taxa(ps, Family!= " f__Mitochondria" )

ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 110K taxa 

####DIVERSITY ####

#overall 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1"))

# Data wrangling fo rdiversity of active microbes in each fraction
metadat.t <- sample_data(ps)
rich<-cbind(rich, sample_data(ps))
rich<-as.data.frame(rich)

# summary
rich %>% group_by(Fraction) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(n_species) %>% summarise(mean(Observed), sd(Observed))

#
rich$Treatment   <- factor(rich$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
rich$Fraction   <- factor(rich$Fraction, levels= c("Active", "Inactive"))




#BCAT_73_S35 is kind of a weird outlier
#rich<-rich[which(rich$SampleID!="BCAT_73_S35"),]


#plots
# Does diversity increase with number of species?

p1<-rich%>%  filter(Fraction=="Active") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
 # scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  #geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "C",
                  x="",
                  y= "Active Shannon Diversity")
 #scale_shape_discrete() 
p1

p2<-rich%>%  filter(Fraction=="Active") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  #geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "D",
       x="",
       y= "Active ASV richness")
  #scale_shape_discrete() 
p2

p3<-rich%>%  filter(Fraction=="Active") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Chao1,  fill=Treatment))+
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values = mycols)+
  #geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("active Chao1 species richness")+
  xlab("")
  #scale_shape_discrete() 
p3

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_diversity")
svg(file="diversity.fc.svg",width = 10, height=4)

require(gridExtra)
#windows(6,8)
grid.arrange(p1, p2, ncol=2)
dev.off()



##DIVERSITY STATS######

# active
# remove soil
rich <- rich %>% filter(Fraction=="Active") %>% filter(Treatment!="Soil")
rich$Treatment   <- factor(rich$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# 1. Shannon ANOVA 
anova1<- aov(Shannon~ Treatment, data = rich)
summary(anova1)
# 2. Run the Tukey HSD test on the ANOVA model
tukey_output <- TukeyHSD(anova1)

# The 'tukey_output' shows the pairwise comparisons and p-values
print(tukey_output)

# 3. Extract the p-values for the factor of interest
p_values <- tukey_output$Treatment[, 4]

# 4. Generate the grouping letters using multcompLetters()
# The 'multcompLetters()' function takes a named vector of p-values
cld <- multcompLetters(p_values)

# 5. Print the results
print(cld)



#observed ANOVA
m1<-lm(Observed ~ Treatment,  data = rich)
summary(m1)

#chao1
m1<-lm(Chao1 ~ Treatment,  data = rich)
summary(m1)

# Analysis of variance 
anova1<- aov(Chao1 ~ Treatment, data = rich)
summary(anova1)
tukey.a1 <- TukeyHSD(anova1)
print(a1) # all difference except LG-LB and LGB-LG
plot(anova1) #homoscedasticity looks fine

# n species
m1<-lm(Shannon ~ rich$n_species,  data = rich)
summary(m1)
m1<-lm(Observed ~ rich$n_species,  data = rich)
summary(m1)
m1<-lm(Chao1 ~ rich$n_species,  data = rich)
summary(m1)





# remove rare taxa ############
  ps<-prune_taxa(taxa_sums(ps) > 0, ps)
  ps
  
  
  # remove true singletons 
  ps<-prune_taxa(taxa_sums(ps) > 1, ps)
  ps
  # 
  
  #remove asvs with a mean of less than 5
  mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
  keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
  ps<-prune_taxa(keep, ps)
  ps
  
    
#####PCOA   ########
  
# subset data  
  ps1 <-subset_samples(ps, Treatment !="Soil" & Treatment!="CTL" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
   metadat2<-filter(metadat, Treatment!="Soil" & Treatment!="CTL")
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
  metadat2$Fraction   <- factor(metadat2$Fraction)
  
  
  # Calculate Bray-Curtis distance between samples
  otus.bray<-vegdist(otu_table(ps1), method = "bray")
  # Perform PCoA analysis of BC distances #
  otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)
  # Store coordinates for first two axes in new variable #
  otus.p <- otus.pcoa$points[,1:2]
  otus.p3 <- otus.pcoa$points[,3:4]
  colnames(otus.p) <- c("PC1", "PC2")
  colnames(otus.p3) <- c("PC3", "PC4")
  
  df.pcoa <- cbind(sample_data(ps1), otus.p)
  df.pcoa<-cbind(df.pcoa, otus.p3)
  df.pcoa
  # Calculate % variance explained by each axis #
  otus.eig<-otus.pcoa$eig
  perc.exp<-otus.eig/(sum(otus.eig))*100
  pe1<-perc.exp[1]
  pe2<-perc.exp[2]
  pe3<-perc.exp[3]
  pe4<-perc.exp[4]
  
  #calculate total variance explained by each principal component
  perc.exp<-otus.eig/(sum(otus.eig))*100
  #scree plot 
  otus.pcoa$eig
  plot(perc.exp[1:8],
       ylab = "percent varience explained",
       xlab = "PC")


  ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCOA ",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
           ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))+
  par(adj = 0)
  title(main= "B")
  par(adj=.5)
  points(otus.p, 
         col= bw[as.factor(metadat2$Fraction)],
         pch= c(22,24)[as.factor(metadat2$Fraction)],
         lwd=1,cex=1.5,
         bg=bw[as.factor(metadat2$Fraction)],)
  ordiellipse(otus.pcoa, as.factor(metadat2$Fraction),  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              col= bw,
              alpha = 30,
              cex=1.5)
  # base r plot
  
  
 
#######PCOA plot active#########
  ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  # 1845 taxa
  # subset metadata
  metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
  
  #factor
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
  metadat2$Fraction   <- factor(metadat2$Fraction)


# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
otus.p3 <- otus.pcoa$points[,3:4]
colnames(otus.p) <- c("PC1", "PC2")
colnames(otus.p3) <- c("PC3", "PC4")

df.pcoa <- cbind(sample_data(ps1), otus.p)
df.pcoa<-cbind(df.pcoa, otus.p3)
df.pcoa$Treatment   <- factor(df.pcoa$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))

df.pcoa
# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]
pe3<-perc.exp[3]
pe4<-perc.exp[4]

#calculate total variance explained by each principal component
perc.exp<-otus.eig/(sum(otus.eig))*100
#scree plot 
otus.pcoa$eig
plot(perc.exp[1:8],
     ylab = "percent varience explained",
     xlab = "PC")
#faceted plot ##

#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols, name="treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% var. explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "PCoA Active ",
       subtitle = "A")+
  coord_fixed()+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~Legume, ncol=3)



#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC3, y = PC4, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols, name="treatment") +
  labs(x = paste("PCoA3 (",round(pe3,2),"% var. explained)"), y = paste("PCoA4 (",round(pe4,2),"% variance explained)"),
       title = "PCoA Active ",
       subtitle = "A")+
  coord_fixed()+
  stat_ellipse(aes(group=Treatment), linetype=2)
  facet_wrap(~Legume, ncol=3)


############CAP ############  
# inactive verse active 
  # subset data  
  ps1 <-subset_samples(ps, Treatment !="Soil" & Treatment!="CTL" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  metadat2<-filter(metadat, Treatment!="Soil" & Treatment!="CTL")
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
  metadat2$Fraction   <- factor(metadat2$Fraction)
  
  
  # 1. Calculate the distance matrix (e.g., Bray-Curtis)
  dist_matrix<-vegdist(otu_table(ps1), method = "bray")
  
  # 2. Run the CAP (db-RDA) analysis
  # Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
  cap_result <- capscale(dist_matrix ~ Treatment*Fraction*Block,
                         data = metadat2,
                         add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA
  
  # 3. Permutation test for significance of constraints
  anova_cap <- anova(cap_result, permutations = 999, by = "term")
  anova_cap

  
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
  
  #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
  #svg("cap.fraction.svg", width = 6 , height = 6)
# split legume present and absence in base r
  
  #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
  #svg("cap.total.fraction.svg", width = 8 , height = 4.5)

  par(mfrow = c(1, 2)) # 1 row, 2 column
 # windows(4,4)
  par(cex.lab = 1.1) # make all fonts in graphs little bigger
  ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
           main="CAP ", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "A")
  par(adj=.5)
  points(sc_si, 
         col= "black",
         pch= c(22,24)[as.factor(metadat2$Fraction)],
         lwd=1,cex=1,
         bg=bw[metadat2$Fraction])
  ordiellipse(sc_si, metadat2$Fraction,  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              col= bw,
              alpha = 30,
              cex=1.2)
  
  legend("bottomleft", legend=c("Active", "Inactive"  ),
         pch=c(22,24 ),
         cex=1,
         title = "",     bty = "n")
  
  #dev.off()
  
  
  
##### CAP active ##################
#  Constrained ordination
  ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  # 1845 taxa
  # subset metadata
  metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
  #factor
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
  metadat2$Fraction   <- factor(metadat2$Fraction)
  
  # 1. Calculate the distance matrix (e.g., Bray-Curtis)
  dist_matrix<-vegdist(otu_table(ps1), method = "bray")
  
  # 2. Run the CAP (db-RDA) analysis
  # Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
  cap_result <- capscale(dist_matrix ~ Treatment*Block,
                         data = metadat2,
                         add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA
  
  anova.cca(cap_result, by="terms")
  
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
  adjusted_pvalues <- p.adjust(raw_pvalues, method = "bonferroni")
  adjusted_pvalues1 <- p.adjust(raw_pvalues, method = "fdr")
  
  #?p.adjust
  # 10. Combine the results for final interpretation
  final_table <- data.frame(
    Pair = rownames(pairwise_results),
    Pseudo_F = pairwise_results[, "F"],
    Raw_P = raw_pvalues,
    fdr_adj_p = adjusted_pvalues1,
    bonferroni_Adj_P = adjusted_pvalues
  )
  
  # 11. Print the final results table
  print(final_table)
  
  
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
  
  #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
  #svg("cap.active.svg", width = 6 , height = 6)
  #windows(6,6)
  par(cex.lab = 1.1) # make all fonts in graphs little bigger
  ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
           main="CAP ", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "B")
  par(adj=.5)
  points(sc_si, 
         col= IBM[metadat2$Treatment],
         pch= 22,
         lwd=1,cex=1,
         bg=IBM[metadat2$Treatment])
  ordiellipse(sc_si, metadat2$Treatment,  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              col= IBM,
              alpha = 40,
              cex=1)
  legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
         fill= IBM,
         cex=1,
         bty = "n")
  
  
  #dev.off()
  
 ############ addtional cap model L +b +G ###  ##########
  
  ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  # 1845 taxa
  # subset metadata
  metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
  #factor
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
  metadat2$Fraction   <- factor(metadat2$Fraction)
  
  # 1. Calculate the distance matrix (e.g., Bray-Curtis)
  dist_matrix<-vegdist(otu_table(ps1), method = "bray")
  
  # 2. Run the CAP (db-RDA) analysis
  # Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
  cap_result <- capscale(dist_matrix ~Brassicae*Grass*Legume,
                         data = metadat2,
                         add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA
  
  anova.cca(cap_result, by="terms")
  
# split legume present and absence in base r
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
svg("cap.facet.svg", width = 8 , height = 4.5)
#windows(8,4.5)  

    par(mfrow = c(1, 2)) # 1 row, 2 columns  
  # legume present 
legume_cols<- c( #IBM colors
    "navy", # dark royal blue L
    "white", # french blue G
    "white", # light purple  B 
    "white", # magenta pink GB
    "#FE6100", # bright orange LB
    "#FFB000", # golden yellow LG 
    "#865338" # medium mocha brown LGB
  )
  par(cex.lab = 1.1) # make all fonts in graphs little bigger
  ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
           main="Legumes present ", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "C")
  par(adj=.5)
  points(sc_si, 
         col= legume_cols[metadat2$Treatment],
         pch= 22,
         lwd=1,cex=1,
         bg=legume_cols[metadat2$Treatment])
  ordiellipse(sc_si, metadat2$Treatment,  
              kind = "ehull", conf=0.95, label=F, 
              draw = "polygon",
              border = 0,
              col= legume_cols,
              alpha = 40,
              cex=1)
  #legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
  #       fill= IBM,
  #       cex=1,
  #       title = "",
  #       bty = "n")
  
  # legume absent
no_legumes <- c( #IBM colors
    "white", # dark royal blue L 
    "#648FFF", # french blue G
    "#785EF0", # light purple B
    "#DC267F", # magenta pink  GB
    "white", # bright orange LB
    "white", # golden yellow LG
    "white" # medium mocha brown LGB
  )
  par(cex.lab = 1.1) # make all fonts in graphs little bigger
  ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
           main="Legumes absent ", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "D")
  par(adj=.5)
  points(sc_si, 
         col= no_legumes[metadat2$Treatment],
         pch= 22,
         lwd=1,cex=1,
         bg=no_legumes[metadat2$Treatment])
  ordiellipse(sc_si, metadat2$Treatment,  
              kind = "ehull", conf=0.95, label=F, 
              draw = "polygon",
              border = 0,
              col= no_legumes,
              alpha = 40,
              cex=1)

  
dev.off()  
  
  
# factor GGplot  
  p1.cap <- ordinate(ps1, method='CAP',distance='bray',formula=~Treatment)      
  # cap plot total
  p1 <-plot_ordination(ps1,  p1.cap, color="Treatment")+
    theme_bw()+
    geom_point( size=2.5)+
    #stat_ellipse(aes(group=Treatment), linetype=1)+
    theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
          legend.position="left")+
    scale_color_manual(values =  mycols, name="composition")+
    ggtitle("model: composition * N")+
    facet_grid(~mixture)
  p1

    
p2<-  plot_ordination(ps1,  p1.cap, color="Treatment")+
    theme_bw()+
    geom_point( size=2.5)+
    #stat_ellipse(aes(group=Treatment), linetype=1)+
    theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
          legend.position="none")+
    scale_color_manual(values =  mycols, name="composition")+
    ggtitle("C")+
    facet_grid(~Legume_label)
  

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
#svg("cap.active.svg", width = 6 , height = 6)
windows(6,4)
p2
dev.off()

###PCOA STATS: BETA DISPERSION#####
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

#########PERMANOVA active ########### 
#between trts
# Constrained ordination
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1845 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")

#factor
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction   <- factor(metadat2$Fraction)


#get asvs table
asvs.clean<-otu_table(ps1)

# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
#
asvs.perm<- adonis2(asvs.clean ~ Treatment, data = metadat2, permutations = 999, method="bray")
asvs.perm

asvs.perm<- adonis2(asvs.clean ~ (Grass + Brassicae + Legume)^2, data = metadat2, permutations = 999, method="bray")
asvs.perm
adonis2(formula = asvs.clean ~ (Grass + Brassicae + Legume)^3, data = metadat2, permutations = 999, method = "bray", by="term")

#### pairwise adonis
library(pairwiseAdonis)

bray_dist <- vegdist(otu_table(ps1), method = "bray")
# Run the pairwise PERMANOVA
pairwise_results <- pairwise.adonis2(bray_dist ~ Treatment, # Use the distance matrix directly
                                     data = metadat2,
                                     permutations = 999,
                                     method = "bray",
                                     p.adjust.m = "bonferroni")
pairwise_results



########## PERMANOVA Fraction############
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


####PERCENT abundance figure by TREATment: ####
#get data

# Transform ASV counts to proportions / relative abundance
total.prop <- transform_sample_counts(subset_samples(total),function(ASV) ASV/sum(ASV))

# Subset by family
total.prop.order <- phyloseq::tax_glom(total.prop, "Order") # 196 taxa

# Take x number of taxa per taxonomic rank based on relative abundance
# Taxa > n will be added to as "other" label
install.packages("fantaxtic")
total.prop.family.top10 <- fantaxtic::nested_top_taxa(physeq_obj = total.prop.family, n=10, relative = TRUE, discard_other = FALSE, other_label = "Other")

# Melt data frame with phyloseq function for plotting
total.family.top10 <- psmelt(p1.ccmAim3.nod.prop.family.top10)

# Reorder levels to put other at the end - otherwise taxa are in alphabetical order
p1.family.nod.top10$Family <- forcats::fct_relevel(as.factor(p1.family.nod.top10$Family), "Other", after = Inf)

# Visualize
ggplot(p1.family.nod.top10, aes(x=Sample, y=Abundance, fill=Family))+
  geom_bar(stat = "identity")+
  facet_wrap(~CommonCoverCropCode, scales="free_x",
             labeller=as_labeller(c(BC='Bush clover', YSC='Yellow sweet clover',
                                    FP='Field pea', CC='Crimson clover',
                                    WC='White clover', CP='Cowpea')))+
  labs(y="Relative abundance", title="Nodule")+
  theme_bw()+
  theme(legend.position='right', axis.title.x=element_blank(),
        axis.text.x=element_blank(), text=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=16), legend.text=element_text
        (size=8))+
  scale_fill_manual(values=c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", 
                             "#526A83FF", "#625377FF", "#68855CFF", "#9C9C5EFF", 
                             "#A06177FF", "#8C785DFF", "#467378FF"))


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


ps1 <-subset_samples(ps, Fraction=="Inactive" & Treatment!="Soil" )
ps1 <-subset_samples(ps, Treatment!="Soil" )

ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1845 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Inactive" & Treatment!="Soil")

#factor
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction   <- factor(metadat2$Fraction)


# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
otus.p3 <- otus.pcoa$points[,3:4]
colnames(otus.p) <- c("PC1", "PC2")
colnames(otus.p3) <- c("PC3", "PC4")

df.pcoa <- cbind(sample_data(ps1), otus.p)
df.pcoa<-cbind(df.pcoa, otus.p3)
df.pcoa$Treatment   <- factor(df.pcoa$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))

df.pcoa
# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]
pe3<-perc.exp[3]
pe4<-perc.exp[4]

#calculate total variance explained by each principal component
perc.exp<-otus.eig/(sum(otus.eig))*100
#scree plot 
otus.pcoa$eig
plot(perc.exp[1:8],
     ylab = "percent varience explained",
     xlab = "PC")

#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols, name="treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% var. explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "PCoA INActive ",
       subtitle = "A")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~mixture, ncol=3)



#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC3, y = PC4, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols, name="treatment") +
  labs(x = paste("PCoA3 (",round(pe3,2),"% var. explained)"), y = paste("PCoA4 (",round(pe4,2),"% variance explained)"),
       title = "PCoA INActive ",
       subtitle = "A")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~Legume, ncol=3)

dev.off()




##### ANCOM###############
# ancom is run on rarefied filtered data

ps1<-subset_samples(ps , Treatment!="Soil")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)

metadat2 <- as.data.frame(sample_data(ps1))
taxon<-as.data.frame(tax_table(ps1))
df<-as.data.frame(t(otu_table(ps1)))

### skip to data import if you already ran ancom
sample_data(ps1)

out1 = ancombc(data = ps1, assay_name = "counts", 
               tax_level = "asv",  
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = .10, lib_cut = 0, 
               group = "Fraction", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.001, global = FALSE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
sample_data(ps1)


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
#log fold change
tab_lfc = res$lfc
head(tab_lfc)
dim(tab_lfc)
col_name = c("asv", "LFC_Intercept", "LFC_FractionViable_Cell")
colnames(tab_lfc) = col_name
head(tab_lfc)
write_delim(as.data.frame(tab_lfc), file = "ancom_Log_fold_change.txt", delim = " ")

# standard error
tab_se = res$se
col_name = c("asv", "se_Intercept", "se_viable")
colnames(tab_se) = col_name
head(tab_se)
tab_se<-as.data.frame(tab_se)
write.table(tab_se, file = "ancom_SE.txt")
tab_se<-read.table("ancom_SE.txt", header = TRUE)
head(tab_se)

#test statistcs W maybe it's willcoxin?
tab_w = res$W
col_name = c("asv", "W_Intercept", "W_viable")
colnames(tab_w) = col_name
head(tab_w)
write.table(as.data.frame(tab_se), file = "ancom_SE.txt")

# P-values from the Primary Result
tab_p = res$p_val
col_name = c("asv", "p_Intercept", "p_viable")
colnames(tab_p) = col_name
head(tab_p)
write.table(as.data.frame(tab_p), file = "ancom_pval.txt")


#Adjusted p-values from the Primary Result"
tab_q = res$q
head(tab_q)
col_name = c("asv", "adj_p_Intercept", "adj_p_viable")
colnames(tab_q) = col_name
head(tab_q)
write.table(as.data.frame(tab_se), file = "ancom_adjpval.txt")

# yes or no is a taxa differentially abundant
tab_diff = res$diff_abn
col_name = c("asv", "DA_Intercept", "DA_Fraction_Viable_cells_Active")
colnames(tab_diff) = col_name
head(tab_diff)
write_delim(as.data.frame(tab_diff), file = "ancom_DA.txt", delim = " ")


# through all togetha nd remove anything with DNA
tab <-tab_lfc %>%
  left_join(., tab_se) %>%
  left_join(., tab_w ) %>%
  left_join(., tab_p) %>%
  left_join(., tab_q) %>%
  left_join(., tab_diff)
write.table(as.data.frame(tab), file = "ancom_table.txt")
head(tab)
######ANCOM import df######
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")

tab<-read.table("ancom_table.txt", header = TRUE)
head(tab)

##add add abundance and taxon info
df<-as.data.frame((otu_table(ps1)))
df$asv<-row.names(df)
df<-left_join(df, tab)
taxon<-as.data.frame(tax_table(ps1))
df<-left_join(taxon, df)
df

### summarise the abundance in active and Viable
df$rhizo.Viable.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
t<-df %>% select(contains("R.SYBR"))
sd_Viable<- apply(t, 1, sd, na.rm=TRUE)
sd_Viable
df$sd_Viable <- sd_Viable

df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
t<-df %>% select(contains("R.POS"))
sd_active<- apply(t, 1, sd, na.rm=TRUE)
sd_active
df$sd_active <- sd_active
head(df)




