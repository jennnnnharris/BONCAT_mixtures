# extract key taxa Figure 5



# clear workspace
rm(list=ls())

## Load required libraries ##

library(tidyverse)
library(vegan)
library(phyloseq)
#library(multcompView)
#library(BiodiversityR)
#ANCOM
#BiocManager::install("ANCOMBC")
#library(ANCOMBC)
#BiocManager::install("microbiome")
#install.packages("microbiome")
library(microbiome)

# set colors

IBM <- c( #IBM colors
  "grey",
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)


#####Import data active #####
## Set the working directory ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")
taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat2.csv",  row.names = 2)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

## order metadata
head(metadat)
metadat$SampleID<-row.names(metadat)
metadat

####filter for just flow cyto samples #
asvs<-asvs[which(metadat$Fraction!="Total"),]
metadat<-metadat[which(metadat$Fraction!="Total"),]
dim(asvs)
dim(metadat)
metadat$Legume   <- factor(metadat$Legume)
metadat$Brassicae   <- factor(metadat$Brassicae)
metadat$Grass   <- factor(metadat$Grass)
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes present", "Legumes absent"))
metadat

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
# remove rare taxa
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


#####extract key taxa

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

commdata = as.data.frame(otu_table(ps1))
# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(dist_matrix ~ Treatment,
                       data = metadat2,
                       comm = commdata,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

# look for key taxa
#plot(cap_result, display = c("sites", "species"))
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
top_spp <- head(species_scores, 20)
top_spp

#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
#svg(filename="active.vectors.svg", height = 6, width = 6)
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
  labs(title = "Top 20 ASV Contributing to Active CAP Variation",
       x = "CAP1", y = "CAP2")
#dev.off()
# plot abundance of these taxa across treatments

# get data
ps1
##add add abundance and taxon infoget
asvkp<-unique(top_spp$asv)


### figure with out soil #######
# make df relative abundance
df<-as.data.frame((otu_table(ps1)))
df<-df/rowSums(df)*100
df<- as.data.frame(t(df))
df$asv<-row.names(df)
# list of taxa
head(df)
# check
tax$asv
active<-df$asv[df$asv %in% tax$asv]
tax[tax$asv %in% active,]


# filter
df <- df[df$asv %in% asvkp, ]
df$asv <-NULL

# list of taxa 
head(df)

# add treatment info
df<-as.data.frame(t(df))
metadat3<-metadat2 %>% select(Treatment, Trt_ID)
df<-cbind(df, metadat3)
df

df<-df %>%
  pivot_longer(cols = c(-Treatment, -Trt_ID) ,
               names_to = "asv",
               values_to = "percent_abundance")
df
# add phyla info
tax<-as.data.frame((tax_table(ps1)))
tax <- tax[tax$asv %in% asvkp, ]
tax<-tax %>% group_by(Phyla, Order, Family, Genus, Species, asv) %>% summarise()
tax<-tax %>% group_by(Phyla, Genus, asv) %>% summarise()
tax
 tax$Blast_ID <- c("Actinomycetes",
                   "Actinomycetes",
                   "Actinomycetes",
                   "Actinomycetes",
                   "Terrabacter sp.",
                   "Microcoleus sp.",
                   "Microcoleus sp.",
                   "Cyanobacterium",
                   "Cyanobacterium",
                   "Cyanobacterium",
                   "Cyanobacterium",
                   "Deinococcus sp.",
                   "Deinococcus sp.",
                   "Deinococcus sp.",
                   "Deinococcus sp.",
                   "Escherichia sp.",
                   "Escherichia sp.",
                   "Escherichia sp.",
                   "Rhizobium sp.",
                   "Rhizobium sp."
 )
tax
df<-left_join(df, tax)
df<-as.data.frame(df)
df

df$percent_abundance <- as.numeric(df$percent_abundance)
df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))




#aggregate by taxa after blasting the sequences
df

df<-df %>% group_by(Phyla, Treatment, Trt_ID, Blast_ID) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
df

df$Treatment
# IBM colors

IBM <- c( #IBM colors
  #"grey",
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
#svg(filename="active.taxa.svg", height = 4, width = 8)
df %>%
  ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
  geom_boxplot(outliers=FALSE)+
  geom_jitter()+
  scale_fill_manual(values= IBM)+
  theme_bw(base_size = 12) +
  facet_grid(~Blast_ID, scales="free", space="free")+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),
        legend.position = "none")+
  labs(title = "A       Active",
       x = "", y = "Percent Abundance")
  

#dev.off()
# 
# #find outliers cutoff
# df1<-df%>% filter(Blast_ID=="Actinomycetes") 
# summary(df1$percent_abundance)
# mean(df1$percent_abundance) + 3*sd(df1$percent_abundance)
# mean(df1$percent_abundance) + 2.5*sd(df1$percent_abundance)
# # Q3 + 1.5 * IQR
# #find outliers cutoff
# df1<-df%>% filter(Blast_ID=="Deinococcus sp.") 
# summary(df1$percent_abundance)
# mean(df1$percent_abundance) + 3*sd(df1$percent_abundance)
# mean(df1$percent_abundance) + 2.5*sd(df1$percent_abundance)
# # Q3 + 1.5 * IQR


# # IBM colors
# # outliers greater than 2.5 sd removed
# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
# #svg(filename="active.taxa.svg", height = 4, width = 8)
# df %>% #filter(percent_abundance<50) %>%
#   ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
#   geom_boxplot(outliers=FALSE)+
#   geom_jitter()+
#   scale_fill_manual(values= IBM)+
#   theme_bw(base_size = 12) +
#   facet_grid(~Blast_ID, scales="free", space="free")+
#   theme(axis.text.x = element_text(angle=60, hjust=1),
#         plot.title = element_text(hjust = 0),
#         legend.position = "none")+
#   labs(title = "A       Active",
#        x = "", y = "Percent Abundance")
#   
#  # dev.off()

########### figure with soil ##############


#  Constrained ordination
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="CTL")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1845 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment!="CTL")
#factor
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction   <- factor(metadat2$Fraction)

# get data
ps1
##add add abundance and taxon infoget
asvkp<-unique(top_spp$asv)

# make df relative abundance
df<-as.data.frame((otu_table(ps1)))
df<-df/rowSums(df)*100
df<- as.data.frame(t(df))
df$asv<-row.names(df)
df <- df[df$asv %in% asvkp, ]
df$asv <-NULL

# add treatment info
df<-as.data.frame(t(df))
metadat3<-metadat2 %>% select(Treatment, Trt_ID)
df<-cbind(df, metadat3)
df

df<-df %>%
  pivot_longer(cols = c(-Treatment, -Trt_ID) ,
               names_to = "asv",
               values_to = "percent_abundance")
df
# add phyla info
tax<-as.data.frame((tax_table(ps1)))
tax <- tax[tax$asv %in% asvkp, ]
tax<-tax %>% group_by(Phyla, Order, Family, Genus, Species, asv) %>% summarise()
tax<-tax %>% group_by(Phyla, Genus, asv) %>% summarise()
tax
tax$Blast_ID <- c("Actinomycetes",
                  "Actinomycetes",
                  "Actinomycetes",
                  "Actinomycetes",
                  "Terrabacter sp.",
                  "Microcoleus sp.",
                  "Microcoleus sp.",
                  "Cyanobacterium",
                  "Cyanobacterium",
                  "Cyanobacterium",
                  "Cyanobacterium",
                  "Deinococcus sp.",
                  "Deinococcus sp.",
                  "Deinococcus sp.",
                  "Deinococcus sp.",
                  "Escherichia sp.",
                  "Escherichia sp.",
                  "Escherichia sp.",
                  "Rhizobium sp.",
                  "Rhizobium sp."
)
tax
df<-left_join(df, tax)
df<-as.data.frame(df)
df

df$percent_abundance <- as.numeric(df$percent_abundance)
df$Treatment <- factor(df$Treatment, levels =c ("Soil","L", "G", "B", "GB", "LB", "LG", "LGB"))




#aggregate by taxa after blasting the sequences
df

df<-df %>% group_by(Phyla, Treatment, Trt_ID, Blast_ID) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
df

# set colors

IBM <- c( #IBM colors
  "grey",
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)
# IBM colors
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
#svg(filename="active.taxa.else.svg", height = 4, width = 9.5)
df %>% filter(Blast_ID!="Actinomycetes") %>%
  ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
  geom_boxplot(outliers=FALSE, alpha=.7)+
  #geom_jitter()+
  scale_fill_manual(values= IBM)+
  theme_bw(base_size = 12) +
  #facet_wrap(~Blast_ID, scales="free", nrow= 1)+
  facet_wrap(~Blast_ID,  nrow= 1)+
  
  
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),
        legend.position = "none")+
  labs(title = "A       Active",
       x = "", y = "Percent Abundance")
dev.off()

# rhizobium
unique(df$Blast_ID)
df %>% filter(Blast_ID=="Rhizobium sp.") %>%
  ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
  geom_boxplot(outliers=FALSE, alpha=.7)+
  #geom_jitter()+
  scale_fill_manual(values= IBM)+
  theme_bw(base_size = 12) +
  #facet_wrap(~Blast_ID, scales="free", nrow= 1)+
  facet_wrap(~Blast_ID,  nrow= 1)+
  
  
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),
        legend.position = "none")+
  labs(title = " Active",
       x = "", y = "Percent Abundance")




# # flip cord
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
svg(filename="active.taxa.tall.svg", height = 9.5, width = 4)
df %>%filter(Blast_ID!="Actinomycetes") %>%
  ggplot(aes(y=Treatment, x=percent_abundance, fill = Treatment))+
  geom_boxplot(outliers=FALSE, alpha=.7)+
  #geom_jitter()+
  scale_fill_manual(values= IBM)+
  theme_bw(base_size = 12) +
  #facet_grid(~Blast_ID, scales="free", space="free")+
  #facet_wrap(~Blast_ID, scales="free", ncol= 1)+
  facet_wrap(~Blast_ID,  ncol= 1)+

  theme(axis.text.x = element_text(angle=0, hjust=1),
        plot.title = element_text(hjust = 0),
        legend.position = "none")+
  labs(title = "A       Active",
       x = "", y = "Percent Abundance")
dev.off()


##### stats DESEQ#####
  library(phyloseq)
  library(DESeq2)
  library(ggplot2)
  
  # asv level ##
  asvkp<-unique(top_spp$asv)
  
  
  #  Constrained ordination
  ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="CTL" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  
  
  Workshop_OTU <-otu_table(ps1)+1
  Workshop_metadat <- sample_data(ps1)
  Workshop_metadat$Grass
  Workshop_metadat$Legume
  Workshop_metadat$Brassicae
  Workshop_metadat$Treatment
  Workshop_taxo <- tax_table(ps1) 
  ps.plusone <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
  ps.plusone
  
  
  
  # Focus only on specific genera
  ps.subset <- subset_taxa(ps.plusone, asv %in% asvkp)
  ps.subset
  
  #  DESeq2 
  ds <- phyloseq_to_deseq2(ps.subset, ~ Legume*Brassicae*Grass)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  resultsNames(diagdds)
  res <- results(diagdds)
  # check results
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  
  
  # effect of L*G * B
  resultsNames(diagdds)
  res <- results(diagdds, name = "Legume1.Brassicae1.Grass1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab # nmo three way interaction
  
  
  
  # effect of L*B
  resultsNames(diagdds)
  res <- results(diagdds, name = "Legume1.Brassicae1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  # effect of L*G
  resultsNames(diagdds)
  res <- results(diagdds, name = "Legume1.Grass1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  # effect of G*B
  resultsNames(diagdds)
  res <- results(diagdds, name = "Brassicae1.Grass1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  # effect of G
  resultsNames(diagdds)
  res <- results(diagdds, name = "Grass_1_vs_0")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  # effect of B
  res <- results(diagdds, name = "Brassicae_1_vs_0")
  # check
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  # L effect
  res <- results(diagdds, name = "Legume_1_vs_0")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
    
  
## Treatment effects
  
  
  #  DESeq2 
  ds <- phyloseq_to_deseq2(ps.subset, ~ Legume)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  resultsNames(diagdds)
  res <- results(diagdds)
  # check results
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  
  #  DESeq2 
  ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  resultsNames(diagdds)
  res <- results(diagdds)
  # check results
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  #  DESeq2 
  ds <- phyloseq_to_deseq2(ps.subset, ~ Grass)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  resultsNames(diagdds)
  res <- results(diagdds)
  # check results
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
    
  ### pairwise contrasts 
  # asv level ##
  asvkp<-unique(top_spp$asv)
  
  #  Constrained ordination
  ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="CTL" & Treatment!="Soil" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  Workshop_OTU <-otu_table(ps1)+1
  Workshop_metadat <- sample_data(ps1)
  Workshop_taxo <- tax_table(ps1) 
  ps.plusone <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
  ps.plusone
  
  
  
  # Focus only on specific genera
  ps.subset <- subset_taxa(ps.plusone, asv %in% asvkp)
  ps.subset
  
  ds <- phyloseq_to_deseq2(ps.subset, ~ Treatment)
  dds <- DESeq(ds)
  resultsNames(diagdds)
  res <- results(diagdds)
 
  design(dds) <- ~ group
  dds <- DESeq(dds)
  resultsNames(dds)
  results(dds, contrast=c("Treatment", "B", "LGB"))
  

  
  
#################Import data total ################

### Clear workspace ###
#rstudioapi::restartSession(clean = TRUE)
rm(list=ls())

# 
# # Load required libraries #
# library(tidyverse)
# library(vegan)
# library(readxl)
# library(lubridate)
# library(phyloseq)
# library(multcompView)
# library(BiodiversityR)

# colors

IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)


#import data#
# Set the working directory 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")
taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat2.csv", row.names = 2)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

## order metadata
head(metadat)
metadat$SampleID<-row.names(metadat)
metadat


# order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
row.names(metadat) <- metadat$SampleID
metadat$N   <- factor(metadat$N)
metadat$Legume   <- factor(metadat$Legume)
metadat$Brassicae   <- factor(metadat$Brassicae)
metadat$Grass   <- factor(metadat$Grass)
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes absent", "Legumes present"))

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


# extract key taxa #

#  Constrained ordination
m<-sample_data(ps)
ps1<-subset_samples(ps, Treatment!="Soil"  & N=="1" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset data

# subset metadata

metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil"  & N=="1" & Treatment!="CTL")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

commdata = as.data.frame(otu_table(ps1))
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

# get site scores
# Extract site scores (WA scores are typically used)
site_scores <- as.data.frame(scores(cap_result, display = "sites"))

# Add your metadata to the scores for mapping (e.g., Treatment groups)
site_scores$Group <- metadat2$Treatment 

# Calculate vector length (distance from origin)
species_scores <- species_scores %>%
  mutate(dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist))

# Select the top 10 species
top_spp <- head(species_scores, 20)
top_spp

# 
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


# plot with sites 
plot(cap_result, display = c("sites", "species"))
  
# plot with sites and colors
ordiplot(cap_result, display= "sites", scaling =1,
         main="", cex.lab =.8)
# par(adj = 0)
#title(main= "A")
par(adj=.5)
points(sc_si, 
       col= IBM[metadat2$Treatment],
       pch= c(16,17)[as.factor(metadat2$Nitrogen_label)],
       lwd=1,cex=.8,
       bg=IBM[metadat2$Treatment])

# arrows
#arrows(0, 0, top_spp[,1], top_spp[,2], length = 0.05, col = "black")

ordiellipse(sc_si, metadat2$Treatment,
            kind = "ehull", conf=0.95, 
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

# 




### plots with arrows

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




# get data
##add add abundance and taxon infoget
asvkp<-unique(top_spp$asv)

### figure without soil ##########

# make df relative abundance
df<-as.data.frame((otu_table(ps1)))
df<-df/rowSums(df)*100
df<- as.data.frame(t(df))
df$asv<-row.names(df)
df <- df[df$asv %in% asvkp, ]
df$asv <-NULL

# add treatment info
df<-as.data.frame(t(df))
metadat3<-metadat2 %>% select(Treatment, Trt_ID)
df<-cbind(df, metadat3)
df

df<-df %>%
  pivot_longer(cols = c(-Treatment, -Trt_ID) ,
               names_to = "asv",
               values_to = "percent_abundance")
df
# add phyla info
tax<-as.data.frame((tax_table(ps1)))
tax <- tax[tax$asv %in% asvkp, ]
tax<-tax %>% group_by(Phyla, Order, Family, Genus, Species, asv) %>% summarise()
tax
tax<-tax %>% group_by(Phyla, Genus, asv) %>% summarise()
tax
tax$Blast_ID <- c("Arthrobacter sp. ",
                  "Limnocylindria",
                  "Limnocylindria",
                  "Saccharimonadota",
                  "Saccharimonadota",
                  "Saccharimonadota",
                  "Saccharimonadota",
                  "Saccharimonadota",
                  "Saccharimonadota",
                  "Saccharimonadota",
                  "Caulobacter sp.",
                  "Caulobacter sp.",
                  "Novosphingobium sp.",
                  "Polaromonas sp.",
                  "Rhizobium sp.",
                  "Rhizobium sp.",
                  "Nitrosocosmicus sp.",
                  "Nitrosocosmicus sp.",
                  "Nitrosocosmicus sp.",
                  "Nitrosocosmicus sp."
                  
)
tax
df<-left_join(df, tax)
df<-as.data.frame(df)
df

df$percent_abundance <- as.numeric(df$percent_abundance)
df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))

df

df<-df %>% group_by(Phyla, Treatment, Trt_ID, Blast_ID) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
df

df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))



#find outliers cutoff
df1<-df%>% filter(Blast_ID=="Rhizobium") 
summary(df1$percent_abundance)
mean(df1$percent_abundance) + 3*sd(df1$percent_abundance)

# remove outlier for easier plotting
  # IBM colors
  setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
  #svg(filename="total.taxa.svg", height = 4, width = 11)
  df %>% filter(percent_abundance<10) %>%
    ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter()+
    scale_fill_manual(values= IBM)+
    theme_bw(base_size = 12) +
    facet_grid(~Blast_ID, scales="free", space="free")+
    theme(axis.text.x = element_text(angle=60, hjust=1),
          plot.title = element_text(hjust = 0),
          legend.position = "none")+
    labs(title = "B       Total",
         x = "", y = "Percent Abundance")
  
#  dev.off()

#####figure with soil ######
  #  Constrained ordination
  m<-sample_data(ps)
  ps1<-subset_samples(ps, N=="1" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  # subset data
  
  # subset metadata
  
  metadat2<-filter(metadat, Fraction=="Total"  & N=="1" & Treatment!="CTL")
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
  
  
  # make df relative abundance
  df<-as.data.frame((otu_table(ps1)))
  df<-df/rowSums(df)*100
  df<- as.data.frame(t(df))
  df$asv<-row.names(df)
  df <- df[df$asv %in% asvkp, ]
  df$asv <-NULL
  
  # add treatment info
  df<-as.data.frame(t(df))
  metadat3<-metadat2 %>% select(Treatment, Trt_ID)
  df<-cbind(df, metadat3)
  df
  
  df<-df %>%
    pivot_longer(cols = c(-Treatment, -Trt_ID) ,
                 names_to = "asv",
                 values_to = "percent_abundance")
  df
  # add phyla info
  tax<-as.data.frame((tax_table(ps1)))
  tax <- tax[tax$asv %in% asvkp, ]
  tax<-tax %>% group_by(Phyla, Order, Family, Genus, Species, asv) %>% summarise()
  tax
  tax<-tax %>% group_by(Phyla, Genus, asv) %>% summarise()
  tax
  tax$Blast_ID <- c("Arthrobacter sp. ",
                    "Limnocylindria",
                    "Limnocylindria",
                    "Saccharimonadota",
                    "Saccharimonadota",
                    "Saccharimonadota",
                    "Saccharimonadota",
                    "Saccharimonadota",
                    "Saccharimonadota",
                    "Saccharimonadota",
                    "Caulobacter sp.",
                    "Caulobacter sp.",
                    "Novosphingobium sp.",
                    "Polaromonas sp.",
                    "Rhizobium sp.",
                    "Rhizobium sp.",
                    "Nitrosocosmicus sp.",
                    "Nitrosocosmicus sp.",
                    "Nitrosocosmicus sp.",
                    "Nitrosocosmicus sp."
                    
  )
  tax
  df<-left_join(df, tax)
  df<-as.data.frame(df)
  df
  
  df$percent_abundance <- as.numeric(df$percent_abundance)
  df$Treatment <- factor(df$Treatment, levels =c ("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
  
  df
  
  df<-df %>% group_by(Phyla, Treatment, Trt_ID, Blast_ID) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
  df
  
  df$Treatment <- factor(df$Treatment, levels =c ("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
  

  
  #find outliers cutoff
  df1<-df%>% filter(Blast_ID=="Rhizobium") 
  summary(df1$percent_abundance)
  mean(df1$percent_abundance) + 3*sd(df1$percent_abundance)
  
  # set colors
  
  IBM <- c( #IBM colors
    "grey",
    "navy", # dark royal blue
    "#648FFF", # french blue
    "#785EF0", # light purple
    "#DC267F", # magenta pink 
    "#FE6100", # bright orange
    "#FFB000", # golden yellow
    "#865338" # medium mocha brown
  )
  
  
  ### order by abundance
  # 1. Calculate column sums
  df 
  total_abundance <- colSums(df1)
  
  # 2. Sort names based on those sums (decreasing = TRUE for most to least)
  ordered_names <- names(sort(total_abundance, decreasing = TRUE))
  
  # 3. Reorder the columns of your data frame
  df_ordered <- df[, ordered_names]
  library(forcats)
  
  # Reorder Blast_ID based on the SUM of percent_abundance across all samples
  df <- df %>%
    mutate(Blast_ID = fct_reorder(Blast_ID, percent_abundance, .fun = sum, .desc = TRUE))
  df$Blast_ID
  # remove outlier for easier plotting
  # IBM colors
  setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
  svg(filename="total.taxa.svg", height = 4, width = 11.5)
  df %>% filter(percent_abundance<10) %>%
    ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
    geom_boxplot(outliers=FALSE, alpha=0.7)+
    #geom_jitter()+
    scale_fill_manual(values= IBM)+
    theme_bw(base_size = 11) +
   # facet_grid(~Blast_ID)+
    #facet_wrap(~Blast_ID, scales="free", nrow= 1)+
    facet_wrap(~Blast_ID,  nrow= 1)+
    theme(axis.text.x = element_text(angle=60, hjust=1),
          plot.title = element_text(hjust = 0),
          legend.position = "none")+
    labs(title = "B       Total",
         x = "", y = "Percent Abundance")
  
    dev.off()
    
    # just rhizobium 
    df %>% filter(percent_abundance<10) %>%
      filter(Blast_ID=="Rhizobium sp.") %>%
      ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
      geom_boxplot(outliers=FALSE, alpha=0.7)+
      #geom_jitter()+
      scale_fill_manual(values= IBM)+
      theme_bw(base_size = 11) +
      facet_wrap(~Blast_ID,  nrow= 1)+
      theme(axis.text.x = element_text(angle=60, hjust=1),
            plot.title = element_text(hjust = 0),
            legend.position = "none")+
      labs(title = "       Total",
           x = "", y = "Percent Abundance")
    
    # just nitrocosmicus
    df %>% filter(percent_abundance<10) %>%
      filter(Blast_ID=="Polaromonas sp.") %>%
      ggplot(aes(x=Treatment, y=percent_abundance, fill = Treatment))+
      geom_boxplot(outliers=FALSE, alpha=0.7)+
      #geom_jitter()+
      scale_fill_manual(values= IBM)+
      theme_bw(base_size = 11) +
      facet_wrap(~Blast_ID,  nrow= 1)+
      theme(axis.text.x = element_text(angle=60, hjust=1),
            plot.title = element_text(hjust = 0),
            legend.position = "none")+
      labs(title = "       Total",
           x = "", y = "Percent Abundance")
    
    
    # cord flip
    setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
    svg(filename="total.taxa.tall.svg", height = 11, width = 4)
    df %>% filter(percent_abundance<10) %>%
      ggplot(aes(x=percent_abundance, y=Treatment, fill = Treatment))+
      geom_boxplot(outliers=FALSE, alpha=0.7)+
      #geom_jitter()+
      scale_fill_manual(values= IBM)+
      theme_bw(base_size = 11) +
      # facet_grid(~Blast_ID)+
      #facet_wrap(~Blast_ID, scales="free", nrow= 1)+
      facet_wrap(~Blast_ID,  ncol= 1)+
      theme(axis.text.x = element_text(angle=60, hjust=1),
            plot.title = element_text(hjust = 0),
            legend.position = "none")+
      labs(title = "B       Total",
           x = "", y = "Percent Abundance")
    
    dev.off()
  
# stats ########
  library(phyloseq)
  library(DESeq2)
  library(ggplot2)
  
  # asv level ##
  asvkp<-unique(top_spp$asv)
  
  #  filter for treatments
  #  Constrained ordination
  
  ps1<-subset_samples(ps,  N=="1" )
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  Workshop_OTU <-otu_table(ps1)+1
  Workshop_metadat <- sample_data(ps1)
  
  # Workshop_metadat$Grass
  # Workshop_metadat$Legume
  # Workshop_metadat$Brassicae
  # Workshop_metadat$Treatment
  Workshop_taxo <- tax_table(ps1) 
  ps.plusone <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
  ps.plusone
  
  # Focus only on specific genera
  ps.subset <- subset_taxa(ps.plusone, asv %in% asvkp)
  ps.subset
  
  #  DESeq2 
  ds <- phyloseq_to_deseq2(ps.subset, ~ Legume*Brassicae*Grass)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  resultsNames(diagdds)
  res <- results(diagdds)
  # check results
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  
  
  
  # effect of L*G * B
  resultsNames(diagdds)
  res <- results(diagdds, name = "Legume1.Brassicae1.Grass1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab # nmo three way interaction
  
  
  
  # effect of L*B
  resultsNames(diagdds)
  res <- results(diagdds, name = "Legume1.Brassicae1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  # effect of L*G
  resultsNames(diagdds)
  res <- results(diagdds, name = "Legume1.Grass1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  # effect of G*B
  resultsNames(diagdds)
  res <- results(diagdds, name = "Brassicae1.Grass1")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  # effect of G
  resultsNames(diagdds)
  res <- results(diagdds, name = "Grass_1_vs_0")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  # effect of B
  res <- results(diagdds, name = "Brassicae_1_vs_0")
  # check
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  # L effect
  res <- results(diagdds, name = "Legume_1_vs_0")
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.subset)[rownames(sigtab), ], "matrix"))  
  sigtab 
  
  
  
  
  
  
  
  
  
  
# treatment effects  
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
  
#### pairwise effect