# 16S analysis
# total
# beta diversity
# BONCAT mixtures
# Jennifer Harris
# Jan 7 2025
# Last Updated: May 2025

# laptop R version 4.2.3 (2023-03-15 ucrt) -- "Shortstop Beagle"

### Initial Setup ###

### Clear workspace ###

#rstudioapi::restartSession(clean = TRUE)
rm(list=ls())



## Load required libraries ##

library(tidyverse)
library(vegan)
#library(readxl)
#library(lubridate)
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
metadat<-read.csv("metadat2.csv", header = T, row.names = 1)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

## order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
row.names(metadat) <- metadat$SampleID
metadat

####filter for just flow cyto samples #
asvs<-asvs[which(metadat$Fraction!="Total"),]
metadat<-metadat[which(metadat$Fraction!="Total"),]
dim(metadat)
metadat$Legume   <- factor(metadat$Legume)
metadat$Brassicae   <- factor(metadat$Brassicae)
metadat$Grass   <- factor(metadat$Grass)

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

# caculate pilou's evenness where formula J= H/ln(S), 
rich$evenness = rich$Shannon/log(rich$Observed)



#plots

p1<-rich%>%  filter(Fraction=="Active") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_fill_manual(values = IBM)+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "A",
                  x="",
                  y= "Active Shannon Diversity")
p1

p2<-rich%>%  filter(Fraction=="Active") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_fill_manual(values=IBM) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "B",
       x="",
       y= "Active ASV richness")
p2

p3<-rich%>%  filter(Fraction=="Active") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=evenness,  fill=Treatment))+
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "C",
       x="",
       y= "Active Pilou's Evenness")
  
p3

# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement/diversity_active_inact")
# svg("active_diversity.svg", width=8, height=3.5)
# #windows(8,3.5)
# grid.arrange(p1, p2, p3, ncol=3)
# dev.off()


# inactive

p4<-rich%>%  filter(Fraction=="Inactive") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_fill_manual(values = IBM)+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "D",
       x="",
       y= "Inactive Shannon Diversity")
p4

p5<-rich%>%  filter(Fraction=="Inactive") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_fill_manual(values=IBM) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "E",
       x="",
       y= "Inactive ASV richness")
p5

p6<-rich%>%  filter(Fraction=="Inactive") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=Treatment, y=evenness,  fill=Treatment))+
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  scale_color_manual(values=IBM) +
  scale_fill_manual(values = IBM)+
  geom_jitter(size=1.5)+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0),legend.position="none")+
  labs(title = "F",
       x="",
       y= "Inactive Pilou's Evenness")

p3

require(gridExtra)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement/diversity_active_inact")
svg("inactive_ac_diversity.svg", width=8, height=7)
#windows(8,3.5)
grid.arrange(p1, p2, p3,p4, p5, p6, ncol=3)
dev.off()


##DIVERSITY STATS######

# active
# remove soil
df <- rich %>% filter(Fraction=="Active") %>% filter(Treatment!="Soil")
df$Treatment   <- factor(df$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# Shannon ANOVA 
anova1<- aov(Shannon~ Treatment+block, data = df)
summary(anova1)

#observed ANOVA
m1<-lm(Observed ~ Treatment+block,  data = df)
summary(m1)

#chao1
m1<-aov(Chao1 ~ Treatment+block,  data = df)
summary(m1)

# Evenness
anova1<- aov(evenness ~ Treatment+block, data = df)
summary(anova1)

# inactive
# remove soil
df <- rich %>% filter(Fraction=="Inactive") %>% filter(Treatment!="Soil")
df$Treatment   <- factor(df$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# Shannon ANOVA 
anova1<- aov(Shannon~ Treatment+block, data = df)
summary(anova1)

#observed ANOVA
m1<-lm(Observed ~ Treatment+block,  data = df)
summary(m1)

#chao1
m1<-aov(Chao1 ~ Treatment+block,  data = df)
summary(m1)

# Evenness
anova1<- aov(evenness ~ Treatment+block, data = df)
summary(anova1)





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
  
    

##### CAP ############  
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
  cap_result <- capscale(dist_matrix ~ Treatment*Fraction+block,
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
  
  setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")    
 # svg("cap.svg", width = 8 , height = 4.5)

  par(mfrow = c(1, 2)) # 1 row, 2 column
 # windows(4,4)
  par(cex.lab = 1.1) # make all fonts in graphs little bigger
  ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
           main="", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "C")
  par(adj=.5)
  points(sc_si, 
         col= "black",
         pch= c(21,22)[as.factor(metadat2$Fraction)],
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
         pch=c(1,15 ),
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
  
  commdata = as.data.frame(otu_table(ps1))
  # 2. Run the CAP (db-RDA) analysis
  # Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
  cap_result <- capscale(dist_matrix ~ Treatment+ Condition(block),
                         data = metadat2,
                         comm = commdata,
                         add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA
  
  anova.cca(cap_result, by="terms")
  
# . Perform all Pairwise Comparisons
 #  # The function will iterate through all pairs of the 'Habitat' factor
 library(BiodiversityR)
   pairwise_results <- multiconstrained(
     formula = dist_matrix ~ Treatment+ Condition(block),  # Same formula as the main CAP model
     data = metadat2,
     constrained = capscale,       # Specify the constrained ordination method
     permutations = 999            # Number of permutations for the test
   )

# view the raw pairwise results
   print(pairwise_results)
# Extract the raw p-values from the results
 raw_pvalues <- pairwise_results[, "Pr(>F)"]
 adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")

 # make table 
    table <- data.frame(
     Group1 = str_split_i(rownames(pairwise_results), "vs. ", 1),
     Group2 = str_split_i(rownames(pairwise_results), "vs. ", 2),
     Pseudo_F = pairwise_results[, "F"],
     Raw_P = raw_pvalues,
     p_value = round(adjusted_pvalues,4)
   )
   print(table)
   
  # order groups
   table$Group1<-factor(table$Group1,levels= c("L ", "G ",  "B ", "GB " , "LB ", "LG " ) )
   table$Group1
   factor(table$Group2)
   table$Group2[order(table$Group2)]
   table$Group2
   table$Group2<-factor(table$Group2,levels= c( "G", "B", "GB", "LB", "LG", "LGB") )
  
   #plot heatmap
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
   
   
   setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_CAPactive")
   pdf("Pval_active.pdf",  width=3.5, height=2.5)
   p_recol_covercrop
   dev.off() 
   
   
   
   
  # 1. Create an empty plot frame (type = "n")
  plot(cap_result, type = "n", scaling = 2)
  points(cap_result, display = "sites", pch = 16)
  
  ### 3. grab info for the plot
  smry <- summary(cap_result)
  smry
  sc_si <- scores(cap_result, display="sites", choices=c(1,2), scaling=2)
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
  ordiplot(cap_result, choices=c(1,2), scaling =2, type="none",
           main="", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "D")
  par(adj=.5)
  points(sc_si, 
         col= IBM[metadat2$Treatment],
         pch= 21,
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
  
  
 
  # dev.off()
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
  
  # plot abundance of these taxa across treatments

  # get data
ps1
  ##add add abundance and taxon infoget
    asvkp<-unique(top_spp$asv)
    
    # make df relative abundance
    df<-as.data.frame((otu_table(ps1)))
    df<-df/rowSums(df)
    df<- as.data.frame(t(df))
    df$asv<-row.names(df)
    df <- df[df$asv %in% asvkp, ]
    
    
    # add treatment info
    df$asv = NULL
    df<-as.data.frame(t(df))
    metadat2<-metadat2 %>% select(Treatment)
    df<-cbind(df, metadat2)
     
      # summarize abundance info by each treatment 
     df<-df %>% group_by(Treatment) %>%
       summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
     df<-as.data.frame(t(df))
     df
     colnames(df) <- c( "L",  "G"  ,"B" ,"GB" ,"LB" ,"LG", "LGB")
     df<-df[-1,] # remove first row which is row names
     df$asv<-row.names(df)
    
       # add taxa info
    tax<-as.data.frame((tax_table(ps1)))
    df<-left_join(df, tax)
    df<-left_join(df, species_scores)
    df

    df=df %>%
    pivot_longer(cols = L:LGB, 
                 names_to = "Treatment",
                 values_to = "percent_abundance")
    df<-as.data.frame(df)
    df$percent_abundance <- as.numeric(df$percent_abundance)
    df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))
    # label
    df$label <- df$Species
    df$label[which(df$label==""  )] <- df$Genus[which(df$label=="" )]
    
    # plot
    mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A",  "#FF7F00")
   

      ggplot(df)+
      geom_bar(aes(x=asv, y=percent_abundance, fill= Phyla), 
               stat="identity", position="dodge")+
      #geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
      #                 ymax=LFC_FractionViable_Cell+se_viable))+ 
      scale_fill_manual(values= mycols)+
      theme_bw(base_size = 12) +
      facet_grid(~Treatment, scales="free", space="free")+
      theme(axis.text.x = element_text(angle=60, hjust=1),
            plot.title = element_text(hjust = 0.5))+
      xlab("ASVS with with highest effect in CAP")
  
      
      # summary 
      df1<-df %>% group_by(Phyla, label, Treatment, ) %>%
        summarise(percent_abundance = sum(percent_abundance))
        
      
      # plot
      mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A",  "#FF7F00")
      
      
      ggplot(df1)+
        geom_bar(aes(x=label, y=percent_abundance, fill= Phyla), 
                 stat="identity", position="dodge")+
        #geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
        #                 ymax=LFC_FractionViable_Cell+se_viable))+ 
        scale_fill_manual(values= mycols)+
        theme_bw(base_size = 12) +
        facet_grid(~Treatment, scales="free", space="free")+
        theme(axis.text.x = element_text(angle=60, hjust=1.1),
              plot.title = element_text(hjust = 0.5))+
        xlab("ASVS with with highest effect in CAP")
##### CAP inactive ##################
  #  Constrained ordination
  ps1 <-subset_samples(ps, Fraction=="Inactive" & Treatment!="Soil" & Treatment!="CTL")
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  ps1
  # 1845 taxa
  # subset metadata
  metadat2<-filter(metadat, Fraction=="Inactive" & Treatment!="Soil" & Treatment!="CTL")
  #factor
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
  metadat2$Fraction   <- factor(metadat2$Fraction)
  
  # 1. Calculate the distance matrix (e.g., Bray-Curtis)
  dist_matrix<-vegdist(otu_table(ps1), method = "bray")
  
  # 2. Run the CAP (db-RDA) analysis
  # Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
  cap_result <- capscale(dist_matrix ~ Treatment + Condition(block),
                         data = metadat2,
                         add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA
  
  anova.cca(cap_result, by="terms")
  
  # 6. Perform all Pairwise Comparisons
  # The function will iterate through all pairs of the 'Habitat' factor
library(BiodiversityR)
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
  #adjusted_pvalues <- p.adjust(raw_pvalues, method = "bonferroni")
  adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")

  #?p.adjust
  # 10. Combine the results for final interpretation
  final_table <- data.frame(
    Pair = rownames(pairwise_results),
    Pseudo_F = pairwise_results[, "F"],
    Raw_P = raw_pvalues,
    fdr_adj_p = adjusted_pvalues
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
  
  setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/supplement") 
  svg("cap.inactive.svg", width = 5 , height = 5)
  par(cex.lab = 1.2) # make all fonts in graphs little bigger
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
           main="", cex = 1.2,
           xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
           ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
  par(adj = 0)
  title(main= "A")
  par(adj=.5)
  points(sc_si, 
         col= IBM[metadat2$Treatment],
         pch= 21,
         lwd=1,cex=1,
         bg=IBM[metadat2$Treatment])
  ordiellipse(sc_si, metadat2$Treatment,  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              col= IBM,
              alpha = 40,
              cex=1.2)
  dev.off()
##### PCOA   ########
  
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

  
 
#######PCOA plot active##
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



########## PERMANOVA
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





####PERCENT abundance figure active: ####
# grab data
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
# remove true singletons 
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
active<-subset_samples(ps, Fraction=="Active" & Treatment!="Soil")
active<-prune_taxa(taxa_sums(active) > 50, active)
active

taxon<- as.data.frame(tax_table(active))
df<-as.data.frame(otu_table(active))

# add metdata
metadat<-as.data.frame(as.matrix(sample_data(active)))
metadat
metadat<- metadat %>% select(SampleID, Trt_ID)
head(metadat)
n<-metadat$Trt_ID
n
# rename columns if you want, not required. 
colnames(df)
length(n)
length(colnames(df))
df

# make it percent
df<-(df/rowSums(df))*100
df<-as.data.frame(t(df))
colnames(df)<-n # edit col names
df.taxa<-cbind(df, taxon)



#rownames(df)<-NULL
# make rownames null
# summarize by phyla
df1<-aggregate(cbind( LG_N1,  LB_N1,  GB_N1,  LGB_N1, L_N2,   G_N2, 
                      B_N2,   LG_N2,  GB_N2,  L_N3  , G_N3,  
                      B_N3,   LG_N3,  GB_N3 , LGB_N3, L_N4 ,
                      G_N4 ,  B_N4,   LG_N4,  LB_N4,  GB_N4, 
                      L_N5,  G_N5,   B_N5,   LG_N5,  LB_N5, 
                      GB_N5,  L_N1,   G_N6,   B_N6,   LB_N6, 
                      GB_N6,  LGB_N6, G_N1  
  
                       
)~ Phyla, data = df.taxa, FUN = sum, na.rm = TRUE)

head(df1)

# gather by sample
df1<-  gather(df1, "sampleID", value, !starts_with("P") )
head(df1)
#remove zeros
df1<-df1[df1$value!=0,]
head(df1)


# make really low abundance taxa other
df1$Phyla[df1$value<1] <- "other"
df1


# add metadat
df1


# sum repeats
df0<-aggregate(cbind(value) ~ sampleID+Phyla, data = df1, FUN = sum, na.rm =TRUE)
head(df0)
df0<-df0[order(df0$sample),]
head(df0)

hist(df0$value)
mycols18<- c( "#1F78B4","#A6CEE3","#E31A1C",  "#FB9A99", "#33A02C","#B2DF8A",  "#FF7F00",  "#FDBF6F", "#6A3D9A" , "#CAB2D6",
               "#B15928", "#FFFF99",  "#eb05db","#edceeb","#1a635a","#9ad6ce" , "#969696", "#232423")

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/")
#svg(file="percent_barplot.svg",width = 12, height=10)
#windows(12,12)
df1%>% 
  ggplot(aes(fill=Phyla, y=value, x=sampleID)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values=mycols18) +
  #scale_fill_viridis(discrete = TRUE) +
  #ggtitle("Top phyla") +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#dev.off()

###### most abundant asvs by treatment



##### ANCOM###############


# remove rare taxa #
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
# remove true singletons 
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 50, ])
ps<-prune_taxa(keep, ps)
# at least 50 reads in active 
active<-subset_samples(ps, Fraction=="Active")
active<-prune_taxa(taxa_sums(active) > 50, active)
keep<-row.names(t(otu_table(active)))
ps<-prune_taxa(keep, ps)
ps 




#Legume
ps1<-subset_samples(ps , Treatment=="L")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1 # 1349
out1 = ancombc(data = ps1, assay_name = "counts", 
               tax_level = "asv",  
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = .10, lib_cut = 0, 
               group = "Fraction", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.001, global = FALSE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
get.table <- function(res, res_global) {
#log fold change
tab_lfc = res$lfc
head(tab_lfc)
dim(tab_lfc)
col_name = c("asv", "LFC_Intercept", "LFC_FractionViable_Cell")
colnames(tab_lfc) = col_name
# standard error
tab_se = res$se
col_name = c("asv", "se_Intercept", "se_viable")
colnames(tab_se) = col_name
head(tab_se)
tab_se<-as.data.frame(tab_se)
#test statistcs W maybe it's willcoxin?
tab_w = res$W
col_name = c("asv", "W_Intercept", "W_viable")
colnames(tab_w) = col_name
# P-values from the Primary Result
tab_p = res$p_val
col_name = c("asv", "p_Intercept", "p_viable")
colnames(tab_p) = col_name
head(tab_p)
#Adjusted p-values from the Primary Result"
tab_q = res$q
head(tab_q)
col_name = c("asv", "adj_p_Intercept", "adj_p_viable")
colnames(tab_q) = col_name
head(tab_q)

# yes or no is a taxa differentially abundant
tab_diff = res$diff_abn
col_name = c("asv", "DA_Intercept", "DA_Fraction_inactive_cells_Active")
colnames(tab_diff) = col_name

# through all togetha nd remove anything with DNA
tab <-tab_lfc %>%
  left_join(., tab_se) %>%
  left_join(., tab_w ) %>%
  left_join(., tab_p) %>%
  left_join(., tab_q) %>%
  left_join(., tab_diff)
return(tab)
}
tab<-get.table(res, res_global)
L<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("L", length(DA_Fraction_inactive_cells_Active)))
head(L)

#brass
ps1<-subset_samples(ps , Treatment=="B")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
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
tab<-get.table(res, res_global)
B<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("B", length(DA_Fraction_inactive_cells_Active)))
head(B)

#grass
ps1<-subset_samples(ps , Treatment=="G")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
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
tab<-get.table(res, res_global)
G<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("G", length(DA_Fraction_inactive_cells_Active)))


# grass brass
ps1<-subset_samples(ps , Treatment=="GB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
out1 = ancombc(data = ps1, assay_name = "counts", 
               tax_level = "asv",  
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = .10, lib_cut = 0, 
               group = "Fraction", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.001, global = FALSE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
tab<-get.table(res, res_global)
GB<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("GB", length(DA_Fraction_inactive_cells_Active)))


# legume brass
ps1<-subset_samples(ps , Treatment=="LB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
out1 = ancombc(data = ps1, assay_name = "counts", 
               tax_level = "asv",  
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = .10, lib_cut = 0, 
               group = "Fraction", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.001, global = FALSE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
tab<-get.table(res, res_global)
LB<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("LB", length(DA_Fraction_inactive_cells_Active)))


# legume grass
ps1<-subset_samples(ps , Treatment=="LG")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
out1 = ancombc(data = ps1, assay_name = "counts", 
               tax_level = "asv",  
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = .10, lib_cut = 0, 
               group = "Fraction", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.001, global = FALSE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
tab<-get.table(res, res_global)
LG<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("LG", length(DA_Fraction_inactive_cells_Active)))


# legume grass brass
ps1<-subset_samples(ps , Treatment=="LGB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
out1 = ancombc(data = ps1, assay_name = "counts", 
               tax_level = "asv",  
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = .10, lib_cut = 0, 
               group = "Fraction", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.001, global = FALSE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
tab<-get.table(res, res_global)
LGB<-tab %>% filter(DA_Fraction_inactive_cells_Active=="TRUE") %>%
  mutate(Treatment=rep("LGB", length(DA_Fraction_inactive_cells_Active)))

#rbind all and save table


tab<-rbind(L, G, B, GB, LB, LG, LGB)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
write.table(tab, "ancom_table.txt")


######ANCOM add abundance info ####
setwd("C:/Users/jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
tab<-read.table("ancom_table.txt", header = TRUE)
head(tab)
library(tidyverse)

# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="G")
ps1<-subset_samples(ps , Treatment=="G")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon infoget
get.taxa<- function(tab1, ps1){
asvkp<-unique(tab1$asv)

# make df relative abundance
df<-as.data.frame((otu_table(ps1)))
df<-df/rowSums(df)
df<- as.data.frame(t(df))

df$asv<-row.names(df)
df <- df[df$asv %in% asvkp, ]
tax<-as.data.frame((tax_table(ps1)))

df<-left_join(df, tax)
df<-left_join(df, tab1)
return(df)
}
df<-get.taxa(tab1,ps1)
df
### summarise the abundance in active and inactive
get.abund<-function(df, trt){
  
  
df$rhizo.inactive.mean <-   rowMeans(df %>% dplyr::select(starts_with("i_"))) %>% glimpse()
t<-df %>% select(.,starts_with("i_"))
sd_inactive<- apply(t, 1, sd, na.rm=TRUE)
sd_inactive
df$sd_inactive <- sd_inactive
# active  
df$rhizo.inactive.mean <-   rowMeans(df %>% dplyr::select(.,starts_with("i_"))) %>% glimpse()
t<-df %>% select(.,starts_with("i_"))
sd_inactive<- apply(t, 1, sd, na.rm=TRUE)
sd_inactive
df$sd_inactive <- sd_inactive
# active
df$active.mean <-   rowMeans(df %>% dplyr::select(contains("BCAT"))) %>%   glimpse()
t<-df %>% select(contains("BCAT"))
sd_active<- apply(t, 1, sd, na.rm=TRUE)
sd_active
df$sd_active <- sd_active
df<-df %>% select(-contains("i_"))%>% select(-contains("BCAT")) %>% select(-contains("ctl")) %>% select(-contains("fc_"))
df<-df%>%mutate(Treatment=rep(trt, length(asv)))
return(df)
}
G<-get.abund(df, "G")
G



# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="B")
ps1<-subset_samples(ps , Treatment=="B")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon
df<-get.taxa(tab1,ps1)
### summarise the abundance 
B<-get.abund(df, "B")
B

# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="L")
ps1<-subset_samples(ps , Treatment=="L")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon
df<-get.taxa(tab1,ps1)
### summarise the abundance 
L<-get.abund(df, "L")
L

# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="LB")
ps1<-subset_samples(ps , Treatment=="LB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon
df<-get.taxa(tab1,ps1)
### summarise the abundance 
LB<-get.abund(df, "LB")
LB


# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="LG")
ps1<-subset_samples(ps , Treatment=="LG")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon
df<-get.taxa(tab1,ps1)
### summarise the abundance 
LG<-get.abund(df, "LG")
LG

# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="GB")
ps1<-subset_samples(ps , Treatment=="GB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon
df<-get.taxa(tab1,ps1)
### summarise the abundance 
GB<-get.abund(df, "GB")
GB

# filter+  get abudance info
tab1<- tab%>% filter(Treatment=="LGB")
ps1<-subset_samples(ps , Treatment=="LGB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
##add add abundance and taxon
df<-get.taxa(tab1,ps1)
### summarise the abundance 
LGB<-get.abund(df, "LGB")
LGB


######ANCOM load data ####
tab<-rbind(L, G, B, GB, LB, LG, LGB)
setwd("C:/Users/jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data")
tab <-  read.csv("ancom_table_abund.csv")
tab$Treatment<-factor(tab$Treatment, levels=c("L", "G", "B", "GB", "LB", "LG", "LGB"))




# add neg values
head(tab)
tablong<-tab %>% mutate(
  "inactive.mean" = -1*rhizo.inactive.mean
) %>% pivot_longer(cols=c(inactive.mean, active.mean), names_to = "fraction", values_to = "mean")

# summarise
tablong<-tablong %>% group_by(Phyla,Treatment, Species, fraction)%>%
  summarise(
    mean= sum(mean)
  )

tablong$Phyla<-gsub("p__", "", tablong$Phyla)


mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A",  "#FF7F00")
mycols18<- c( "#1F78B4","#A6CEE3","#E31A1C",  "#FB9A99", "#33A02C","#B2DF8A",  "#FF7F00",  "#FDBF6F", "#6A3D9A" , "#CAB2D6",
              "grey", "#FFFF99",  "#eb05db","#edceeb","#1a635a","#9ad6ce" , "#969696", "#232423")
#plot


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_predict_mbiome")
svg("DA.svg", width = 9, height = 6)
ggplot(tablong)+
  geom_bar(aes(x=Species, y=mean, fill= Phyla), 
           stat="identity", position="dodge")+
  #geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
   #                 ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values= mycols18)+
  theme_bw(base_size = 12) +
  facet_grid(~Treatment, scales="free", space="free")+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  #coord_flip()+
   ylim(c(-500,900)) +
  xlab("Differentially Abundant Asvs")+
  ylab("number of reads")
dev.off()




##############predicting microbial community #############

# filter ps object only abundant taxa and active 
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
# remove true singletons 
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps
#only active 
ps<-subset_samples(ps , Fraction=="Active")
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps # 1804 taxa
tax_table<-tax_table(ps)

# import biomass data 

#load libraries
library(readxl)
library(tidyverse)
setwd("C:/Users/jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("percent.biomass.csv")
# only n+ 
biomass<-biomass %>% filter(N==1)


# LG
# filter for L community
ps1<-subset_samples(ps , Treatment=="L")
#ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # no rep 6
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LG") %>% filter(Species=="legume") %>% arrange(Pot_ID) %>% 
  filter(Rep!=6)# remove rep 6 
vector<-sp1$percent/100
data_frame_multiplied1 <- df %>%
  mutate(
    across(
     everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
# filter for G community
ps1<-subset_samples(ps , Treatment=="G" & Rep!="6")
#ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # no rep 6 for L so we will skip for G
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by G biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LG") %>% filter(Species=="grass") %>% arrange(Pot_ID) %>% 
  filter(Rep!=6)# remove rep 6 
vector<-sp1$percent/100
sp1
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
row.names(LG)<-c("predict_LG_N1", "predict_LG_N2", "predict_LG_N3", "predict_LG_N4", "predict_LG_N5")








# GB
# filter for G community
ps1<-subset_samples(ps , Treatment=="G" & Rep!="1") # no rep 1 for brass so skip that one 
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by G biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="GB") %>% filter(Species=="grass") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1)# remove rep 1
vector<-sp1$percent/100
vector
data_frame_multiplied1 <- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )

print(data_frame_multiplied1)

# filter for B community
ps1<-subset_samples(ps , Treatment=="B" )
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # no rep 1 for L so we will skip for 1
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by G biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="GB") %>% filter(Species=="brassica") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1)# remove rep 1
vector<-sp1$percent/100
sp1
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
row.names(GB)<-c("predict_GB_N2", "predict_GB_N3", "predict_GB_N4", "predict_GB_N5", "predict_GB_N6")




# LB
# filter for L community
ps1<-subset_samples(ps , Treatment=="L" & Rep!="1" & Rep!="6") # no rep 1 for brass so skip that one 
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LB") %>% filter(Species=="legume") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1 )# remove rep 1 a LB doesn't have a rep 5, so use 6 insteasd
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
# filter for B community
ps1<-subset_samples(ps , Treatment=="B" & Rep!="6" )
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # no rep 1 for B so we will skip for 1 and also skip rep 6
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by B biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LB") %>% filter(Species=="brassica") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1)# remove rep 1
vector<-sp1$percent/100
sp1
vector
data_frame_multiplied2<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
#print(data_frame_multiplied2)
# add together
LB<-data_frame_multiplied1 + data_frame_multiplied2
row.names(LB)<-c("predict_LB_N2", "predict_LB_N3", "predict_LB_N4", "predict_LB_N5")




# LBG 
# filter for L community
ps1<-subset_samples(ps , Treatment=="L" & Rep!="1" & Rep!="6") # no rep 1 for brass so skip that one 
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by L biomass
head(biomass)
sp1<-biomass %>% filter(Treatment=="LGB") %>% filter(Species=="legume") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1 & Rep!=6)# remove rep 1  and 6
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

# filter for B community
ps1<-subset_samples(ps , Treatment=="B" & Rep!="6" )
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # no rep 1 for B so we will skip for 1 and also skip rep 6
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
#names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by B biomass
#head(biomass)
sp1<-biomass %>% filter(Treatment=="LGB") %>% filter(Species=="brassica") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1 & Rep!=6)# remove rep 1 and 6
vector<-sp1$percent/100
#sp1
#vector
data_frame_multiplied2<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )
#print(data_frame_multiplied2)

# filter for G community
ps1<-subset_samples(ps , Treatment=="G" & Rep!="6" & Rep!="1")
df<-as.data.frame(otu_table(ps1))
#as.data.frame(sample_data(ps1))   # no rep 1 and no rep 6
# arrange by rep
names<-row.names(df)
names<-str_split(names, "_", simplify = TRUE)
#names
df$Pot_ID<-as.numeric(names[,2])
df<-df %>% arrange(Pot_ID)
#df$Pot_ID
df$Pot_ID=NULL
# get biomass info and multiply by B biomass
#head(biomass)
sp1<-biomass %>% filter(Treatment=="LGB") %>% filter(Species=="brassica") %>% arrange(Pot_ID) %>% 
  filter(Rep!=1 & Rep!=6)# remove rep 1
vector<-sp1$percent/100
#sp1
#vector
data_frame_multiplied3<- df %>%
  mutate(
    across(
      everything(),  # Selects columns starting with "Value"
      .fns = ~ .x * vector # The function to apply: multiply the current column (.x) by the vector
    )
  )

# add together
LGB<-data_frame_multiplied1 + data_frame_multiplied2 + data_frame_multiplied3
row.names(LGB)<-c("predict_LGB_N2", "predict_LGB_N3", "predict_LGB_N4", "predict_LGB_N5")

df1<-rbind(LG, GB, LB, LGB)
row.names(df1)



#### write files
df<-as.data.frame(otu_table(ps))
otus<-rbind(df1, df)
otus<-t(otus)
setwd("C:/Users/jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/predicted")
write.csv(otus, "feature.table.csv")


#taxon
tax_table<-as.data.frame(tax_table(ps))
write.csv(tax_table, "taxonomy.csv")


##### import predicted #####
IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)
## Set the working directory ###
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/predicted")
taxon <- read.csv("taxonomy.csv", row.names = 1)
asvs <- read.csv("feature.table.csv", row.names = 1)
metadat<-read.csv("metadata_predicted16S.csv", header = T)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-as.data.frame(t(asvs))
asvs[1:5,1:5]#taxa are columns

## order asvs
n<-rownames(asvs)
asvs$SampleID <- n
asvs<-as.data.frame(asvs[order(asvs$SampleID),])
asvs[1:5,1:5]
asvs$SampleID <- NULL


## order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
row.names(metadat) <- metadat$SampleID
metadat
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
#factor(metadat$Treatment)
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes present", "Legumes absent"))

# import it phyloseq

Workshop_OTU <- otu_table(as.matrix(asvs), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 1804 taxa when rarefied 


##### predicted plot #####
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1745 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL" )
#factor
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
#metadat2$Fraction   <- factor(metadat2$Fraction)

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
colnames(otus.p) <- c("PC1", "PC2")

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-round(perc.exp[1],2)
pe2<-round(perc.exp[2],2)
pe2

#calculate total variance explained by each principal component
perc.exp<-otus.eig/(sum(otus.eig))*100
#scree plot 
plot(otus.pcoa$eig)


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
svg("pcoa.predicted1.svg", height = 4, width = 4)
par(adj=.5)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="",
         xlab=paste("PCoA1 (",pe1,"% var. explained)"),
         ylab=paste("PCoA2 (",pe2,"% var. explained)"))
  par(adj = 0)
title(main= "E")
par(adj=.5)
points(otus.p, 
       col= IBM[as.factor(metadat2$Treatment)],
       pch= c(16,8)[as.factor(metadat2$Measurement)],
       
       lwd=2,cex=1.2,
       bg=IBM[as.factor(metadat2$Treatment)],)
# ordiellipse(otus.pcoa, as.factor(metadat2$Measurement),  
#             kind = "ehull", conf=0.95, label=T, 
#             draw = "polygon",
#             border = 0,
#             col= IBM,
#             alpha = 50,
#             cex=1.5)
# legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
#        fill= IBM,
#        cex=1,
#        bty = "n")
legend("bottomright", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=1,
       title = "",     bty = "o")

dev.off()
# permanova
adonis2(otus.bray ~ Treatment, data = metadat2)
adonis2(otus.bray ~ Treatment+Measurement+Treatment*Measurement, data = metadat2, by="terms")

# seperatated out by treatment
#GB
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment=="GB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1745 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment=="GB")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)

# permanova
adonis2(otus.bray ~ Measurement, data = metadat2)


#LB
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment=="LB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1745 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment=="LB")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2)



#LG
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment=="LG")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1745 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment=="LG")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(5-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2) # marginal



#LGB
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment=="LGB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1745 taxa
# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment=="LGB")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(5-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2) # significantly different

#### CAP ##### 

# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(otus.bray ~ (Treatment+Measurement)^2,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

anova.cca(cap_result, by="terms")

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
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
svg("cap.active.svg", width = 5 , height = 5)
#windows(6,6)
par(cex.lab = 1.1) # make all fonts in graphs little bigger
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         main="", cex = 1.2,
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj = 0)
title(main= "E")
par(adj=.5)
points(sc_si, 
       col= IBM[metadat2$Treatment],
       pch= c(16,8)[as.factor(metadat2$Measurement)],
       lwd=1,cex=2,
       bg=IBM[metadat2$Treatment])
ordiellipse(sc_si, metadat2$Treatment,
            kind = "ehull", conf=0.95, label=F,
            draw = "polygon",
            border = 0,
            col= IBM,
            alpha = 40,
            cex=1)
# legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
#        fill= IBM,
#        cex=1,
#        bty = "n")
legend("bottomleft", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=1,
       title = "",     bty = "o")

dev.off()







#### extract PC1 ##########


# filter data 
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps
# filter for active fraction 
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1 #1785 taxa
# filter metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil" & Treatment!="CTL")
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction   <- factor(metadat2$Fraction)

# 1. Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# 1. Perform PCA
# scale. = TRUE is critical to ensure all variables are treated equally
pca_result <- prcomp(dist_matrix, center = TRUE, scale. = TRUE)

# 2. Extract the first principal component (PC1)
# This is your new one-dimensional variable
pc1_variable <- pca_result$x[, 1]
pc2_variable <- pca_result$x[,2]

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

# # n fix 
# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
# nfix<-read.csv("Nfix.csv")
# head(nfix)
# nfix<-nfix %>% 
#   filter(duplicate=="N") %>%select(Trt_ID, perc.Ndfa, n_fix_per_legume )
# df<-left_join(df, nfix)

# weeds
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed<-read.csv("weed_seed_decay.csv")
head(weed)
weed<-weed %>% select(Trt_ID, foxtail_prop_nongerm, pigweed_prop_nongerm)
df<-left_join(df, weed)

# z transform
df$z_Root.Biomass<-as.vector(scale(df$Root.Biomass))
df$z_Shoot.Biomass<-as.vector(scale(df$Shoot.Biomass))
df$z_foxtail_prop_nongerm<-as.vector(scale(log(df$foxtail_prop_nongerm+1)))
df$z_pigweed_prop_nongerm<-as.vector(scale(log(df$pigweed_prop_nongerm+1)))
df$z_boncat_freq<-as.vector(scale(df$boncat_freq))
df$z_active_cel_per_g<-as.vector(scale(log(df$active_cel_per_g+1)))
#df$z_perc.Ndfa<-as.vector(scale(df$perc.Ndfa))
#df$z_n_fix_per_legume<-as.vector(scale(df$n_fix_per_legume))

# hist 
hist(df$z_Shoot.Biomass) # okay 
hist(df$z_pigweed_prop_nongerm) # okay
hist(df$z_boncat_freq) # okay
hist(df$z_active_cel_per_g) # okay after log transform
#hist(df$z_perc.Ndfa) # okay
#hist(df$z_n_fix_per_legume)


# 4. Add it to your data frame
df$PC1 <- pc1_variable
df$PC2 <- pc2_variable

# corrplots

par(mfrow=c(2,2))
plot(df$PC1, df$z_Shoot.Biomass, main = "active, not sig")
plot(df$PC1, df$z_Root.Biomass, main = "acitve, p=.07, rsq =.07")
plot(df$PC1, df$z_boncat_freq, main = "active, not sig")
plot(df$PC1, df$z_active_cel_per_g, main = "active, p=.01, rsq=.25")

par(mfrow=c(1,2))
plot(df$PC1, df$z_pigweed_prop_nongerm, main= "active p=.06, Rsq=.07")


###ggplot 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig6_mbiome_functional")
svg("shoot_pc1.svg", width=6, height=3)
ggplot(df, aes(y=z_Shoot.Biomass, x=PC1, col=Treatment))+
  geom_point()+
  theme_bw(base_size = 12)+
  scale_color_manual(values=IBM)+
  labs(y = "Z transformed shoot biomass",
       x = "PC1 active microbiome")
dev.off()

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig6_mbiome_functional")
svg("shoot_pc1.svg", width=6, height=3)
ggplot(df, aes(y=z_Shoot.Biomass, x=PC2, col=Treatment))+
  geom_point()+
  theme_bw(base_size = 12)+
  scale_color_manual(values=IBM)+
  labs(y = "Z transformed shoot biomass",
       x = "PC1 active microbiome")
dev.off()



svg("weed_pc1.svg", width=4.5, height=3.5)
ggplot(df, aes(y=z_pigweed_prop_nongerm, x=PC2))+
  geom_point()+
  geom_smooth(method=lm)+
  theme_bw(base_size = 12)+
  labs(y = "non germinating pigweed",
       x = "PC1 active microbiome")
  #annotate("text", x = 1, y = .5, label = "p=0.05, Rsq=.11", 
  #         color = "black", size = 5, fontface = "bold")
  
dev.off()



svg("weed_pc1.svg", width=4.5, height=3.5)
ggplot(df, aes(y=z_pigweed_prop_nongerm, x=PC1))+
  geom_point()+
  geom_smooth(method=lm)+
  theme_bw(base_size = 12)+
  labs(y = "non germinating pigweed",
       x = "PC1 active microbiome")
#annotate("text", x = 1, y = .5, label = "p=0.05, Rsq=.11", 
#         color = "black", size = 5, fontface = "bold")

dev.off()


# lm 
m1<-lm(df$PC1~df$z_Root.Biomass)
m<-summary(m1) # trend
m
m<-m$coefficients
m[4]



m1<-lm(df$PC1~z_Shoot.Biomass*Treatment, data=df)
summary(m1)

m1<-lm(df$PC1~df$z_boncat_freq)
summary(m1)

m1<-lm(df$PC1~df$z_active_cel_per_g)
summary(m1) # sig

m1<-lm(df$PC1~df$z_foxtail_prop_nongerm)
summary(m1)

m1<-lm(df$PC1~df$z_pigweed_prop_nongerm)
summary(m1) # trend

m1<-lm(df$z_pigweed_prop_nongerm~df$PC1)
summary(m1) # trend

m1<-lm(df$PC1~df$z_n_fix_per_legume)
summary(m1)

m1<-lm(df$z_perc.Ndfa)
summary(m1) 


######extract key taxa ####

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

setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
svg(filename="active.vectors.svg", height = 6, width = 6)
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
  labs(title = "Top 15 ASV Contributing to Active CAP Variation",
       x = "CAP1", y = "CAP2")
dev.off()
# plot abundance of these taxa across treatments

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
# tax
#
# tax<-tax %>% group_by(Phyla, Genus, asv) %>% summarise()
# # tax
# # tax$Blast_ID <- c(
#   "Uncultured Actinomycetes",
#   "Uncultured Actinomycetes" ,
#   "Uncultured Actinomycetes" ,
#   "Uncultured Cyanobacterium",
#   "Uncultured Cyanobacterium",
#   "Uncultured Cyanobacterium",
#   "Uncultured Cyanobacterium",
#   "Uncultured Cyanobacterium",
#   "Deinococcus sp.",
#   "Deinococcus sp.",
#   "Deinococcus sp.",
#   "Deinococcus sp.",
#   "Escherichia sp.",
#   "Escherichia sp.",
# #   "Rhizobium Leguminosarum")
# tax
# df<-left_join(df, tax)
# df<-as.data.frame(df)
# df
#
# df$percent_abundance <- as.numeric(df$percent_abundance)
# df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))
#
# #
# # # plot
# # mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A",  "#FF7F00")
# #
# # df$Phyla
# # df %>%
#   ggplot(aes(x=Treatment, y=percent_abundance, fill = Phyla))+
#   geom_boxplot(outliers=FALSE)+
#   geom_jitter()+
#   scale_fill_manual(values= mycols)+
#   theme_bw(base_size = 12) +
#   facet_grid( ~asv, scales="free", space="free")+
#   theme(axis.text.x = element_text(angle=60, hjust=1),
#         plot.title = element_text(hjust = 0.5))
#
#
# # aggregate by taxa after blasting the sequences
# df
#
# df<-df %>% group_by(Phyla, Treatment, Trt_ID, Blast_ID) %>%
#   summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
# df
#
# df %>%
#   ggplot(aes(x=Treatment, y=percent_abundance, fill = Phyla))+
#   geom_boxplot(outliers=FALSE)+
#   geom_jitter()+
#   scale_fill_manual(values= mycols)+
#   theme_bw(base_size = 12) +
#   facet_grid(~Blast_ID, scales="free", space="free")+
#   theme(axis.text.x = element_text(angle=60, hjust=1),
#         plot.title = element_text(hjust = 0.5))
#
#
#
# # top 20 asvs #############
#
#
# # Select the top 10 species
# top_spp <- head(species_scores, 20)
# top_spp
#
# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
# svg(filename="active.vectors20.svg", height = 6, width = 6)
# ggplot() +
#   # Draw a circle/origin cross for reference
#   geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
#   geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
#   # Add the vectors (arrows)
#   geom_segment(data = top_spp,
#                aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
#                arrow = arrow(length = unit(0.2, "cm")), color = "darkred") +
#
#   # Add labels with some padding
#   geom_text(data = top_spp,
#             aes(x = CAP1, y = CAP2, label = asv),
#             color = "black", fontface = "italic", vjust = -0.5) +
#
#   theme_bw() +
#   labs(title = "Top 15 ASV Contributing to Active CAP Variation",
#        x = "CAP1", y = "CAP2")
# dev.off()
# # plot abundance of these taxa across treatments
#
# # get data
# ps1
# ##add add abundance and taxon infoget
# asvkp<-unique(top_spp$asv)
#
# # make df relative abundance
# df<-as.data.frame((otu_table(ps1)))
# df<-df/rowSums(df)*100
# df<- as.data.frame(t(df))
# df$asv<-row.names(df)
# df <- df[df$asv %in% asvkp, ]
# df$asv <-NULL
#
# # add treatment info
# df<-as.data.frame(t(df))
# metadat3<-metadat2 %>% select(Treatment, Trt_ID)
# df<-cbind(df, metadat3)
# df
#
# df<-df %>%
#   pivot_longer(cols = c(-Treatment, -Trt_ID) ,
#                names_to = "asv",
#                values_to = "percent_abundance")
# df
# # add phyla info
# tax<-as.data.frame((tax_table(ps1)))
# tax <- tax[tax$asv %in% asvkp, ]
# tax<-tax %>% group_by(Phyla, Order, Family, Genus, Species, asv) %>% summarise()
#
#
# tax<-tax %>% group_by(Phyla, Genus, asv) %>% summarise()
# print(tax)
# tax$Blast_ID <- c("Actinomycetes",
#                   "Actinomycetes",
#                   "Actinomycetes",
#                   "Actinomycetes",
#                   "Terrabacter sp.",
#                   "Cyanobacterium",
#                   "Cyanobacterium",
#                   "Cyanobacterium",
#                   "Cyanobacterium",
#                   "Cyanobacterium",
#                   "Cyanobacterium",
#                   "Deinococcus sp.",
#                   "Deinococcus sp.",
#                   "Deinococcus sp.",
#                   "Deinococcus sp.",
#                   "Escherichia sp.",
#                   "Escherichia sp.",
#                   "Escherichia sp.",
#                   "Rhizobium sp.",
#                   "Rhizobium sp."
# )
# tax
# df<-left_join(df, tax)
# df<-as.data.frame(df)
# df
#
# df$percent_abundance <- as.numeric(df$percent_abundance)
# df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))
#
#
# # plot
# mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A",  "#FF7F00")
#
#
#
# # aggregate by taxa after blasting the sequences
# df
#
# df<-df %>% group_by(Phyla, Treatment, Trt_ID, Blast_ID) %>%
#   summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
# df
#
#
# setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_taxa")
# svg(filename="active.taxa20.svg", height = 4, width = 9)
# df %>%
#   ggplot(aes(x=Treatment, y=percent_abundance, fill = Phyla))+
#   geom_boxplot(outliers=FALSE)+
#   geom_jitter()+
#   scale_fill_manual(values= mycols)+
#   theme_bw(base_size = 12) +
#   facet_grid(~Blast_ID, scales="free", space="free")+
#   theme(axis.text.x = element_text(angle=60, hjust=1),
#         plot.title = element_text(hjust = 0.5))
# dev.off()
#
#

### genus overall ########
tax<-as.data.frame((tax_table(ps1)))

tax <- tax[tax$asv %in% asvkp, ]
tax<-tax %>% group_by(Phyla, Order, Family, Genus, Species, asv) %>% summarise()
print(tax)
target<-unique(tax$Genus)
  # Subset the phyloseq object
  ps1 <- subset_taxa(ps1, Genus %in% target)
  ps1
  
  
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
  df<-df %>% group_by(Genus) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
  df
  
  # save genus info
  Genus <- df$Genus
  df$Genus = NULL
  df
  # add treatment info
  df<-as.data.frame(t(df))
  colnames(df) <- Genus
  df
  
  metadat3<-metadat2 %>% select(Treatment, Trt_ID)
  df<-cbind(df, metadat3)
  df
  
  # pivot
  df<-df %>%
    pivot_longer(cols = starts_with(" g_"), 
                 names_to = "Genus",
                 values_to = "percent_abundance")
  df
  
  # add phyla info
  tax<-as.data.frame((tax_table(ps1)))
    tax <- tax[tax$Genus %in% target, ]
  tax
  tax<-tax %>% group_by(Phyla, Genus) %>% summarise()
  tax
  df<-left_join(df, tax)
  df<-as.data.frame(df)
  df
  
  df$Genus<-gsub(" g__", "", df$Genus)
  df$Genus<-gsub("_501058", "", df$Genus)
  df$Genus<-gsub("_C", "", df$Genus)
  
  
  # plot
  mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A",  "#FF7F00")
  df$percent_abundance <- as.numeric(df$percent_abundance)
  df$Treatment <- factor(df$Treatment, levels =c ("L", "G", "B", "GB", "LB", "LG", "LGB"))
  
  setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_predict_mbiome")
  #svg("active.genus.svg", height = 4, width=8)
  df %>%
    ggplot(aes(x=Treatment, y=percent_abundance, fill = Phyla))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter()+
    scale_fill_manual(values= mycols)+
    theme_classic(base_size = 12) +
    facet_grid(~Genus, scales="free", space="free")+
    theme(axis.text.x = element_text(angle=60, hjust=1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  #### global genus tests
  df1<-df %>% filter(Genus=="AC-14")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n1<-m1$p.value
  
  df1<-df %>% filter(Genus=="Deinococcus")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1)
  n2<-m1$p.value
    
  df1<-df %>% filter(Genus=="Escherichia")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) 
  n3<-m1$p.value
  
  df1<-df %>% filter(Genus=="Microcoleus")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n4<-m1$p.value

  df1<-df %>% filter(Genus=="Rhizobium") 
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # sig
  n5<-m1$p.value
  
  df1<-df %>% filter(Genus=="Vampirovibrio")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # sig
  n6<-m1$p.value
  
  
  # raw results
  raw_pvalues<-c( n1, n2, n3, n4,  n5, n6)
  taxa <-  unique(df$Genus) 
  test <- "Global"
  # pvalue Adjustment
  adjusted_pvalues <- p.adjust(raw_pvalues, method = "fdr")
  adjusted_pvalues
  final_table <- data.frame(
    taxa = taxa,
    test = test,
    Raw_P = raw_pvalues,
    fdr_adj_p = adjusted_pvalues
    #fdr_adj_p_nonsci = format(adjusted_pvalues, scientific = FALSE)
    
  )
  print(final_table)
  
    
## pairwise stats
  #   "Deinococcus_C"
  n<-unique(df$Genus)
  n
  df1<-df %>% filter(Genus=="Deinococcus")  %>% filter(Treatment=="L"| Treatment=="B")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n1<-m1$p.value
  
  n<-unique(df$Genus)
  df1<-df %>% filter(Genus=="Deinococcus")  %>% filter(Treatment=="L"| Treatment=="G")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n2<-m1$p.value 
  
  n<-unique(df$Genus)
  df1<-df %>% filter(Genus=="Deinococcus")  %>% filter(Treatment=="B"| Treatment=="G")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n3<-m1$p.value 
    
  df1<-df %>% filter(Genus=="Deinococcus")  %>% filter(Treatment=="B"| Treatment=="GB")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n4<-m1$p.value 
  
  
  df1<-df %>% filter(Genus=="Deinococcus")  %>% filter(Treatment=="B"| Treatment=="LB")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n5<-m1$p.value 
  
  # raw results
  raw_pvalues<-c( n1, n2, n3, n4,  n5)
  taxa <-  "Deinococcus_C"
  pair <- c("L-B", "L-G", "B-G", "B-GB", "B-LB")
  #  # 9. Apply the Holm (Holm-Bonferroni) Adjustment
  
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
  n
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
  
  n<-unique(df$Genus)
  df1<-df %>% filter(Genus=="Rhizobium")  %>% filter(Treatment=="L"| Treatment=="LG")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n4<-m1$p.value 
  
  n<-unique(df$Genus)
  df1<-df %>% filter(Genus=="Rhizobium")  %>% filter(Treatment=="L"| Treatment=="LB")
  m1<-kruskal.test(percent_abundance ~ Treatment, data = df1) # 
  n5<-m1$p.value 
  
  
  # raw results
  raw_pvalues<-c( n1, n2, n3, n4, n5)
  taxa <-  "Rhizobium_C"
  pair <- c("L-B", "L-G", "B-G", "L-LG", "L-LB")
  #  # 9. Apply the Holm (Holm-Bonferroni) Adjustment
  
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
  
  
### deseq ##########
  library(phyloseq)
  library(DESeq2)
  library(ggplot2)
  
  
  # 1. Convert phyloseq object to DESeq2 format
  # We add 1 to all counts to handle zeros (Laplace smoothing)
   # import it phyloseq
   Workshop_OTU <-otu_table(ps1)+1
   Workshop_metadat <- sample_data(ps1)
   Workshop_taxo <- tax_table(ps1) 
   ps.plusone <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
   ps.plusone
  
  asvkp
  
  # Focus only on specific genera
  ps.subset <- subset_taxa(ps.plusone, asv %in% asvkp)
 
    
  #  DESeq2 
  ds <- phyloseq_to_deseq2(ps.subset, ~ factor(Legume))
  # 'test="Wald"' is standard for two-group comparisons
  # 'fitType="local"' is often better for sparse microbial data 
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  res <- results(diagdds)
  # check results
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))  
  sigtab  
    
  
  # DEseq2
  ds <- phyloseq_to_deseq2(ps.subset, ~ Grass)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  res <- results(diagdds)
  # check
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  
  # DEseq
  ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  res <- results(diagdds)
  # check
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  
  # DEseq
  ds <- phyloseq_to_deseq2(ps.subset, ~ Brassicae*Grass)
  diagdds <- DESeq(ds, test="Wald", fitType="local")
  res <- results(diagdds)
  # check
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))  
  sigtab  
  
  
  
### aggregate 2 genus level
  # 1. Convert phyloseq object to DESeq2 format
  # We add 1 to all counts to handle zeros (Laplace smoothing)
  # import it phyloseq
  
  ps.subset <- subset_taxa(ps1, Genus %in% target)
  
  
  # aggregate by genus
  df<-as.data.frame(t(otu_table(ps.subset)))
  tax<- as.data.frame(tax_table(ps.subset))
  tax <- tax %>% select(-Confidence)
  df<-cbind(df,tax)
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
    group_by( Phyla,Class, Order, Family , Genus) %>%
    summarise()
  tax<-as.data.frame(tax)
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
  