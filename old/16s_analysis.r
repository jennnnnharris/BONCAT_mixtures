# 16S analysis
# total
# beta diversity
# BONCAT mixtures
# Jennifer Harris
# Jan 7 2025
# Last Updated: May 2025

#R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

### 1. Initial Setup ###

### Clear workspace ###

rm(list=ls())


################## Load required libraries ############

#basic
library(tidyverse)
library(vegan)
library(readxl)
library(lubridate)

library(phyloseq)
library(MicEco)
#library(fantaxtic)
#library(microbiome)
#library(DT)



# set colors
blues<-c( "#9CA9BAFF", "#5480B5FF", "#3D619DFF", "#405A95FF", "#345084FF")
mycols7<-c( "#4B2D4BFF", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
mycols4 <- c("#4B2D4BFF",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )
#df$Treatment   <- factor(df$Treatment, levels= c( "L", "LB", "LG", "LGB"))
mycols3<- c(  "#f4f1bb", "#ed6a5a","#9bc1bc")
mycols3<- c( "#006d77",  "#f4d35e", "#e94f37")



#####Import data#####
## Set the working directory ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")

taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat.csv", header = T)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns

## check the reads per sample ##
#rowSums(asvs)[order(rowSums(asvs))]

#T_DNA_23_S153 has really few reads so I am omitting it.
asvs<-asvs[which(row.names(asvs)!= "T_DNA_23_S153"),]

#get min number of reads in a sample
min.s<-min(rowSums(asvs))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs, min.s)

## order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
row.names(metadat) <- metadat$SampleID
metadat

# ck names in metadata file.
# length(intersect(colnames(asvs.raw) , metadat$SampleID)) # Apply setdiff function to see what's missing from the tree
# mynames<- setdiff(metadat$SampleID , colnames(asvs.raw) )
# length(mynames)
# mynames

#make taxon matrix row names OTUs
#taxon[1:5,1:5]
row.names(taxon) <- taxon$Feature.ID

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs.r), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps
# 329K taxa when rarefied 

####### rarefaction curve ######

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")

svg("observed.vs.rarefied.svg" )
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
dev.off()

svg("rarefaction.svg", height=8, width=8 )
rarecurve(i, step = 10000, sample = raremax, col = "blue", cex = 0.6)
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



# many very rare taxa


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
# 218 k asvs

# remove true singletons 
total<-prune_taxa(taxa_sums(total) > 1, total)
total
# 183095 

#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(total)))/nsamples(total)
keep<-row.names(t(otu_table(total))[ mean.reads > 5, ])
total<-prune_taxa(keep, total)
total
#1696 asvs

# remove asvs that are in less than 5 samples
total <-ps_prune(total, min.samples = 3)
total
#1696 asvs

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



###PCOA set shapes and cols######
  

  #mycols8<- c( "grey", "#1F78B4",  "#eb05db", "#33A02C", "#FF7F00","#1a635a", "#6A3D9A",  "gold")
  #mycols7<- c( "#1F78B4",  "#eb05db", "#33A02C", "#FF7F00","#1a635a", "#6A3D9A",  "gold")
  #antique<-c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", "#526A83FF", "#625377FF", "#68855CFF","#9C9C5EFF", "#A06177FF", "#8C785DFF", "#467378FF", "#7C7C7CFF")
  #mycols7 <- c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", "#526A83FF", "#625377FF", "#68855CFF")
  #mycols7 <-c("#4B2D4BFF", "#3C3C5AFF", "#4B6987FF", "#789696FF", "#968787FF", "#D2C3C3FF", "#875A2DFF", "#873C3CFF")
  #mycols7<-c("#78A5C3FF",  "#3C3C5AFF", "#5A784BFF", "#E1D2D2FF", "#694B5AFF", "#A55A2DFF", "#873C3CFF")
  #30025c
  #690061
  #970860
  #bd2d5b
  #db5255
  #f17951
  #ffa251
  #mycols7<-c("#E3C1CBFF", "#AD5A6BFF", "#C993A2FF", "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")
  
  
  myshapes <- c(1, 12,15 ,21, 22, 23 , 24)
  myshapes2 <- c(21 , 12, 24,1, 15 , 22, 23 )
  mycols3<- c(  "#f4f1bb", "#ed6a5a","#9bc1bc")
  mycols3.1<- c(  "#f4d35e", "#e94f37","#006d77")
  setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
  
  
#####PCOA overall singals  ########  
  # remove control
  ps2<-subset_samples(ps, Fraction !="CTL" )
  ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
  ps2
  
  # Calculate Bray-Curtis distance between samples
  otus.bray<-vegdist(otu_table(ps2), method = "bray")
  # Perform PCoA analysis of BC distances #
  otus.pcoa <- cmdscale(otus.bray, k=(44-1), eig=TRUE)
  
  # Store coordinates for first two axes in new variable #
  otus.p <- otus.pcoa$points[,1:2]
  
  # Calculate % variance explained by each axis #
  otus.eig<-otus.pcoa$eig
  perc.exp<-otus.eig/(sum(otus.eig))*100
  pe1<-perc.exp[1]
  pe2<-perc.exp[2]
  
    
  # subset metadata
  metadat2<-filter(metadat, Fraction!="CTL")
  
  #factor
  metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
  mycols3<- c(  "#f4f1bb", "#ed6a5a","#9bc1bc")
  mycols3.1<- c(  "#f4d35e", "#e94f37","#006d77")
  metadat2$Fraction   <- factor(metadat2$Fraction)
  
  #plot
  #windows(6,6)
  #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
  
  svg("pcoa.allfractions.svg",  width = 5, height = 5 )
  ordiplot(otus.pcoa,choices=c(1,2), type="none", main="all fractions ",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
           ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
  points(otus.p, 
         col= mycols3.1[metadat2$Fraction],
         pch= myshapes[metadat2$Rep],
         lwd=1,cex=1.5,
         bg=mycols3.1[metadat2$Fraction])
  
  ordiellipse(otus.pcoa, metadat2$Fraction,  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              #lwd=.1,
              col= mycols3,
              alpha = 30)
  dev.off()
  
  
  
   
  
#####PCOA  of total ######

# Calculate Bray-Curtis distance between samples
# Perform PCoA analysis of BC distances #
total <-subset_samples(total, Treatment !="Soil" )
otus.bray<-vegdist(otu_table(total), method = "bray")
otus.pcoa <- cmdscale(otus.bray, k=(44-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
colnames(otus.p) <- c("PC1", "PC2")
df.pcoa <- cbind(sample_data(total), otus.p)
# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")

#set factors 
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))


svg("pcoa.total.svg", height = 5, width =5)
windows(6,6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Total DNA",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, 
       col= mycols7[metadat2$Treatment],
       pch= myshapes[as.factor(metadat2$Rep)],
       lwd=2,cex=2,
       bg=mycols7[metadat2$Treatment])

legend("topleft", legend=c( "L", "G", "B", "GB", "LB", "LG", "LGB") ,
       fill= mycols7,
       cex=1,
       title = "Treatment",
       bty = "n")

ordiellipse(otus.pcoa, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            lwd=.1,
            col= mycols7,
            alpha = 40)
dev.off()


#####PLOTS PCOA nitrogen + no nitrogen  ##########

# Calculate Bray-Curtis distance between samples
# Perform PCoA analysis of BC distances #

total <-subset_samples(total, Treatment !="Soil" )
otus.bray<-vegdist(otu_table(total), method = "bray")
otus.pcoa <- cmdscale(otus.bray, k=(44-1), eig=TRUE)
otus.p <- otus.pcoa$points[,1:2]
colnames(otus.p) <- c("PC1", "PC2")

df.pcoa <- cbind(sample_data(total), otus.p)
# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")

#set factors 
#metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

df.pcoa$Treatment   <- factor(df.pcoa$Treatment, levels= c( "L", "G", "B", "GB", "LB", "LG", "LGB"))
factor(metadat2$Fraction)
df.pcoa$Rep <- factor(df.pcoa$Rep)
metadat2$N <- factor(metadat2$N)
myshapes2 <- c(21 , 12 )

# scree plot
otus.pcoa$eig
plot(perc.exp[1:8],
     ylab = "percent varience explained",
     xlab = "PC")

# Visualize ordination

ggplot(df.pcoa, aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  # Replace 'Group' with the metadata column you want to use for coloring
  geom_point(size = 3, alpha=.7) +
  theme_minimal() +
  scale_color_manual(values=mycols7, name= "Treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"), title = "Bray-Curtis PCoA Total DNA") +
  facet_wrap(~N, ncol=2)
  
dev.off()

### not with facet wrap
p1<-df.pcoa %>% filter(N=="0") %>%
ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal() +
  scale_color_manual(values=mycols7, name= "Treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Bray-Curtis PCoA  Total DNA ",
       subtitle = "A  -Nitrogen")+
  stat_ellipse(aes(group=Treatment), linetype=2)

p2<-df.pcoa %>% filter(N=="1") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal() +
  scale_color_manual(values=mycols7, name= "Treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Bray-Curtis PCoA  Total DNA ",
       subtitle = "B  +Nitrogen")+
  stat_ellipse(aes(group=Treatment), linetype=2)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="pcoa.total.Nitrogen.svg",width = 8, height=4)

require(gridExtra)
#windows(8,4)
grid.arrange(p1, p2, ncol=2)
dev.off()


### color by Nitrogen
svg("total.pcoa.Nonly.svg", height = 4, width =4)
df.pcoa %>% #filter(N=="0") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(N))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=mycols8, name= "Nitrogen") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Bray-Curtis PCoA  Total DNA ",
       subtitle = "")+
  stat_ellipse(aes(group=N), linetype=2)
dev.off()

#both treatments
svg("total.pcoa.full.svg", height = 4, width =4)
df.pcoa %>% 
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=mycols7, name= "Treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Bray-Curtis PCoA  Total DNA ",
       subtitle = "")+
  stat_ellipse(aes(group=Treatment), linetype=2)

dev.off()

# base R below just colors by N and no N
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$N   <- factor(metadat2$N)

svg("total.pcoa.Nonly.svg", height = 5, width =5)

ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Total DNA",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, 
         col= mycols8[metadat2$N],
         pch= myshapes[as.factor(metadat2$Rep)],
         lwd=1,cex=1.5,
         bg=mycols8[metadat2$N])
  
  legend("topleft", legend=c( 0, 1) ,
         fill= mycols8,
         cex=1,
         title = "Nitrogen",
         bty = "n")
  
  ordiellipse(otus.pcoa, metadat2$N,  
              kind = "ehull", conf=0.95, label=T, 
              draw = "polygon",
              border = 0,
              lwd=.1,
              col= mycols8,
              alpha = 60)
dev.off()



#####PLOTS  facet by number of species#####

p1<-df.pcoa %>% filter(N=="0") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols7, guide= "none") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% var. explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "PCoA  Total DNA ",
       subtitle = "A -Nitrogen")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)
  
  


p2<-df.pcoa %>% filter(N=="1") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal( base_size = 14) +
  scale_color_manual(values=mycols7, guide = "none") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% var. explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "PCoA  Total DNA ",
       subtitle = "B  +Nitrogen")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="pcoa.total.nspecies.svg",width = 8, height=6)

require(gridExtra)
windows(9,7)
grid.arrange(p1, p2, ncol=2)
dev.off()


#####PLOTS by group#####

mycols3<- c(  "#f4d35e", "#e94f37","#006d77")


p1<-df.pcoa %>% filter(N=="0") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Legume))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=c("grey",  "#006d77"), name= "Legume") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Total DNA -N",
       subtitle = "Legume")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)


p2<-df.pcoa %>% filter(N=="0") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Grass))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=c("#5e615f",  "#e6b802"), name= "Grass") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Total DNA -N ",
       subtitle = "Grass")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)


p3<-df.pcoa %>% filter(N=="0") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Brassicae))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=c("grey",  "#e94f37"), name= "Brassica") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Total DNA -N  ",
       subtitle = "Brassica")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="pcoa.total.functionalgroup.svg",width = 10, height=6)

#require(gridExtra)
#windows(12,8)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#### with nitrogen



p1<-df.pcoa %>% filter(N=="1") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Legume))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=c("grey",  "#006d77"), name= "Legume") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Total DNA +N",
       subtitle = "Legume")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)


p2<-df.pcoa %>% filter(N=="1") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Grass))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=c("#5e615f",  "#e6b802"), name= "Grass") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Total DNA +N ",
       subtitle = "Grass")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)


p3<-df.pcoa %>% filter(N=="1") %>%
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Brassicae))) +  
  geom_point(size = 3, alpha=.9) +
  theme_minimal() +
  scale_color_manual(values=c("grey",  "#e94f37"), name= "Brassica") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% variance explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "Total DNA +N  ",
       subtitle = "Brassica")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=1)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")
svg(file="pcoa.total.functionalgroupN.svg",width = 10, height=6)

#require(gridExtra)
#windows(12,8)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



#### total CAP #####
# Constrained ordination
# Perform vegdist analysis of BC distances #
total <-subset_samples(total, Treatment !="Soil" )
otus.bray<-vegdist(otu_table(total), method = "bray")
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")

#set factors 
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

p1.cap <- ordinate(total, method='CAP',distance='bray',formula=~N*Legume*Grass*Brassicae)
anova.cca(p1.cap, by="terms)


p1.cap <- ordinate(total, method='CAP',distance='bray',formula=~N*Treatment)
anova.cca(p1.cap, by="terms")

#mycols7
#B G GB L LB LG LBG
mycols7<-c( "#4B2D4BFF", "#AD5A6BFF", "#E3C1CBFF", "#365C83FF", "#384351FF", "#4D8F8BFF", "#CDD6ADFF")

# Visualize ordination
windows(6,6)

svg("cap.total.svg", width = 6, height = 4)
plot_ordination(total,  p1.cap,color="Treatment")+
  #facet_wrap(~n_species)+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
        legend.position="left")+
  scale_color_manual(values =  mycols7, name="Treatment")+
  scale_shape_discrete(name= "Nitrogen")
dev.off()                     

p1.cap$CCA
p1.cap$terms

#########cap by treatments subsetted#####

t1<-subset_samples(total, n_species=="1" & Legume=="0")
#asvs.clean<-otu_table(t1)
metadat2<-as.data.frame(as.matrix(sample_data(t1)))
head(metadat2)
p1.cap <- ordinate(t1, method='CAP',distance='bray',formula=~N*Brassicae)
anova.cca(p1.cap, by="terms")

####total permanova ######

asvs.clean<-otu_table(total)
metadat2<-as.data.frame(as.matrix(sample_data(total)))
head(metadat2)
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
head(asvs.bray)
asvs.perm<- adonis2(asvs.clean ~ N, data = metadat2, permutations = 999, method="bray")
asvs.perm

asvs.perm<- adonis2(asvs.clean ~ (Grass + Brassicae + Legume + N)^2, data = metadat2, permutations = 999, method="bray")
asvs.perm
adonis2(formula = asvs.clean ~ (Grass + Brassicae + Legume)^3, data = metadat2, permutations = 999, method = "bray")

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

###total Permanova  between treatments #######

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


##############PCOA of inactive vs active###########

# Calculate Bray-Curtis distance between samples
# Perform PCoA analysis of BC distances #
fc <-subset_samples(fc, Treatment !="Soil" )
otus.bray<-vegdist(otu_table(fc), method = "bray")
otus.pcoa <- cmdscale(otus.bray, k=(44-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
colnames(otus.p) <- c("PC1", "PC2")
df.pcoa <- cbind(sample_data(fc), otus.p)
# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2 <- as.data.frame(sample_data(fc))

#set factors 
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction <- factor(metadat2$Fraction)

svg("pcoa.total.svg", height = 5, width =5)
windows(6,6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Sorted Cells",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, 
       col= mycols8[metadat2$Fraction],
       pch= myshapes[as.factor(metadat2$Rep)],
       lwd=2,cex=2,
       bg=mycols8[metadat2$Fraction])

legend("topleft", legend=c( "Active", "Inactive") ,
       fill= mycols8,
       cex=1,
       title = "Treatment",
       bty = "n")

ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            lwd=.1,
            col= mycols7,
            alpha = 40)
dev.off()



#######PCOA plot active#########

act<-subset_samples(fc, Fraction=="Active")
act<-prune_taxa(taxa_sums(act) > 0, act)
metadat2<-as.data.frame(sample_data(act))

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(act), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
otus.p3 <- otus.pcoa$points[,3:4]
colnames(otus.p) <- c("PC1", "PC2")
colnames(otus.p3) <- c("PC3", "PC4")

df.pcoa <- cbind(sample_data(act), otus.p)
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

# subset metadata

#faction
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction   <- factor(metadat2$Fraction)

myshapes <- c(1, 12, 15 ,21, 22, 23 , 24, 25)

#windows(6,6)
svg("pcoa.3.svg",  width = 6, height = 6 )
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="active ",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, 
       col= mycols8[metadat2$Treatment],
       pch= myshapes[metadat2$Rep],
       lwd=1,cex=1.5,     
       bg=mycols8[metadat2$Treatment])

ordiellipse(otus.pcoa, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= mycols8,
            alpha = 30)
dev.off()

### pc 3 and 4 
ordiplot(otus.pcoa,choices=c(3,4), type="none", main="active  ",xlab=paste("PCoA3 (",round(pe3,2),"% variance explained)"),
         ylab=paste("PCoA4 (",round(pe4,2),"% variance explained)"))
points(otus.p, 
       col= mycols8[metadat2$Treatment],
       pch= myshapes[metadat2$Treatment],
       lwd=1,cex=1.5,     
       bg=mycols8[metadat2$Treatment])

ordiellipse(otus.pcoa, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= mycols8,
            alpha = 30)
dev.off()


####faceted plot ######
svg("active.nspecies.svg", width=7, height=4)
#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols7, name="treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% var. explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "PCoA Active ",
       subtitle = "A")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=3)

dev.off()

svg("active.nspeciesP34.svg", width=7, height=4)
#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC3, y = PC4, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols7, name="treatment") +
  labs(x = paste("PCoA3 (",round(pe3,2),"% var. explained)"), y = paste("PCoA4 (",round(pe4,2),"% variance explained)"),
       title = "PCoA Active ",
       subtitle = "A")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=3)

dev.off()



####PCOA inactive######

iact<-subset_samples(fc, Fraction=="Inactive")
iact<-prune_taxa(taxa_sums(iact) > 0, iact)
metadat2<-as.data.frame(sample_data(iact))

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(iact), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
otus.p <- otus.pcoa$points[,1:2]
colnames(otus.p) <- c("PC1", "PC2")
df.pcoa <- cbind(sample_data(iact), otus.p)

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

# subset metadata

#faction
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat2$Fraction   <- factor(metadat2$Fraction)

myshapes <- c(1, 12, 15 ,21, 22, 23 , 24, 25)

#windows(6,6)
svg("pcoa.3.svg",  width = 6, height = 6 )
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Inactive ",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, 
       col= mycols8[metadat2$Treatment],
       pch= myshapes[metadat2$Rep],
       lwd=1,cex=1.5,     
       bg=mycols8[metadat2$Treatment])

ordiellipse(otus.pcoa, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= mycols8,
            alpha = 30)
dev.off()

### pc 3 and 4 
ordiplot(otus.pcoa,choices=c(3,4), type="none", main="Inactive  ",xlab=paste("PCoA3 (",round(pe3,2),"% variance explained)"),
         ylab=paste("PCoA4 (",round(pe4,2),"% variance explained)"))
points(otus.p, 
       col= mycols8[metadat2$Treatment],
       pch= myshapes[metadat2$Treatment],
       lwd=1,cex=1.5,     
       bg=mycols8[metadat2$Treatment])

ordiellipse(otus.pcoa, metadat2$Treatment,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= mycols8,
            alpha = 30)
dev.off()


####faceted plot ######
svg("Inactive.nspecies.svg", width=8, height=4)
#windows(4,8)
df.pcoa %>% 
  ggplot( aes(x = PC1, y = PC2, color= as.factor(Treatment))) +  
  geom_point(size = 3, alpha=.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values=mycols7, name="treatment") +
  labs(x = paste("PCoA1 (",round(pe1,2),"% var. explained)"), y = paste("PCoA2 (",round(pe2,2),"% variance explained)"),
       title = "PCoA Inactive ",
       subtitle = "A")+
  stat_ellipse(aes(group=Treatment), linetype=2)+
  facet_wrap(~n_species, ncol=3)

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

########PERMANOVA##############

## total 
#full dna  between trts
ps2<-subset_samples(ps, Fraction =="Total"& Treatment!="Soil")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
ps2

#get asvs table
asvs.clean<-otu_table(ps2)
metadat2<-filter(metadat, Fraction =="Total" & Treatment!="Soil")

# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
#
asvs.perm<- adonis2(asvs.clean ~ Treatment, data = metadat2, permutations = 999, method="bray")
asvs.perm
#adonis2(formula = asvs.clean ~ Treatment, data = metadat2, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Treatment  6   3.0996 0.13358 1.9528  0.001 ***
##  Residual  76  20.1049 0.86642                  
##Total     82  23.2045 1.00000  

asvs.perm<- adonis2(asvs.clean ~ (Grass + Brassicae + Legume)^2, data = metadat2, permutations = 999, method="bray")
asvs.perm
adonis2(formula = asvs.clean ~ (Grass + Brassicae + Legume)^3, data = metadat2, permutations = 999, method = "bray")

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



########PERMANOVA active ########### 
#full dna  between trts
ps2<-subset_samples(ps, Fraction =="Active"& Treatment!="Soil")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
ps2

#get asvs table
asvs.clean<-otu_table(ps2)
metadat2<-filter(metadat, Fraction =="Active"& Treatment!="Soil")

# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(asvs.clean, method = "bray")
#
asvs.perm<- adonis2(asvs.clean ~ Treatment, data = metadat2, permutations = 999, method="bray")
asvs.perm
#adonis2(formula = asvs.clean ~ Treatment, data = metadat2, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Treatment  6   3.0996 0.13358 1.9528  0.001 ***
##  Residual  76  20.1049 0.86642                  
##Total     82  23.2045 1.00000  

asvs.perm<- adonis2(asvs.clean ~ (Grass + Brassicae + Legume)^2, data = metadat2, permutations = 999, method="bray")
asvs.perm
adonis2(formula = asvs.clean ~ (Grass + Brassicae + Legume)^3, data = metadat2, permutations = 999, method = "bray")

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







#####BARPLOT TOP ASVS FIG 5######
# this is on rarefied ASVS
# filter for total and active

df<-subset_samples(ps.r, Compartment!="ctl"& Compartment!="Bulksoil" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df<-as.data.frame(t(otu_table(df)))
df<-cbind(taxon,df)

#get means
df$nodule.total.mean <-   rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
df$nodule.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root.total.mean <-   rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
df$root.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
colnames(df)
#remove extra columns
asv <- dplyr::select(df, -contains("POS") )  %>% dplyr::select(., -contains("SYBR"))
head(asv)

# select top 50 ASVs viable rhizosphere
asv<-asv[order(asv$rhizo.total.mean, decreasing = TRUE),]
top<-asv[1:50,]
top$asv<-row.names(top)
top

# what percent of the whole community are the top asvs
all<-colSums(asv[,12:13])
t<-colSums(top[,12:13])
t/all
# top 50 otus is 27% of the total population;

#edit labels
top$Phyla<-sub("p__", "", top$Phyla)
top<-top[order(top$rhizo.total.mean, decreasing = TRUE),]
top$otu1<-c(1:50)
top$otu1<-paste0("ASVS ",top$otu1)
top

#colors
# acido - dk blue
# actino - light blue
# bacteriodota - pink
# chlorofexi light red
# proteobacteria, light green
# Methylomirabilota
# verrucomicrobiota - pale gold


mycols8<- c( "#1F78B4","#A6CEE3", "#75026d",  "#ed6361",  "#6A3D9A", "#B2DF8A", "#FF7F00","#FDBF6F")

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/")

svg(file="Fig5_scatter/barplot_toptaxa.svg",width = 8, height=5)
windows(6,5)
ggplot(top)+
  geom_bar(aes(x=reorder(otu1, -rhizo.total.mean) , y=rhizo.bcat.mean, fill= Phyla), stat="identity", position="identity") +
  geom_line(aes(x=as.numeric(reorder(otu1, -rhizo.total.mean)) , y=rhizo.total.mean), linewidth=1) +
  xlab("top 50 Viable Cell ASVS")+
  ylab("Average rarefied reads per sample")+
  scale_fill_manual(values=mycols8)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_blank(), legend.position = c(.8, .7))+
  guides(fill=guide_legend(title="Abundance Active Cells"))
dev.off()


## with labs
svg(file="Fig5_scatter/barplot_toptaxa_lab.svg",width = 8, height=5)
ggplot(top)+
  geom_bar(aes(x=reorder(otu1, -rhizo.total.mean) , y=rhizo.bcat.mean, fill= Phyla), stat="identity", position="identity") +
  geom_line(aes(x=as.numeric(reorder(otu1, -rhizo.total.mean)) , y=rhizo.total.mean), linewidth=1) +
  xlab("top 50 Viable Cell ASVS")+
  ylab("Average rarefied reads per sample")+
  scale_fill_manual(values=mycols8)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_text(angle=90), legend.position = c(.8, .7))+
  guides(fill=guide_legend(title="Abundance Active Cells"))
dev.off()
    
#####BARPLOT PHYLA LEVEL FIG 6 #######
#use non rarefied taxa because we are analyzing the proportion of read we will normalize by Number of reads.
# filter for total and active rhizosphere
    ps
    sample_data(ps)
    df<-subset_samples(ps, Compartment!="ctl"& Compartment=="Rhizosphere" & Fraction!="Total_DNA")
    df<-prune_taxa(taxa_sums(df) > 0, df)
    taxon<-tax_table(df)
    df1<- as.data.frame(otu_table(df)) # this is for later :)
    df<-as.data.frame(otu_table(df))
    dim(df)
    
# normalize by number of reads
    df<-df/rowSums(df)
    df<-as.data.frame(t(df))
    head(df)
    # 4854 taxa
    df<-cbind(taxon,df)
 
    #aggregate
    df<-aggregate(cbind(C10R.POS_S30, C10R.SYBR_S20, C1R.POS_S27,   C1R.SYBR_S16, C2R.POS_S28 , 
           C2R.SYBR_S17 , C5R.POS_S29 ,  C5R.SYBR_S18 , C7R.SYBR_S19  ) ~ Phyla, data = df, FUN = sum, na.rm = TRUE)
    colnames(df)<- c("Phyla", "BONCAT_1" , "Total_1", "BONCAT_2",   "Total_2" , "BONCAT_3"  , "Total_3" ,
                    "BONCAT_4"  , "Total_4",  "Total_5") 
  
    head(df)
    dim(df)
    
# grab top phyla
    row.names(df) <- df$Phyla
    df$sum <- rowSums(df[,2:9])
    df<-df[order(-df$sum),]    
    top<-df[c(1:5,7),]
  # rm sum
  top
  top<-top%>%
  pivot_longer(2:10, values_to = "abundance", names_to= "rep_fraction" )
  top<-top%>%mutate(fraction = str_split_i(top$rep_fraction, "_", 1)) %>%
  mutate(Phyla= str_split_i(top$Phyla, "_", 3))
  head(top)
  
#label
top<-top %>% mutate(fraction=recode(fraction, 'Total'='Viable Cells'))
top$fraction<-factor(top$fraction, levels = c('Viable Cells', 'BONCAT'))
#colors
# acido - dk blue
# actino - light blue
# bacteriodota - pink
# chlorofexi light red
# proteobacteria, light green
# Methylomirabilota
# verrucomicrobiota - pale gold


mycols8<- c( "#1F78B4","#A6CEE3", "#75026d",  "#ed6361",  "#6A3D9A", "#B2DF8A", "#FF7F00","#FDBF6F")

phycols7<-c("#1F78B4","#A6CEE3", "#75026d",  "#FB9A99", "#33A02C","#FF7F00",  "#FB9A99")

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")
svg(file="barplot_phyla.svg",width = 8, height=4)
windows(7,3)
ggplot(top) +
  geom_boxplot(aes(x= reorder(fraction, -sum), y= abundance, fill=Phyla), outlier.shape = NA, alpha=.5, show.legend = FALSE)+
  geom_jitter(aes(x= reorder(fraction, -sum), y= abundance, fill=Phyla), position = position_jitter(width = .2), show.legend = FALSE)+
  facet_wrap(~reorder(Phyla, -sum), nrow=1)+
  theme_minimal(base_size = 12)+
  scale_fill_manual(values=phycols7)+
  labs(y="proportion of reads in sample")+
  theme(axis.text.x=element_text(angle=45, hjust=0.9), axis.title.x = element_blank())
dev.off()
#####BARPLOT STATS ########
# mixed model with binomial regression
# we need a data frame of # success to # failure


df<-subset_samples(ps, Compartment!="ctl"& Compartment=="Rhizosphere" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df1<- as.data.frame(otu_table(df)) # this is for later :)
# make column that is the number of reads
reads<-rowSums(df1)
df1<-as.data.frame(t(df1))
# add taxa info
df1<-cbind(taxon,df1)
#aggregate
df1<-aggregate(cbind(C10R.POS_S30, C10R.SYBR_S20, C1R.POS_S27,   C1R.SYBR_S16, C2R.POS_S28 , 
                    C2R.SYBR_S17 , C5R.POS_S29 ,  C5R.SYBR_S18 , C7R.SYBR_S19  ) ~ Phyla, data = df1, FUN = sum, na.rm = TRUE)

# grab top phyla
row.names(df1) <- df1$tPhyla
df1$sum <- rowSums(df1[,2:9])
df1<-df1[order(-df1$sum),]    
top<-df1[c(1:5,6,7),]
top
# pivot
top<-top%>%
  pivot_longer(2:10, values_to = "abundance", names_to= "rep_fraction" )
top<-top%>%mutate(fraction = str_split_i(top$rep_fraction, "_", 1)) %>%
  mutate(Phyla= str_split_i(top$Phyla, "_", 3))

#put No. reads in
top$reads =rep(reads, 7 )
top<-select(top , -sum, -fraction, -rep_fraction)
top$failures = top$reads-top$abundance
trt<-c(rep(c("Active", "Total"), 4), "Total")
top$trt = rep(trt, 7)
top$rep<-rep(c(1,1, 2, 2, 3, 3,4, 4,5),7)
top


#proteobacteria model
df<-filter(top, Phyla=="Proteobacteria") %>% select(abundance,failures, trt, rep )

m1<-glm(data= df, cbind(abundance,failures)~trt+rep, family = binomial)
anova(m1, test= "LRT")  
#
#acido model
df<-filter(top, Phyla=="Acidobacteriota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df, cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  


#Verrucomicrobiota model
df<-filter(top, Phyla=="Verrucomicrobiota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  
#not sig

# Actinobacteriota model
df<-filter(top, Phyla=="Actinobacteriota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  


#Bacteroidota model
df<-filter(top, Phyla=="Bacteroidota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  
#not sig 

#Planctomycetota model
df<-filter(top, Phyla=="Planctomycetota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  

#Chloroflexota model
df<-filter(top, Phyla=="Chloroflexota") %>%select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
anova(m1, test= "LRT")  



##### ANCOM ###############
# ancom is run on unrarefied data
# subset to rhizosphere active and viable
ps1<-subset_samples(ps, Compartment=="Rhizosphere" & BONCAT!="DNA")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# add 1 to everything... b/c ancom can't deal with structural zeros.
df<-as.data.frame(otu_table(ps1))
df<-df+1
# remove other column
df<-select(df,-Others)

# add asv column
taxon<-as.data.frame(tax_table(ps1))
taxon$asv<-row.names(taxon)
#remove others
taxon<-filter(taxon, asv!="Others")

Workshop_OTU <- otu_table(as.matrix(df), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(ps1)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps1 <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps1

### skip to data import if you already ran ancom
sample_data(ps1)

ps2 = mia::makeTreeSummarizedExperimentFromPhyloseq(ps1)

ps2

out1 = ancombc(data = ps2, assay_name = "counts", 
               tax_level = "asv", phyloseq = NULL, 
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
               group = "Fraction", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
sample_data(ps1)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/tables")
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
col_name = c("asv", "DA_Intercept", "DA_Fraction_Total_cells_Active")
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

######import df from ANCOM######
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/tables")
tab<-read.table("ancom_table.txt", header = TRUE)
head(tab)

##add add abundance and taxon info
df<-as.data.frame(t(otu_table(ps1)))
df$asv<-row.names(df)
df<-left_join(df, tab)
taxon<-as.data.frame(tax_table(ps1))
df<-left_join(taxon, df)
df

### summarise the abundance in active and total
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
t<-df %>% select(contains("R.SYBR"))
sd_total<- apply(t, 1, sd, na.rm=TRUE)
sd_total
df$sd_total <- sd_total
  
df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
t<-df %>% select(contains("R.POS"))
sd_active<- apply(t, 1, sd, na.rm=TRUE)
sd_active
df$sd_active <- sd_active
head(df)

#remove some columns to make things simpler
df<-df %>% select(-contains("R."))
df<-df %>% select(-contains("Intercept"))
df<-df %>% select(-contains("DNA"))
colnames(df)

#quick histogream
hist(df$rhizo.bcat.mean, breaks = 100)
hist(df$rhizo.total.mean, breaks=100)

# make a column the lowest taxa rank assigned 
df$label<-df$Species
df$label[which(df$label==" s__"  )] <- df$Genus[which(df$label==" s__" )]
df$label[which(df$label==""  )] <- df$Genus[which(df$label=="" )]
df$label[which(df$label==" g__"  )] <- df$Family[which(df$label==" g__" )]
df$label[which(df$label==""  )] <- df$Family[which(df$label=="" )]
df$label[which(df$label==" f__"  )] <- df$Order[which(df$label==" f__" )]
df$label[which(df$label== " g__SCN-69-37" )] <- df$Order[which(df$label==" g__SCN-69-37" )]
df$label[which(df$label== " f__UBA2999" )] <- df$Order[which(df$label== " f__UBA2999"  )]
df$label[which(df$label== " g__WHSN01" )] <- df$Order[which(df$label== " g__WHSN01"  )]
df$label[which(df$label== " s__UBA11740 sp003168335"    )] <- df$Order[which(df$label== " s__UBA11740 sp003168335"    )]
df$label[which(df$label== " g__PSRF01"      )] <- df$Order[which(df$label== " g__PSRF01"    )]
df$label[which(df$label==    " s__VFJQ01 sp009885995"    )] <- df$Family[which(df$label==  " s__VFJQ01 sp009885995"  )]
df$label[which(df$label==   " s__OLB17 sp001567505"       )] <- df$Family[which(df$label==  " s__OLB17 sp001567505"     )]

# right now the active community is the reference. It's a little confusing. 
# multiply by -1 to make it so the viable community is the reference! 
# so a negative number would be depleted in active
df$LFC_FractionViable_Cell<-df$LFC_FractionViable_Cell*-1

#rm "p__" in phyla
df$Phyla<-sub(" p__", "", df$Phyla)
# number ASvs
df$asv_no <- paste0("ASV", row.names(df))

# quick overall volcano plot to check distribution
ggplot(df, aes(x=LFC_FractionViable_Cell , y=p_viable)) + 
  geom_jitter()+ 
  theme_bw()
# lfc verse the pvalue
ggplot(df, aes(x=LFC_FractionViable_Cell , y=-log(p_viable), col=DA_Fraction_Total_cells_Active)) + 
  geom_jitter()+ 
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  theme_bw( )

#LOOK AT DA taxa
DA<-df %>% filter(DA_Fraction_Total_cells_Active==TRUE)
dim(DA)
head(DA)
#save DA taxa file
write.csv(DA, file="ancom_DA_taxa.csv")

# sort for negative slope - more dormant taxa
dormant<-DA[DA$LFC_FractionViable_Cell <0,]
dormant<-dormant[order(dormant$rhizo.total.mean, decreasing = TRUE),]
dormant
write.csv(dormant, file ="ancom_dormanttaxa.csv")

# sort for more active taxa - postive slope
active<-DA[DA$LFC_FractionViable_Cell > 0,]
active<-active[order(active$rhizo.bcat.mean, decreasing = TRUE),]
head(active)
write.csv(active, file ="ancom_activetaxa.csv")

######DA TAXA figure ############ 
# filter for at least 50 reads
DA<-filter(DA, DA$rhizo.bcat.mean>50 |  DA$rhizo.total.mean>50  )
dim(DA)

#remove numbers as the end of labels to make them easier to read in the figure
DA$label<-gsub("_48326", "", DA$label)
DA$label<-gsub("_48670", "", DA$label)
DA$label<-gsub("_A_50105", "", DA$label)

#add asvs No. label
DA$label<-paste(DA$asv_no, DA$label)
DA$label

DA$type <- "Differentially Abundant ASVS"

#add * for DA taxa
DA$sig<-ifelse( DA$DA_Fraction_Total_cells_Active=="TRUE", "*", "")
DA$label<-paste(DA$label, DA$sig)
DA$label

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")

### cols for this plot
# actinobacteria - 
# acidobacteria - dk blue 
# actino - light blue
# proteobaceria - light green
# gemma - 


svg(file="barplot_top20_DA.svg",width = 8, height=7.5)
windows(8,7.5)
ggplot(DA)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla),
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=LFC_FractionViable_Cell-se_viable,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values=c("#A6CEE3", "#75026d","#B2DF8A"))+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab("Differentially Abundant Asvs")+
  ylab("Log fold change viable to active")
  #theme(legend.position = "none")
  
dev.off()  


###LFC top 10 taxa ##
## filter for top 50 most abundant taxa or DA taxa
df<-df %>%
  arrange( -rhizo.total.mean)
top<-df[1:10,]

#remove the " f__ " part
top$label <-substr(top$label, 5, nchar(top$label)-1)

#remove numbers as the end of labels to make them easier to read in the figure
top$label<-gsub("_48326", "", top$label)
top$label<-gsub("_48670", "", top$label)
top$label<-gsub("_A_50105", "", top$label)
top$label<-gsub("_E_64746", "", top$label)
top$label<-gsub("_58024", "", top$label)
top$label<-gsub("_A_58049", "", top$label)
top$label

#add asvs label
top$label<-paste(top$asv_no, top$label)

#add * for DA taxa
top$sig<-ifelse( top$DA_Fraction_Total_cells_Active=="TRUE", "*", "")
top$label<-paste(top$label, top$sig)
top$label

#group label
top$type <- "Most Abundant Asvs"
### cols for this plots
# acidobacteria - dk blue 
# proteobaceria - light green
head(top)
#plot

svg(file="barplot_top10.svg",width = 8, height=3.5)
windows(8,8)
ggplot(top)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla), 
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=LFC_FractionViable_Cell-se_viable,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values=c("#1F78B4", "#B2DF8A"))+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab("Most abudant Asvs")+
  ylab("Log fold change viable to active")+
  # facet_grid( scales = "free", space = "free",  rows=vars(Phyla)) 
  theme(legend.position = "none")
dev.off()  

##DA fig combined###
#combine
combine<-full_join(top, DA)
head(combine)

#colors
"#06568c" #  acidobacteria - dk blue 
"#B2DF8A" # proteobaceria - light green
"#52b8d1" # actinobacteria - light blue
"#d40d63" # Gemmatimonadota raspberry

mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A" )
#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")
svg(file="barplot_big.svg",width = 8, height=7)
windows(8,8)
ggplot(combine)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla), 
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values= mycols)+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab(" Asvs")+
  ylab("Log fold change viable to active")
  #facet_grid( space = "free",  rows=vars(Type))
  #theme(legend.position = "none")
dev.off()  

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")
svg(file="barplot_big_label.svg",width = 8, height=7)
windows(8,8)
ggplot(combine)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla), 
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values= mycols)+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab(" Asvs")+
  ylab("Log fold change viable to active")+
  #facet_grid( space = "free",  rows=vars(Type))+
  theme(legend.position = "none")
dev.off()  

#####VENN DIAGRAM Viable#####
#on rarefied data
sample_data(ps.r)
ps1<-subset_samples(ps.r, Fraction=="Viable_Cell")
# at least in 3 samples min reads is 50
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 542 taxa

# grab taxonomy + asv column
taxon<-as.data.frame(tax_table(ps1))
taxon$asv<-row.names(taxon)
taxon<-filter(taxon, asv!="Others")

#grab data and rename rows by compartment
df<-as.data.frame(otu_table(ps1))

head(df)
n<-row.names(df)
n[grepl("N" , n)]="nodule"
n[grepl("E" , n)]="roots"
n[grepl("R" , n)]="rhizo"
n

# sum by compartment
df<-rowsum(df, n)
#transform
df<-as.data.frame(t(df))
#remove others row
df<-df[-which(row.names(df)=="Others"),]
head(df)

# present = 1
df[df>1] <- 1

# make 3 groups
nodule<-rownames(df[df$nodule==1,])
roots<-rownames(df[df$roots==1,])
rhizo<-rownames(df[df$rhizo==1,])

x <- list(
  nodule = nodule, 
  roots = roots, 
  rhizo = rhizo
)


#### venn diagram #
library(ggvenn)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig7_heatmap")
svg(file="OTU_level_total_venn.svg",width = 4, height=4 )
windows(4,4)
mycols= c( "#a3c9fa", "#4e97ed", "#045bc2")
ggvenn(
  x, 
  fill_color = mycols,
  stroke_size = 2, set_name_size = 4, text_size = 3, digits = 1, fill_alpha=.6
  #auto_scale = TRUE
) +
  ggtitle("Viable Cells")
dev.off()



#####VENN DIAGRAM Active######
### need 50 reads + in 3 samples
sample_data(ps)
ps1<-subset_samples(ps, Fraction=="Active_Cell")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 324 taxa

# grab taxonomy + asv column
taxon<-as.data.frame(tax_table(ps1))
taxon$asv<-row.names(taxon)
taxon<-filter(taxon, asv!="Others")

#grab data and rename rows by compartment
df<-as.data.frame((otu_table(ps1)))
head(df)
n<-row.names(df)
n[grepl("N" , n)]="nodule"
n[grepl("E" , n)]="roots"
n[grepl("R" , n)]="rhizo"
n

# summ by compartment
df<-rowsum(df, n)
df<-as.data.frame(t(df))
head(df)
#remove others
df<-df[-which(row.names(df)=="Others"),]

#make present =1 and split into groups
df[df>1] <- 1
head(df)
nodule<-rownames(df[df$nodule==1,])
roots<-rownames(df[df$roots==1,])
rhizo<-rownames(df[df$rhizo==1,])
df
x <- list(
  nodule = nodule, 
  roots = roots, 
  rhizo = rhizo
)
x

#######venn diagram #
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig7_heatmap")
library(ggvenn)

svg(file="active_venn.svg",width = 4, height=4 )

mycols = c( "#e89c6f", "#d1663f", "#992600")

windows(4,4)
ggvenn(
  x, 
  fill_color = mycols,
  stroke_size = 2, set_name_size = 4, text_size = 3, digits = 1, fill_alpha=.7
  #auto_scale = TRUE
 ) +
  ggtitle("Active Otus")
dev.off()



#####HEATMAP#########
##big heatmap, no phylogeny ##
# df of values in total and active
# at least in 3 samples min reads is 50
sample_data(ps.r)
df<-subset_samples(ps.r, Fraction=="Active_Cell" | Fraction=="Viable_Cell")
df<-ps_prune(df, min.samples = 3, min.reads = 50)
df # 1018 taxa
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
head(df)

####### agregate to the family level
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)
head(df)

#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"
df$Family[grepl(" f__Bacillaceae_H", df$Family)] <- " Bacillaceae_H"
df$Family[grepl(" f__Streptomycetaceae", df$Family)] <-  " f__Streptomycetaceae"

#aggregate
df<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,  C1R.SYBR_S16,  C2E.SYBR_S22,
                    C2N.SYBR_S15,  C2R.SYBR_S17 ,  C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25,  C7R.SYBR_S19,
                    C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32, C2N.POS_S62,  C5E.POS_S33,
                    C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 156 families

### summarize by compartment 
df$Viable_rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB"))) %>%   glimpse()
df$Viable_root <-   rowMeans(df %>% dplyr::select(contains("E.SYB"))) %>%   glimpse()
df$Viable_nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB"))) %>%   glimpse()

#### summarize by compartment 
df$Active_rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$Active_root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df$Active_nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()

df<-df%>% select(c(Active_rhizo, Active_root, Active_nodule, Viable_rhizo,Viable_root, Viable_nodule )) %>% mutate(otu= row.names(df))
head(df)

#rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>%   glimpse()
df<-mutate(df, Family=otu) %>% select(., -otu)
head(df)
#rm families that are unknown

df<-df[ !grepl(" f__UBA", df$Family) , ]
df<-df[ !grepl(" f__SG8", df$Family) , ]
df<-df[ !grepl(" f__SCT", df$Family) , ]
df<-df[ !grepl(" f__AC-14", df$Family) , ]
df<-df[ !grepl(" f__SCN", df$Family) , ]
df<-df[ !grepl(" f__RSA", df$Family) , ]
df<-df[ !grepl(" f__WHT", df$Family) , ]
df<-df[ !grepl(" f__RBG", df$Family) , ]
df<-df[ !grepl(" f__CSP", df$Family) , ]
df<-df[ !grepl(" f__DSM", df$Family) , ]
df<-df[ !grepl(" f__TK", df$Family) , ]
df<-df[ !grepl(" f__QHB", df$Family) , ]
df<-df[ !grepl(" f__2013", df$Family) , ]
df<-df[ !grepl(" f__B-17", df$Family) , ]
df<-df[ !grepl(" f__JA", df$Family) , ]
df<-df[ !grepl(" f__J0", df$Family) , ]
df<-df[ !grepl(" f__Gp", df$Family) , ]
df<-df[ !grepl(" f__GWC", df$Family) , ]
df<-df[ !grepl(" f__Fen", df$Family) , ]
df<-df[ !grepl(" f__FW", df$Family) , ]
df<-df[ !grepl(" f__HR", df$Family) , ]



m <- df
# remove those "F__
m$Family<-sub(" f__", "", m$Family)

row.names(m) <- m$Family

#rm family column
m<-select(m, -Family)

#make matrix

m<-log10(m)
m[m== "-Inf"] <- 0
m
#m<-m[,c(3,2,1)]
m<-as.matrix(m)
m

# clustering agrorithm
distance = dist(m, method = "euclidean")
distance
cluster = hclust(distance, method = "ward.D2")
my_list<-cluster$labels
#### use this order for heatmaps with out phylogeny
m<-m[match(row.names(m),as.character(my_list)),]

#p<-ggtree(cluster) + 
#  geom_tiplab(size=3, align=TRUE, linesize=0, offset = -.2) + 
#  theme_tree2()+
#  xlim_tree(1) 

#gheatmap(p, m, 
#         colnames=FALSE,
#         legend_title="active taxa", offset = .5) 
  #scale_x_ggtree() + 
  #scale_fill_gradient(low= "#fffaa2", high =  "#bb0000", aesthetics = "fill", na.value = "white",
  #                    name="Abundance in Active")+

   #ggtitle("Families in Active Community")
#dev.off()


# creates a own color palette
my_palette <- colorRampPalette(c("white", "#95cefc", "#04063b"))(n = 99)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.01,length=2),  # for red
               seq(0.1,0.8,length=48),           # for yellow
               seq(0.81,5.1,length=50))             # for green

windows(7,7)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")

svg(filename = "big_blue_heatmap.svg", width = 9, height = 12)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = 1,
         lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
         lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
         lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
        #key.title = "log 10 Abundance")
        #key.xlab = "")
         key = FALSE) 
#title("Active", line= -4)
dev.off()  
#key
svg(filename = "key.svg", width = 10, height = 10)

windows(10,10)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          #labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = .00001,
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
          lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
          key.title = "log 10 Abundance")
dev.off()  
#key = FALSE

#### hehe make a red one too
# creates a own color palette
my_palette <- colorRampPalette(c("white", "#fffaa2", "#bb0000"))(n = 99)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.01,length=2),  # for red
               seq(0.1,0.8,length=48),           # for yellow
               seq(0.81,5.1,length=50))             # for green

windows(7,7)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")

svg(filename = "big_red_heatmap.svg", width = 9, height = 12)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = 1,
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
          lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
          #key.title = "log 10 Abundance")
          #key.xlab = "")
          key = FALSE) 
#title("Active", line= -4)
dev.off()  

#key
svg(filename = "red_key.svg", width = 8, height = 8)

windows(10,10)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          #labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = .00001,
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
          lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
          key.title = "log 10 Abundance")
dev.off()  
          #key = FALSE) 

#####HEATMAP with phylogeny:SUPPLEMENT######
##select active taxa 
# at least in 3 samples min reads is 50
df<-subset_samples(ps.r, Fraction=="Active_Cell")
df<-ps_prune(df, min.samples = 3, min.reads = 50)
df # 510 taxa
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
#remove "other" ASVS where rare taxa went
df<-filter(df, !grepl("Others", row.names(df)))
# remove "other" row where rare taxa went
taxon<-taxon[c(which(row.names(taxon)!="Others")),]

####### agregate to the family level
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)
head(df)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"
df$Family[grepl(" f__Blastocatellaceae", df$Family)] <- " f__Blastocatellaceae"

#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#aggregate
df<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 86 families

#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)

#rm unknown families
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")
df<-filter(df, !grepl("UBA", df$otu))
df<-filter(df, !grepl("SG", df$otu))
df<-filter(df, !grepl("SCN", df$otu))

# rename column
df<-mutate(df, Family=otu) %>% select(., -otu)

#removed 'f__'
df$Family<-sub(" f__", "", df$Family)
row.names(df) <- df$Family

#importtree
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")

#match tree to my taxa
length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
mynames<-(setdiff(unique(df$Family), tree$tip.label))
length(mynames)
mynames

#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, df$Family) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE,  cex = .5)
# grab correct order 
target<-tree.short$tip.label
df<-df[match(target, df$Family),]
head(df)
dim(df)
tree.short

# rm family column from df
df<-subset(df, select = c( -Family))
dim(df) # 46 families
# make matrix

m <- df
summary(m)
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")
svg(file="heatmap_phyl_active.svg",width = 6, height=6)
windows(6,6)

p<- ggtree(tree.short, branch.length = .001) + 
 geom_tiplab(size=4, align=TRUE, linesize=0, offset = 0) + 
 theme_tree2()+
 xlim_tree(1) 
gheatmap(p, m, 
         colnames=FALSE,
         legend_title="active taxa", offset = 1.5, font.size = 10,) +
  scale_x_ggtree() + 
  scale_fill_gradient(low= "#fffaa2", high =  "#bb0000", aesthetics = "fill", na.value = "white",
                      name="Abundance in Active")+
  ggtitle("Families in Active Community")
dev.off()


##viable###
#select viable cell communtiy, ASVs with at least in 3 samples min reads is 50
df<-subset_samples(ps, Fraction=="Viable_Cell")
df<-ps_prune(df, min.samples = 3, min.reads = 50)
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
#remove "other" ASVS where rare taxa went
df<-filter(df, !grepl("Others", row.names(df)))
# remove "other" row where rare taxa went
taxon<-taxon[c(which(row.names(taxon)!="Others")),]


####### aggregate to the family level
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"
df$Family[grepl(" f__Blastocatellaceae", df$Family)] <- " f__Blastocatellaceae"

#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

colnames(df)
head(df)
#aggregate
df<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,
                    C1R.SYBR_S16,  C2E.SYBR_S22, C2N.SYBR_S15,  C2R.SYBR_S17,
                    C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25, C7R.SYBR_S19)
                     ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) #103
head(df)

#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.SYB"))) %>%   glimpse()
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)

#rm unknown families
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")
df<-filter(df, !grepl("UBA", df$otu))
df<-filter(df, !grepl("SG", df$otu))
df<-filter(df, !grepl("SCN", df$otu))
df$otu

# rename column
df<-mutate(df, Family=otu) %>% select(., -otu)

#removed 'f__'
df$Family<-sub(" f__", "", df$Family)
row.names(df) <- df$Family


#importtree
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")
tree
# 1057 tips
#length(tree$tip.label) # look at the tip labels 
# modify tip labels
#tree$tip.label <- paste0(" f__", tree$tip.label)
length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
mynames<-(setdiff(unique(df$Family), tree$tip.label))
length(mynames)
mynames

#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, df$Family) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE,  cex = .5)
# grab correct order 
target<-tree.short$tip.label
df<-df[match(target, df$Family),]
head(df)
dim(df)
tree.short

# rm family column from df
df<-subset(df, select = c( -Family))
dim(df) # 49 families
# make matrix

m <- df
summary(m)
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")
svg(file="heatmap_phy_viable.svg",width = 6, height=6)
windows(6,6)

p<- ggtree(tree.short, branch.length = .001) + 
  geom_tiplab(size=4, align=TRUE, linesize=0, offset = 0) + 
  theme_tree2()+
  xlim_tree(1) 
gheatmap(p, m, 
         colnames=FALSE,
         legend_title="Viable taxa", offset = 1.5, font.size = 10) +
  scale_x_ggtree() + 
  scale_fill_gradient(low= "#96bceb", high =  "#033d85", aesthetics = "fill", na.value = "white",
                      name="Abundance in Viable")+
  ggtitle("Families in Viable Community")
dev.off()




###BINOMIAL model####
# rrarefied data to be able to compare samples 
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/")
### Import Data ###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)

## Transpose ASVS table ##
asvs.t <- t(asvs.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs.t<-asvs.t[order(row.names(asvs.t)),]

## Determine minimum available reads per sample ##
min.s<-min(rowSums(asvs.t))
min.s
### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs.t, min.s)
dim(asvs.t)
dim(asvs.r)

###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Viable_Cell', 'POS'='Active_Cell', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Fraction, metadat$Compartment))

##---make phyloseq object with rarefied data -------#
asvs.phyloseq<- (asvs.r)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_ASVS <- otu_table(asvs.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

# taxa that are in the plant and are unassigned
#afb244d96b70a4c948b205f7f7eea5c5 259366
#9bf55fb48ef29780111f9e54dd204793  33352
df<-subset_samples(ps, Compartment=="Roots" | Compartment=="Nodule" )
df<-prune_taxa(taxa_sums(df) > 0, df)
remove<-subset_taxa(df, Domain=="Unassigned" |  Phyla=="" | Phyla==" p__"  ) 
remove
# 53 taxa
unique(taxon$Domain)
########## remove these guys
badtaxa<-taxa_names(remove)
alltaxa<-taxa_names(ps)
mytaxa <- alltaxa[!(alltaxa %in% badtaxa)]
ps<-prune_taxa(mytaxa, ps )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 12652 taxa

###make data frame##
#presence in plant ~ abundance in rhizosphere
# don't include taxa that are super rare less than 10
# in 3 samples

ps1<-subset_samples(ps, Compartment != "ctl"& Fraction != "Total_DNA")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# 927 asvs
taxon<- as.data.frame(tax_table(ps1))
df<-as.data.frame(otu_table(ps1))
#remove "other" column where rare taxa went
df<-select(df, -Others)
# remove "other" row where rare taxa went
taxon<-taxon[c(which(row.names(taxon)!="Others")),]

#can't use C67 (rep4) because it doesn't have active in the rhizosphere.

#### seperate out each rep
#               abundance_active  abundance_total present_in_plant
#otu #1 rep 1
#otu #1 rep 2
n<-row.names.data.frame(df)
n[grepl("C10" , n)]="5"
n[grepl("C2" , n)]="2"
n[grepl("C5" , n)]="3"
n[grepl("C7" , n)]="4"
n[grepl("C1R" , n)]="1"
n[grepl("C1N" , n)]="1"
n[grepl("C1E" , n)]="1"

df$rep <- n

#rep 1
df1<-filter(df, rep=="1")
dim(df1)
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
head(df1)
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
head(df1)
row.names(df1)<-paste0(row.names(df1), "Rep1")
colnames(df1) <- c("Active", "Total", "inplant")
rep1<-df1
head(rep1)

##rep 2
df1<-filter(df, rep=="2")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep2")
colnames(df1) <- c("Active", "Total", "inplant")
df1
rep2<-df1

##rep 3
df1<-filter(df, rep=="3")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep3")
colnames(df1) <- c("Active", "Total", "inplant")
df1
rep3<-df1

##rep 4 only has sybr
df1<-filter(df, rep=="4")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep4")
colnames(df1) <- c("Total", "inplant")
df1
rep4<-df1
dim(rep4)



##rep 5
df1<-filter(df, rep=="5")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep5")
#change column names
colnames(df1) <- c("Active", "Total", "inplant")
rep5<-df1
dim(rep5)
rep5
### put them together
df1<-rbind(rep1, rep2)  
df1<-rbind(df1, rep3)
df1<-rbind(df1, rep5)

df1$inplant <- as.numeric(df1$inplant)
colnames(df1)
dim(df1)
head(df1)


# Models:
fit <- glm(df1$inplant ~ df1$Total, family = binomial)
summary(fit)

fit <- glm(df1$inplant ~ df1$Total+ df1$Active + df1$Active*df1$Total, family = binomial)
summary(fit)

fit <- glm(df1$inplant ~ df1$Active , family = binomial)
summary(fit)

#plots :
ggplot(df1)+
  geom_point(aes(Total, inplant))+
  stat_smooth(aes(Total, inplant), method="glm", color="#045bc2", se=FALSE, 
                method.args = list(family=binomial))+
  xlim(0,3800)+
  theme_bw()


ggplot(df1)+
  geom_point(aes(Active, inplant))+
    stat_smooth(aes(Active, inplant), method="glm", color="#045bc2", se=FALSE, 
              method.args = list(family=binomial))+
  xlim(0,3800)+
  theme_bw()

##########TOP TAXA effecting the model #######
df1<-df1[order(-df1$Active),]
df2<-filter(df1, inplant==1)


# make rep column and average
head(df2)
otu<-row.names(df2)
rep<-str_sub(otu, -4, -1)

otu<-gsub("Rep1", "", otu)
otu<-gsub("Rep2", "", otu)
otu<-gsub("Rep3", "", otu)
otu<-gsub("Rep4", "", otu)
otu<-gsub("Rep5", "", otu)

df2$otu <- otu
df2$rep <- rep

# summarise by otu
df2<-df2 %>% group_by(otu) %>% summarise(mean_active = mean(Active), sd_active=sd(Active), mean_viable=mean(Total), sd_viable=sd(Total))

df2 <- filter(df2, mean_active>0)
df2

taxon$otu <-row.names(taxon) 
top_taxa<-left_join(df2, taxon)
head(top_taxa)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/tables")
write.csv(top_taxa, "activetaxa_inplant.csv")

