## figure 6 combining data sets total dna


### Initial Setup ###
rm(list=ls())

## Load required libraries ##

library(tidyverse)
library(vegan)
library(phyloseq)
library(multcompView)
library(BiodiversityR)

## set colors

IBM <- c( #IBM colors
  
  "navy", # dark royal blue
  "#a4bdfc",  # "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" ,# medium mocha brown,
  "grey"
)


####### import activity data/ weed seed / bio masss #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto")
fc <-read.csv("processed_flow_cyto.csv")
head(fc)

#make dates be dates
fc$Date_Sorted<-mdy(fc$Date_Sorted)

# remove outliers# filter out day were pos ctl didn't work
# fc <- filter(fc, Trt_ID!="B+N6")
# fc<-filter(fc, Date_Sorted != "2023-05-25")
# head(fc)

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)

fc<-fc %>% select(	Trt_ID, boncat_freq, active_cel_per_g)


# import weed seed
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed<-read.csv("weed_seed_decay.csv")
weed

df<- left_join(weed, fc)
head(df)

# import biomass
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("biomass_potlevel.csv")
colnames(biomass)
biomass<-biomass %>% select(Trt_ID, Stem.Biomass.g, Root.Biomass.g, Total.Root.g)

df<-left_join(df, biomass)
df


####Import 16S data #####
# Set the working directory 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")

taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat2.csv", header = T, row.names = 1)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)


#metadat<-as.data.frame(metadat[order(metadat$SampleID),])
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
row.names(metadat) <- metadat$SampleID
head(metadat)

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
# 

#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps #1902 taxa 



#### extract PC1 
# filter data 
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps
# filter for active fraction 
ps1 <-subset_samples(ps, N!=0)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1 #1785 taxa
# filter metadata
metadat2<-filter(metadat,   N!=0)
metadat2
# Calculate the distance matrix (e.g., Bray-Curtis)
dist_matrix<-vegdist(otu_table(ps1), method = "bray")

# Perform PCA
# scale. = TRUE is critical to ensure all variables are treated equally
pca_result <- prcomp(dist_matrix, center = TRUE, scale. = TRUE)


#######add pcs
# Extract the first principal component (PC1)
# This is your new one-dimensional variable
pc1_variable <- pca_result$x[, 1]
pc2_variable <- pca_result$x[,2]

# get variables of interest - weed seed decay, biomass, nfix, activity
# clean up metadat2
metadat2<-metadat2 %>% select(Trt_ID,)
#metadat2$Pot_ID<-as.numeric(metadat2$Pot_ID)
metadat2$PC1 <- pc1_variable
metadat2$PC2 <- pc2_variable
metadat2
df<-df %>%
  mutate(df, foxtail_prop_germ = foxtail_num_germ/foxtail_total_seeds,
         pigweed_prop_germ= pigweed_num_germ/pigweed_total_seeds)




df<-left_join(df, metadat2)


# z transform
df$z_Shoot.Biomass<-as.vector(scale(df$Stem.Biomass.g))
df$z_active_cel_per_g<-as.vector(scale(log(df$active_cel_per_g)))
df$z_Root.Biomass<-as.vector(scale(df$Total.Root.g))
df$z_boncat_freq <- as.vector(scale(df$boncat_freq))
df$z_pigweed <- as.vector(scale(df$pigweed_prop_germ))
df$z_foxtail <- as.vector(scale(log(df$foxtail_prop_germ)))

hist(df$z_active_cel_per_g)
hist(df$z_boncat_freq)
hist(df$z_pigweed)
hist(df$PC1)
hist(df$PC2)
# Weed Seed decay #

# check
m1<-lm(data= df, z_pigweed~z_active_cel_per_g)
summary(m1)
plot(m1)

m1<-lm(data= df, z_pigweed~z_boncat_freq)
summary(m1)
plot(m1)

m1<-lm(data= df, z_pigweed~PC1)
summary(m1)

m1<-lm(data= df, z_pigweed~PC2)
summary(m1)
plot(m1)

quad<-lm(data= df, z_pigweed~PC1+ I(PC1^2))
summary(quad)

vif_values <- vif(m1)
print(vif_values)
plot(m1)

anova(m1, quad)

df$pigweed_perc<-df$pigweed_prop_germ*100

p1<-df %>% 
  ggplot( aes(y=z_pigweed, x=z_boncat_freq, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  #geom_smooth(method = lm, col="grey")+
  theme_bw(base_size = 12)+
  labs(y= "z transformed pigweed seed decay",
       x= "z transformed % active microbes")
p1

p2<-df %>% 
  ggplot( aes(y=pigweed_prop_germ, x=z_active_cel_per_g, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  #geom_smooth(method = lm, col="grey")+
  theme_bw(base_size = 12)+
  labs(y= "z transformed pigweed seed decay",
       x= "z transformed active cells/ g")
p2

p2<-df %>% 
  ggplot( aes(y=pigweed_prop_germ, x=PC1, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  #geom_smooth(method = lm, col="grey")+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col="grey") +
  theme_bw(base_size = 12)+
  labs(y= "z transformed pigweed seed decay",
       x= "Total rhizosphere community PC1 ")
p2

p2<-df %>% 
  ggplot( aes(y=pigweed_prop_germ, x=PC1, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  #geom_smooth(method = lm, col="grey")+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col="grey") +
  theme_bw(base_size = 12)+
  labs(y= "z transformed pigweed seed decay",
       x= "Total rhizosphere community PC1 ")
p2

p2<-df %>% 
  ggplot( aes(y=pigweed_prop_germ, x=PC2, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  #geom_smooth(method = lm, col="grey")+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col="grey") +
  theme_bw(base_size = 12)
p2



setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig6_mbiome_correlations")
svg("weed_seed_total.svg", width=5, height=6)
require(gridExtra)
grid.arrange(p1, p2,ncol=1)
dev.off()






# root biomass  ################

### model selection
numeric_data<-df %>% select(PC1, PC2, boncat_freq, active_cel_per_g, Stem.Biomass.g, Root.Biomass.g)
# Select only numeric columns
numeric_data <- numeric_data[sapply(numeric_data, is.numeric)]

# Calculate the correlation matrix
# 'use = "complete.obs"' handles missing values by ignoring those rows
cor_matrix <- cor(numeric_data, use = "complete.obs")
# Basic circle plot
corrplot(cor_matrix, method = "circle")

# A more professional "mixed" plot
corrplot.mixed(cor_matrix, 
               lower = "number", 
               upper = "circle", 
               tl.col = "black")


# biomass verse activity 
#
m1<-lm(data= df, z_Root.Biomass~PC1)
summary(m1) 

m1<-lm(data= df, z_Root.Biomass~z_boncat_freq)
summary(m1) 

m1<-lm(data= df, z_Root.Biomass~z_active_cel_per_g)
summary(m1) 

m1<-lm(data= df, z_Root.Biomass~z_boncat_freq*PC1)
summary(m1) 


# ggplot
df<-df %>% filter(Treatment!="Soil")
p1<-ggplot(df, aes(y=z_Root.Biomass, x=z_boncat_freq, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  labs(y = "Z transformed root biomass",
       x = "z transformed %active microbes")+
  theme_bw(base_size = 12)+
  geom_smooth(method = "lm", se = FALSE, color = "grey")
p1


p2<-ggplot(df, aes(y=z_Root.Biomass, x=z_active_cel_per_g, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  labs(y = "Z transformed root biomass",
       x = "z transformed active cells/g")+
  theme_bw(base_size = 12)
  #geom_smooth(method = "lm", se = FALSE, color = "grey")
p2



df<-df %>% filter(Treatment!="Soil")
p2<-ggplot(df, aes(y=z_Root.Biomass, x=PC1, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  labs(y = "Z transformed root biomass",
       x = "Total rhizosphere community PC1")+
  theme_bw(base_size = 12)+
geom_smooth(method = "lm", se = FALSE, color = "grey")
p2





setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig6_mbiome_correlations")

svg("root.active.total.svg", width=5, height = 6)
require(gridExtra)
grid.arrange(p1, p2 , ncol=1)
dev.off()

####stats #######
# full model
df$z_active_cel_per_g
m1<-lm(data= df, z_Shoot.Biomass~(PC1+PC2+z_active_cel_per_g)*Legume )
summary(m1)  
# main effect is legume 


### legume present
df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~(PC1*PC2*z_active_cel_per_g) )
summary(m1)

df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~(PC1) )
summary(m1)

df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~(PC2) )
summary(m1)

df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~z_active_cel_per_g)
summary(m1)


### legume absent
df1<-df %>% filter(Legume=="0")
m1<-lm(data= df1, z_Shoot.Biomass~(PC1*PC2*z_active_cel_per_g) )
summary(m1)

df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~(PC1) )
summary(m1)

df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~(PC2) )
summary(m1)

df1<-df %>% filter(Legume=="1")
m1<-lm(data= df1, z_Shoot.Biomass~z_active_cel_per_g)
summary(m1)


# don't focus on legumes

df$z_active_cel_per_g
m1<-lm(data= df, z_Shoot.Biomass~(PC1+PC2+z_active_cel_per_g))
summary(m1)  
# main effe


### root?
# full model
df$z_active_cel_per_g
m1<-lm(data= df, z_Root.Biomass~(PC1+PC2+z_active_cel_per_g)*Legume )
summary(m1)  
# main effect is legume 




