## figure 6 combining data sets

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
  "#865338" # medium mocha brown
)


####### import activity data/ weed seed / bio masss #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto")
fc <-read.csv("processed_flow_cyto.csv")
head(fc)

#make dates be dates
fc$Date_Sorted<-mdy(fc$Date_Sorted)

# remove outliers# filter out day were pos ctl didn't work
fc <- filter(fc, Trt_ID!="B+N6")
fc<-filter(fc, Date_Sorted != "2023-05-25")
head(fc)

# make treatment and day factors
fc$Treatment   <- factor(fc$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
fc$Date_Sorted   <- factor(fc$Date_Sorted)

fc<-fc %>% select(Rep,	Group,	Treatment,	Pot_ID,	Nitrogen,	Block,	Trt_ID,	Legume_long,	Legume,	Brassicae,
              Grass, boncat_freq, active_cel_per_g)


# import weed seed
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed<-read.csv("weed_seed_decay.csv", row.names = 1)
weed
fc
df<- left_join(fc, weed)
head(df)

# import biomass
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("biomass_potlevel.csv")
colnames(biomass)
df<-left_join(df, biomass)
df


####Import 16S data #####
# Set the working directory 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")
taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.table("metadat.txt",sep="\t", header=T)

# Transpose ASVS table 
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

# order metadata
head(metadat)
row.names(metadat) <- metadat$SampleID
metadat

# filter for just flow cyto samples 
asvs<-asvs[which(metadat$Fraction!="Total"),]
metadat<-metadat[which(metadat$Fraction!="Total"),]
dim(metadat)
metadat$Legume   <- factor(metadat$Legume)
metadat$Brassicae   <- factor(metadat$Brassicae)
metadat$Grass   <- factor(metadat$Grass)

metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
metadat$Legume_label <- factor(metadat$Legume_label, levels= c("Legumes present", "Legumes absent"))

# make taxon matrix row names OTUs
row.names(taxon) <- taxon$asv

# get min number of reads in a sample
min.s<-min(rowSums(asvs))

# Rarefy to obtain even numbers of reads by sample
set.seed(336)
asvs.r<-rrarefy(asvs, min.s)

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs.r), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps# 110K taxa when rarefied 

# remove plant contamination
ps<-subset_taxa(ps, Class!="c__Chloroplast" )
ps<-subset_taxa(ps, Class!=" c__Chloroplast" )
ps<-subset_taxa(ps, Family!= " f__Mitochondria" )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps

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
metadat2$Pot_ID<-as.numeric(metadat2$Pot_ID)
metadat2

df<-left_join( metadat2, df)

df<-df %>% select(-Nitrogen, -Nitrogen_label, -X, -n_species, -Rep, -Block, -Group, -SampleID, -Pot_ID)



# z transform
df$z_Shoot.Biomass<-as.vector(scale(df$Stem.Biomass.g))
df$z_active_cel_per_g<-as.vector(scale(log(df$active_cel_per_g)))
df$z_Root.Biomass<-as.vector(scale(df$Total.Root.g))

hist(df$z_Root.Biomass)
hist(df$z_Shoot.Biomass)
hist(df$z_active_cel_per_g)

# 4. Add it to your data frame
df$PC1 <- pc1_variable
df$PC2 <- pc2_variable



# Weed Seed decay ################

### model selection
numeric_data<-df %>% select(PC1, PC2, z_active_cel_per_g, z_Root.Biomass, z_Shoot.Biomass, Legume)
# Select only numeric columns
numeric_data <- numeric_data[sapply(numeric_data, is.numeric)]
# Calculate the correlation matrix
cor_matrix <- cor(numeric_data, use = "complete.obs")
library(corrplot)
# Basic circle plot
corrplot(cor_matrix, method = "circle")

# A more professional "mixed" plot
corrplot.mixed(cor_matrix, 
               lower = "number", 
               upper = "circle", 
               tl.col = "black")



library(car)

#full model
y<-cbind(df$pigweed_num_nongerm, df$pigweed_num_germ)
m1<-glm(data= df, y~(PC1+PC2+z_active_cel_per_g), family = binomial)
summary(m1)#
vif_values <- vif(m1)
print(vif_values)
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue", las = 1)
abline(v = 5, lwd = 3, lty = 2, col = "red") # Add a threshold line at 5
# check vif 


#full model
y<-cbind(df$pigweed_num_nongerm, df$pigweed_num_germ)
m1<-glm(data= df, y~(PC1+PC2+z_active_cel_per_g)*Legume, family = binomial)
summary(m1)

#reduced model with interactions
y<-cbind(df$pigweed_num_nongerm, df$pigweed_num_germ)
m1<-glm(data= df, y~(PC2+z_active_cel_per_g)*Legume, family = binomial)
summary(m1)# PC1 and active cell per g are correlated so maybe don't have them both in the model


# legume present
df1<-df %>% filter(Legume=="1")
y<-cbind(df1$pigweed_num_nongerm, df1$pigweed_num_germ)
m1<-glm(data= df1, y~(PC1*PC2*z_active_cel_per_g), family = binomial)
summary(m1) 
# reduced model 
df1<-df %>% filter(Legume=="1")
y<-cbind(df1$pigweed_num_nongerm, df1$pigweed_num_germ)
m1<-glm(data= df1, y~(z_active_cel_per_g), family = binomial)
summary(m1)

# legume absent
df1<-df %>% filter(Legume=="0")
y<-cbind(df1$pigweed_num_nongerm, df1$pigweed_num_germ)
m1<-glm(data= df1, y~PC1*PC2*z_active_cel_per_g, family = binomial)
summary(m1) 
df1<-df %>% filter(Legume=="0")
y<-cbind(df1$pigweed_num_nongerm, df1$pigweed_num_germ)
m1<-glm(data= df1, y~PC2*z_active_cel_per_g, family = binomial)
summary(m1)




# reduced model 1
y<-cbind(df$pigweed_num_nongerm, df$pigweed_num_germ)
m1<-glm(data= df, y~PC1*PC2*z_active_cel_per_g, family = binomial)
summary(m1)

# reduced model 2
y<-cbind(df$pigweed_num_nongerm, df$pigweed_num_germ)
m1<-glm(data= df, y~PC2*z_active_cel_per_g, family = binomial)
summary(m1)


df$pigweed_perc<-df$pigweed_prop_nongerm*100

p1<-df %>% filter(Treatment!="Soil" & Legume!="NA") %>%
  ggplot( aes(y=pigweed_perc, x=active_cel_per_g, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  theme_bw(base_size = 12) +  facet_wrap(~Legume_long)+
  scale_x_continuous(limits=c(0,1800))+
  labs(y= "non-germinating pigweed (%)",
       x= "Active cells / g rhizosphere")
  
p1

p2<-df %>% filter(Treatment!="Soil" & Legume!="NA") %>%
  ggplot( aes(y=pigweed_perc, x=PC2, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  theme_bw(base_size = 12)   + facet_wrap(~Legume_long) +
  labs(y= "non-germinating pigweed (%)",
       x= "Active PC2")
p2



setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig6_mbiome_functional")
svg("weed_seeed.svg", width=6, height=6)
require(gridExtra)
grid.arrange(p1, p2,ncol=1)
dev.off()






# Shoot biomass  ################

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
m1<-lm(data= df, Stem.Biomass.g~PC1*PC2*z_active_cel_per_g)
summary(m1) 

m1<-lm(data= df, Stem.Biomass.g~PC1*z_active_cel_per_g)
summary(m1) 

m1<-lm(data= df, Stem.Biomass.g~PC1)
summary(m1) 


# ggplot
df<-df %>% filter(Treatment!="Soil")
p1<-ggplot(df, aes(y=Stem.Biomass.g, x=active_cel_per_g, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  #annotate("text", x = 1500, y = 5, label = "p<0.001, Rsq=.38", 
  #         color = "black", size = 5, fontface = "bold")+
  labs(y = "Shoot biomass (g)",
       x = "Active cells / g rhizosphere")+
  theme_bw(base_size = 12)#+
  #geom_smooth(method = "lm", se = FALSE, color = "black")
p1

p2<-df %>%
  ggplot(aes(y=Stem.Biomass.g, x=PC1, col=Treatment))+
  geom_point()+
  scale_color_manual(values=IBM)+
  theme_bw()+
  labs(y = "Shoot biomass (g)",
       x= "Active PC1")
 # geom_smooth(method = "lm", color = "black")
p2


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig6_mbiome_functional")
 svg("shoot.active.svg", width=6, height = 6)
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




