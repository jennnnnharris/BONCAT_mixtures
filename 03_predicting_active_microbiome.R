#predicted active microbiome
# may 2026
# j harris

## Clear workspace ##
rm(list=ls())

## Load required libraries ##
library(tidyverse)
library(vegan)
library(phyloseq)
library(multcompView)
library(compositions)
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
# load biomass info
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/biomass")
biomass<-read.csv("biomass_percent.csv")
# only n+ 
biomass<-biomass %>% filter(N==1)
biomass 

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
## Set the working directory ###
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/predicted")
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


#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
#svg("pcoa.predicted1.svg", height = 4, width = 4)
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
legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
        fill= IBM,
        cex=.5,
        bty = "n")
legend("bottomright", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=.5,
       title = "",     bty = "o")

#dev.off()
# permanova
adonis2(otus.bray ~ Treatment, data = metadat2)
adonis2(otus.bray ~ Treatment+Measurement+Treatment*Measurement, data = metadat2, by="terms")


### pairwise adonis test #####
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

# Run the CAP (db-RDA) analysis
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
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAPactive")
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_CAP_active")
#svg("cap.active.svg", width = 5 , height = 5)
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
 legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
        fill= IBM,
        cex=.5,
        bty = "n")
legend("bottomleft", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=.5,
       title = "",     bty = "o")

#dev.off()


########## Aitchison distance PLOT  #####
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil"  )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Active" & Treatment!="Soil") %>%
  select(Treatment, Measurement, Rep)
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

metadat2
# CLR transform
clr_data <- clr(otu_table(ps1)+1)
clr_data<-as.matrix(clr_data)
clr_data[1:5, 1:5]

# caculate distance matrix
dist_mat <- vegan::vegdist(clr_data, method = "euclidean")
dist_mat<-as.matrix(dist_mat) 
dist_mat

# PCA and plot
pca1<- prcomp(clr_data)

sc_si <-scores(pca1, display="sites", choices=c(1,2), scaling=1)

#load col
IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)

# Calculate percentage of variance explained for axis labels
var_explained <- (pca1$sdev^2) / sum(pca1$sdev^2) * 100
perc1 <- paste0("PC1 (", round(var_explained[1], 2), "%)")
perc2 <- paste0("PC2 (", round(var_explained[2], 2), "%)")


# save 
# plot
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_predicted_active")
svg("active_PCA.svg",  width=6, height=6)
par(adj=.5)
ordiplot(pca1, choices=c(1,2),
         type="none",
         cex.lab = 1,
         xlab=perc1,
         ylab=perc2)

points(sc_si, 
       col= IBM[as.factor(metadat2$Treatment)],
       pch= c(16,8)[as.factor(metadat2$Measurement)],
       lwd=2,cex=1.2,
       bg=IBM[as.factor(metadat2$Treatment)],)
legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
       fill= IBM,
       cex=1,
       bty = "o")
legend("bottomright", legend=c("measured", "expected"  ),
       pch=c(16,8 ),
       cex=1,
       bty = "o")

dev.off()


###scalar projection GB #####

# clr matrix
ps1 <-subset_samples(ps, Fraction=="Active" & Treatment!="Soil"  )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Active"  & Treatment!="Soil") %>%
  select(Treatment, Measurement, Rep)
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# CLR transform
clr_data <- clr(otu_table(ps1)+1)
clr_data<-as.matrix(clr_data)
clr_data[1:5, 1:5]

# # add metadata
# key <-cbind(metadat2, clr_data)
# key[1:5, 1:5]


# 1. Calculate the centroids (mean vector) for your baseline monocultures
# (Assuming 'clr_matrix' contains only your numeric CLR-transformed columns)
centroid_G <- colMeans(clr_data[ which(metadat2$Treatment=="G"), ])
centroid_B <- colMeans(clr_data[which(metadat2$Treatment=="B"), ])

# 2. Define the axis vector from G to B
v <- centroid_B - centroid_G
v_length_sq <- sum(v^2)


# 3.  Create a function to project a sample vector onto the G->B axis  symmetric [-1, +1] axis
project_to_axis <- function(sample_vector, start_centroid, axis_vector, axis_len_sq) {
  sample_adj <- sample_vector - start_centroid
  dot_product <- sum(sample_adj * axis_vector)
  
  # Calculate 0 to 1 index
  original_index <- dot_product / axis_len_sq
  
  # Rescale to -1 to +1
  symmetric_index <- (2 * original_index) - 1
  return(symmetric_index)
}
# 4. Apply this to your target samples (GB observed and GB predicted)
# Let's say you target a data frame of your mixtures
mix_data <- clr_data[which(metadat2$Treatment=="GB"), ]
mix_data[1:5, 1:5]
dim(mix_data)



# Calculate the index for each row
GB_index <- apply(mix_data, 1, function(row) {
  project_to_axis(row, centroid_G, v, v_length_sq)
})

GB_index  
mix_data<-cbind(GB_index, mix_data)
metadat3<-metadat2 %>% filter(Treatment=="GB") 
mix_data<-cbind(metadat3, mix_data)
mix_data

# # plot raw -1 to 1
 ggplot(mix_data, aes(x = Measurement, y = GB_index, fill = Measurement)) +
   geom_boxplot(alpha = 0.6, outlier.shape = NA) +
   scale_fill_manual(values = c("#DC267F", "grey70"))+
   geom_jitter(width = 0.1, size = 2) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
   annotate("text", x = 0.5, y = -1, label = "G", hjust = 0) +
   annotate("text", x = 0.5, y = 1, label = "B", hjust = 0) +
   labs(title = "Mixture Composition Along G-to-B Axis",
        y = "Projection Index (G → B)",
        x = "Treatment") +
   theme_minimal()+
   coord_flip()
### plot with predicted as zero ###
# subtract each rep 
mix_data<-mix_data %>% select(GB_index, Treatment, Measurement, Rep) %>% filter(Rep!="1")
mix_data
P<-rep(mix_data$GB_index[which(mix_data$Measurement == "predicted")],2)
mix_data$diff_predict <- mix_data$GB_index-P 

# plot
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_index")
svg("GB_active.svg",  width=4, height=1.5)
mix_data %>% filter(Measurement=="measured") %>%
  ggplot(aes(x = Treatment, y = diff_predict, fill = Measurement)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.7, y = -1, label = "G", hjust = 0) +
  annotate("text", x = 0.7, y = 1, label = "B", hjust = 0) +
  labs(title = " ",
       y = "Projection Index (G → B)",
       x = "Mixture") +
  scale_fill_manual(values = c("#DC267F", "grey70"))+
  coord_flip()+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
dev.off()  

# stats #############

mix_data

# Run an independent t-test (Welch's t-test is default, which handles unequal variance safely)
t_test_result <- t.test(GB_index ~ Measurement, data = mix_data, alternative = "two.sided")

# Print results
print(t_test_result)



###scalar projection LB #####

# # clr matrix
# ps1 <-subset_samples(ps, Fraction=="Total" &  Rep!=2 & Treatment!="Soil"  )
# ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
# ps1
# 
# # subset metadata
# metadat2<-filter(metadat, Fraction=="Total" &  Rep!=2 & Treatment!="Soil") %>%
#   select(Treatment, Measurement, Rep)
# metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
# 
# # CLR transform
# clr_data <- clr(otu_table(ps1)+1)
# clr_data<-as.matrix(clr_data)
# clr_data[1:5, 1:5]

# 1. Calculate the centroids (mean vector) for your baseline monocultures
# (Assuming 'clr_matrix' contains only your numeric CLR-transformed columns)
centroid_L <- colMeans(clr_data[ which(metadat2$Treatment=="L"), ])
centroid_B <- colMeans(clr_data[which(metadat2$Treatment=="B"), ])

# Define the axis vector from L to B
v <- centroid_B - centroid_L
v_length_sq <- sum(v^2)


#  Apply this to your target samples (LB observed and LB predicted)
# Let's say you target a data frame of your mixtures
mix_data <- clr_data[which(metadat2$Treatment=="LB"), ]
mix_data[1:5, 1:5]
dim(mix_data)
mix_data[1:5,1:5]


# Calculate the index for each row
LB_index <- apply(mix_data, 1, function(row) {
  project_to_axis(row, centroid_L, v, v_length_sq)
})

LB_index  
mix_data<-cbind(LB_index, mix_data)
metadat3<-metadat2 %>% filter(Treatment=="LB") 
mix_data<-cbind(metadat3, mix_data)

# raw plot
# # plot raw -1 to 1
ggplot(mix_data, aes(x = Measurement, y = LB_index, fill = Measurement)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("#FE6100", "grey70"))+
  geom_jitter(width = 0.1, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.5, y = -1, label = "L", hjust = 0) +
  annotate("text", x = 0.5, y = 1, label = "B", hjust = 0) +
  labs(title = "Mixture Composition Along G-to-B Axis",
       y = "Projection Index (L → B)",
       x = "Treatment") +
  theme_minimal()+
  coord_flip()


#### plot with predicted as zero ###
# subtract each rep #
mix_data<-mix_data %>% select(LB_index, Treatment, Measurement, Rep)
mix_data
P<-rep(mix_data$LB_index[which(mix_data$Measurement == "predicted")],2)
mix_data$diff_predict <- mix_data$LB_index-P 
mix_data

# plot
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_index")
#svg("LB_active.svg",  width=4, height=1.5)
mix_data %>% filter(Measurement=="measured") %>%
  ggplot(aes(x = Treatment, y = diff_predict, fill = Measurement)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.7, y = -1, label = "L", hjust = 0) +
  annotate("text", x = 0.7, y = 1, label = "B", hjust = 0) +
  labs(title = "",
       y = "Projection Index (L → B)",
       x = "Mixture") +
  scale_fill_manual(values = c("#FE6100", "grey70"))+
  coord_flip()+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
#dev.off()  


mix_data

# Run an independent t-test (Welch's t-test is default, which handles unequal variance safely)
t_test_result <- t.test(LB_index ~ Measurement, data = mix_data, alternative = "two.sided")

# Print results
print(t_test_result)



###scalar projection LG #####

# # clr matrix
# ps1 <-subset_samples(ps, Fraction=="Total" &  Rep!=2 & Treatment!="Soil"  )
# ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
# ps1
# 
# # subset metadata
# metadat2<-filter(metadat, Fraction=="Total" &  Rep!=2 & Treatment!="Soil") %>%
#   select(Treatment, Measurement, Rep)
# metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
# 
# # CLR transform
# clr_data <- clr(otu_table(ps1)+1)
# clr_data<-as.matrix(clr_data)
# clr_data[1:5, 1:5]

# 1. Calculate the centroids (mean vector) for your baseline monocultures
# (Assuming 'clr_matrix' contains only your numeric CLR-transformed columns)
centroid_L <- colMeans(clr_data[ which(metadat2$Treatment=="L"), ])
centroid_G <- colMeans(clr_data[which(metadat2$Treatment=="G"), ])

# 2. Define the axis vector from L to G
v <- centroid_G - centroid_L
v_length_sq <- sum(v^2)


# 4. Apply this to your target samples (LB observed and LB predicted)
# Let's say you target a data frame of your mixtures
mix_data <- clr_data[which(metadat2$Treatment=="LG"), ]
mix_data[1:5, 1:5]
dim(mix_data)
mix_data[1:5,1:5]


# Calculate the index for each row
index <- apply(mix_data, 1, function(row) {
  project_to_axis(row, centroid_L, v, v_length_sq)
})

index  
mix_data<-cbind(index, mix_data)
meta3<-metadat2 %>% filter(Treatment=="LG") 
mix_data<-cbind(meta3, mix_data)


#### plot with predicted as zero ###
# subtract each rep #
mix_data<-mix_data %>% select(index, Treatment, Measurement, Rep)
mix_data

P<-rep(mix_data$index[which(mix_data$Measurement == "predicted")],2)
mix_data$diff_predict <- mix_data$index-P 
mix_data



#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_index")
#svg("LG_active.svg", width=4, height=1.5)
mix_data %>% filter(Measurement=="measured") %>%
  ggplot(aes(x = Treatment, y = diff_predict, fill = Measurement)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.7, y = -1, label = "L", hjust = 0) +
  annotate("text", x = 0.7, y = 1, label = "G", hjust = 0) +
  labs(title = "",
       y = "Projection Index (L → G)",
       x = "Mixture") +
  scale_fill_manual(values = c("#FFB000", "grey70"))+
  coord_flip()+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

#dev.off()  


# Run an independent t-test (Welch's t-test is default, which handles unequal variance safely)
t_test_result <- t.test(index ~ Measurement, data = mix_data, alternative = "two.sided")

# Print results
print(t_test_result)


######simple ordination three species ####

# LGB, G, B, G example
ps1 <-subset_samples(ps) %>% subset_samples(Treatment=="L" | Treatment=="G" | Treatment=="B" | Treatment=="LGB")    
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-metadat %>% filter(Treatment=="L" | Treatment=="G" | Treatment=="B" | Treatment=="LGB") %>%
  select(Treatment, Measurement, Rep)
metadat2
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "L", "G", "B", "LGB"))
#metadat2$Trt_measure <- paste0(metadat2$Treatment, metadat2$Measurement)



# CLR transform
clr_data <- clr(otu_table(ps1)+1)
clr_data<-as.matrix(clr_data)
clr_data[1:5, 1:5]

# caculate distance matrix
dist_mat <- vegan::vegdist(clr_data, method = "euclidean")
dist_mat<-as.matrix(dist_mat) 
dist_mat

# add metadata
key <-cbind(metadat2, dist_mat)
key

# PCA and plot
pca1<- prcomp(clr_data)

sc_si <-scores(pca1, display="sites", choices=c(1,2), scaling=1)

#load col
IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  #"#DC267F", # magenta pink 
  #"#FE6100", # bright orange
  # "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)

# Calculate percentage of variance explained for axis labels
var_explained <- (pca1$sdev^2) / sum(pca1$sdev^2) * 100
perc1 <- paste0("PC1 (", round(var_explained[1], 2), "%)")
perc2 <- paste0("PC2 (", round(var_explained[2], 2), "%)")

par(adj=.5)
ordiplot(pca1, choices=c(1,2),
         type="none",
         cex.lab = 1,
         xlab=perc1,
         ylab=perc2)

points(sc_si, 
       col= IBM[as.factor(metadat2$Treatment)],
       pch= c(16,8)[as.factor(metadat2$Measurement)],
       lwd=2,cex=1.2,
       bg=IBM[as.factor(metadat2$Treatment)],)
legend("topright", legend=c("L", "G", "B", "LGB"),
       fill= IBM,
       cex=1,
       bty = "n")
legend("bottomright", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=1,
       bty = "o")


##### triplot scalar projection LGB ########

# data
clr_data

# 1. Calculate the centroids (mean vector) for your baseline monocultures
# (Assuming 'clr_matrix' contains only your numeric CLR-transformed columns)
centroid_L <- colMeans(clr_data[ which(metadat2$Treatment=="L"), ])
centroid_G <- colMeans(clr_data[which(metadat2$Treatment=="G"), ])
centroid_B <- colMeans(clr_data[which(metadat2$Treatment=="B"), ])


# 2. Define the axis vector from L to G
# v1 <- end - start
# L to G
v1 <- centroid_G - centroid_L
# G to B
v2 <- centroid_B - centroid_G
# B to L
v3 <- centroid_L - centroid_B

v1_length_sq <- sum(v1^2)
v2_length_sq <- sum(v2^2)
v3_length_sq <- sum(v3^2)

# Apply this to your target samples (LGB observed and LGB predicted)
# Let's say you target a data frame of your mixtures
mix_data <- clr_data[which(metadat2$Treatment=="LGB"), ]
mix_data[1:5, 1:5]
dim(mix_data)


# Create a function to project a sample vector onto the G->B axis  symmetric [-1, +1] axis
project_to_axis <- function(sample_vector, start_centroid, axis_vector, axis_len_sq) {
  sample_adj <- sample_vector - start_centroid
  dot_product <- sum(sample_adj * axis_vector)
  
  # Calculate 0 to 1 index
  original_index <- dot_product / axis_len_sq
  
  # Rescale to -1 to +1
  #symmetric_index <- (2 * original_index) - 1
  return(original_index)
}

# row, start, v, length
LG_index <- apply(mix_data, 1, function(row) {
  project_to_axis(row, centroid_L, v1, v1_length_sq)
})

GB_index <- apply(mix_data, 1, function(row) {
  project_to_axis(row, centroid_G, v2, v2_length_sq)
})

BL_index <- apply(mix_data, 1, function(row) {
  project_to_axis(row, centroid_B, v3, v3_length_sq)
})

# add metadata
mix_data<-cbind(LG_index, mix_data)
mix_data<-cbind(GB_index, mix_data)
mix_data<-cbind(BL_index, mix_data)
meta3<-metadat2 %>% filter(Treatment=="LGB") 
mix_data<-cbind(meta3, mix_data)
mix_data<-mix_data %>% select(LG_index, GB_index, BL_index, Treatment, Rep, Measurement)
mix_data



library(ggtern)
#library(ggplot2)


# View one of the triads
# Use ggtern which is a specialized package using ggplot2 tools.
# You will define the data source (demodata1) and then provide the three column names where ggtern
# will find the x, y, and z coordinate data (T1A, T1B, and T1C)
# Then you tell ggtern what kind of geometry to use to display the data, in this case "geom_point"


p4<-ggtern(data=mix_data, aes(x=LG_index, y=GB_index, z=BL_index, colour = Measurement, shape = factor(Rep))) +
  geom_point(size=2)+
  theme_minimal()+
  scale_color_manual(values = c( "#865338", "grey70"), labels=c("measured", "expectation"))+
  scale_shape_manual(values=c(8, 15, 17, 19, 9, 18), name="Rep")+
  xlab("L")  +                  
  ylab("G") +
  zlab("B")   


p4


setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_index")
svg("LGB_active.svg", width=4, height=4) 
p4  
dev.off()




