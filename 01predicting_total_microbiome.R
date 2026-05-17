# predicting total composition in mixtures from monocultures
# J. Harris
# May 2026


########Initial Setup ##################

### Clear workspace ###
rm(list=ls())

# Load required libraries #
library(tidyverse)
library(vegan)
library(readxl)
library(phyloseq)
#library(multcompView)
#library(BiodiversityR)


#import data# ###############
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing")
taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat2.csv", header = T, row.names = 1)

# Transpose ASVS table #
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

# order metadata #
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

#T_DNA_23_S153 is omitted. It had few reads (13,000) and low DNA concentration.
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
ps<-subset_taxa(ps, Class!="c__Chloroplast" )
ps<-subset_taxa(ps, Class!=" c__Chloroplast" )
ps<-subset_taxa(ps, Family!= " f__Mitochondria" )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 220 taxa and 86 samples when mitochondria removed. 

#total
ps<-subset_samples(ps, Fraction=="Total")
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 220 k asvs




##############predicting microbial community #############


# filtering #
# remove singletons and remove asvs with a mean of less than 5, and filter for N+ samples only 
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps<-subset_samples(ps, N==1)
ps

# import biomass data  #
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
biomass<-read.csv("percent.biomass.csv")
# filter for only n+ and required columns
biomass<-biomass %>% filter(N==1) %>% select(-Total.Root.g, -Total.withbulk, -Stem.Biomass.g, -Root.Biomass.g, -Bulk.Root.g)

# Note: interpolate the mean is for missing biomass observation LB mixture rep 5 #




# predict for LG #################

# 1. filter for L community and save as a data frame
ps1<-subset_samples(ps , Treatment=="L")
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # view sample data, appears to be no rep 2 

# 2. arrange by rep
sort_by_rep <- function(data, part_index = 3, split_char = "_") {
  # Get the names 
  row_names_vec <- row.names(data)
  
  # Split the strings and extract the specific part
  split_matrix <- str_split(row_names_vec, split_char, simplify = TRUE)
  
  # Add the temporary Pot_ID column, ensuring it is numeric for proper sorting
  data$temp_sort_id <- as.numeric(split_matrix[, part_index])
  
  # Sort the dataframe and then remove the temporary column
  data_sorted <- data %>%
    arrange(temp_sort_id) %>%
    select(-temp_sort_id)
  
  return(data_sorted)
}

# use function
df<-sort_by_rep(df)

# 3. get biomass info 
sp1<-biomass %>%
  filter(Treatment=="LG") %>%
  filter(Species=="legume") %>%
  arrange(Rep) %>% 
  filter(Rep!=2) # remove rep 2 because we are missing the microbial composition

# 4. multiply by Legume biomass in LG
multiply_by_percent <- function(data, percent_df, percent_col = "percent") {
  library(dplyr)
  
  # 1. Extract and transform the vector (convert 0-100 to 0-1)
  # This assumes percent_df has a column named 'percent'
  vec <- percent_df[[percent_col]] / 100
  
  # 2. Safety Check: Does the vector length match the data row count?
  if (length(vec) != nrow(data)) {
    stop(paste("Row mismatch! Data has", nrow(data), 
               "rows, but percent vector has", length(vec)))
  }
  
  # 3. Apply the multiplication across all columns
  data_multiplied <- data %>%
    mutate(
      across(
        everything(),
        .fns = ~ .x * vec
      )
    )
  
  return(data_multiplied)
}

# use function
data_frame_multiplied1 <- multiply_by_percent(df, sp1)


# 5. filter for G community
ps1<- subset_samples(ps , Treatment=="G" & Rep!=2)
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1)) # all reps present, but filter out rep 2 because that is missing from L

# 6. arrange by rep 
df<-sort_by_rep(df)

# 7. get biomass info 
head(biomass)
sp1<-biomass %>%
  filter(Treatment=="LG") %>%
  filter(Species=="grass") %>%
  arrange(Rep) %>% 
  filter(Rep!=2) # remove rep 6 

# 8 multiply by biomass
data_frame_multiplied2<- multiply_by_percent(df, sp1)

# 9. add together
LG<-data_frame_multiplied1 + data_frame_multiplied2
row.names(LG)<-c("predict_LG_N1", "predict_LG_N3", "predict_LG_N4", "predict_LG_N5", "predict_LG_N6")
LG[1:5, 1:5]

# predict for GB ###################

# 1. filter for G community and save as a data frame
ps1<-subset_samples(ps , Treatment=="G")
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # view sample data, all reps present

# 2. arrange by rep
df<-sort_by_rep(df)

# 3. get biomass info 
sp1<-biomass %>%
  filter(Treatment=="GB") %>%
  filter(Species=="grass") %>%
  arrange(Rep) 

# 4. multiply by biomass 
data_frame_multiplied1 <- multiply_by_percent(df, sp1)


# 5. filter for B community
ps1<- subset_samples(ps , Treatment=="B" )
df<-as.data.frame(otu_table(ps1))
#as.data.frame(sample_data(ps1)) # all reps present

# 6. arrange by rep 
df<-sort_by_rep(df)

# 7. get biomass info 
head(biomass)
sp1<-biomass %>%
  filter(Treatment=="GB") %>%
  filter(Species=="brassica") %>%
  arrange(Rep) 

# 8 multiply by  biomass
data_frame_multiplied2<- multiply_by_percent(df, sp1)
# add together
GB<-data_frame_multiplied1 + data_frame_multiplied2
row.names(GB)<-c("predict_GB_N1","predict_GB_N2", "predict_GB_N3", "predict_GB_N4", "predict_GB_N5", "predict_GB_N6")


# predict for LB ##########
# 1. filter for L community and save as a data frame
ps1<-subset_samples(ps , Treatment=="L")
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # view sample data, one rep is missing (rep 2)

# 2. arrange by rep
df<-sort_by_rep(df)

# 3. get biomass info 
sp1<-biomass %>%
  filter(Treatment=="LB") %>%
  filter(Species=="legume") %>%
  arrange(Rep) %>%
  filter(Rep!=2)# remove rep 2 because we are missing the microbial composition

# 4. multiply by L  biomass in LB mix 
data_frame_multiplied1 <- multiply_by_percent(df, sp1)


# 5. filter for B community
ps1<- subset_samples(ps , Treatment=="B" & Rep!=2 )
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1)) # all reps present, but filter out rep 2 because it's missing in L

# 6. arrange by rep 
df<-sort_by_rep(df)

# 7. get biomass info 
head(biomass)
sp1<-biomass %>%
  filter(Treatment=="GB") %>%
  filter(Species=="brassica") %>%
  arrange(Rep) %>%
  filter(Rep!=2)# remove rep 2 because we are missing the microbial composition in L

# 8. multiply by  biomass
data_frame_multiplied2<- multiply_by_percent(df, sp1)
# add together
LB<-data_frame_multiplied1 + data_frame_multiplied2
row.names(LB)<-c("predict_LB_N1", "predict_LB_N3", "predict_LB_N4", "predict_LB_N5", "predict_LB_N6")



# predict for LBG ########## 
# 1. filter for L community and save as a data frame
ps1<-subset_samples(ps , Treatment=="L")
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1))   # view sample data, one rep is missing (rep 2)

# 2. arrange by rep
df<-sort_by_rep(df)

# 3. get biomass info 
sp1<-biomass %>%
  filter(Treatment=="LGB") %>%
  filter(Species=="legume") %>%
  arrange(Rep) %>%
  filter(Rep!=2)# remove rep 2 because we are missing the microbial composition

# 4. multiply by biomass
data_frame_multiplied1 <- multiply_by_percent(df, sp1)


# 5. filter for G community
ps1<- subset_samples(ps , Treatment=="G" & Rep!=2 )
df<-as.data.frame(otu_table(ps1))
as.data.frame(sample_data(ps1)) # all reps present, but filter out rep 2 because it's missing in L

# 6. arrange by rep 
df<-sort_by_rep(df)

# 7. get biomass info 
sp1<-biomass %>%
  filter(Treatment=="LGB") %>%
  filter(Species=="grass") %>%
  arrange(Rep) %>%
  filter(Rep!=2)# remove rep 2 because we are missing the microbial composition in L

# 8 multiply by biomass
data_frame_multiplied2<- multiply_by_percent(df, sp1)

# 9. filter for B community
ps1<- subset_samples(ps , Treatment=="B" & Rep!=2 )
df<-as.data.frame(otu_table(ps1))

# 10. arrange by rep 
df<-sort_by_rep(df)

# 11. get biomass info 
sp1<-biomass %>%
  filter(Treatment=="LGB") %>%
  filter(Species=="brassica") %>%
  arrange(Rep) %>%
  filter(Rep!=2)# remove rep 2 because we are missing the microbial composition in L

# 12. multiply by  biomass
data_frame_multiplied3<- multiply_by_percent(df, sp1)

# make data frame
LGB<-data_frame_multiplied1 + data_frame_multiplied2 + data_frame_multiplied3
row.names(LGB)<-c("predict_LGB_N1", "predict_LGB_N3", "predict_LGB_N4", "predict_LGB_N5", "predict_LGB_N6")


df1<-rbind(LG, GB, LB, LGB)
row.names(df1)
df1[1:5, 1:5]



#### write files
df<-as.data.frame(otu_table(ps))
otus<-rbind(df, df1)
otus<-t(otus)
row.names(otus)
colnames(otus)
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/predicted")
write.csv(otus, "total.feature.table.csv")

#taxon
tax_table<-as.data.frame(tax_table(ps))
write.csv(tax_table, "total.taxonomy.csv")



##### import predicted #####
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/predicted")
taxon <- read.csv("total.taxonomy.csv", row.names = 1)
asvs <- read.csv("total.feature.table.csv", row.names = 1)
metadat<-read.csv("total.metadata_predicted16S.csv", header = T)

## Transpose ASVS table ##
asvs[1:5,1:5] #taxa are rows
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
head(metadat)

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 1604 taxa when rarefied 


##### predicted plot #####
ps1 <-subset_samples(ps, Fraction=="Total" & Treatment!="Soil" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil" )

#factor
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))


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

# load colors

IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)

# PCoA 
par(adj=.5)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="",
         xlab=paste("PCoA1 (",pe1,"% var. explained)"),
         ylab=paste("PCoA2 (",pe2,"% var. explained)"))
par(adj = 0)
#title(main= "E")
par(adj=.5)
points(otus.p, 
       col= IBM[as.factor(metadat2$Treatment)],
       pch= c(16,8)[as.factor(metadat2$Measurement)],
       
       lwd=2,cex=1.2,
       bg=IBM[as.factor(metadat2$Treatment)],)

legend("bottomright", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=1,
       title = "",     bty = "o")


# permanova
adonis2(otus.bray ~ Treatment, data = metadat2)
adonis2(otus.bray ~ Treatment+Measurement+Treatment*Measurement, data = metadat2, by="terms")


#### CAP ##### 

# 2. Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(otus.bray ~ Treatment,
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
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/")
svg("cap.total.predicted.svg", width = 6 , height = 6)

#windows(6,6)
par(cex.lab = 1.8) # make all fonts in graphs little bigger
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"))
par(adj = 0)
#title(main= "E")
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
            cex=2)
 # legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
 #        fill= IBM,
 #        cex=1,
 #        bty = "n")
legend("topright", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=1.8,
      bty = "o")

dev.off()

### extra plot ####
ps1 <-subset_samples(ps, Fraction=="Total" & Treatment!="Soil" & Measurement!="predicted")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil" & Measurement!="predicted" )

#factor
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))


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

# load colors

IBM <- c( #IBM colors
  "navy", # dark royal blue
  "#648FFF", # french blue
  "#785EF0", # light purple
  "#DC267F", # magenta pink 
  "#FE6100", # bright orange
  "#FFB000", # golden yellow
  "#865338" # medium mocha brown
)

# PCoA 
par(adj=.5)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="",
         xlab=paste("PCoA1 (",pe1,"% var. explained)"),
         ylab=paste("PCoA2 (",pe2,"% var. explained)"))
par(adj = 0)
#title(main= "E")
par(adj=.5)
points(otus.p, 
       col= IBM[as.factor(metadat2$Treatment)],
       pch= c(16,8)[as.factor(metadat2$Measurement)],
       
       lwd=2,cex=1.2,
       bg=IBM[as.factor(metadat2$Treatment)],)

legend("bottomright", legend=c("measured", "predicted"  ),
       pch=c(16,8 ),
       cex=1,
       title = "",     bty = "o")


# permanova
adonis2(otus.bray ~ Treatment, data = metadat2)


#### CAP  

#  Run the CAP (db-RDA) analysis
# Formula: distance_matrix ~ environmental_variable_1 + environmental_variable_2
cap_result <- capscale(otus.bray ~ Treatment,
                       data = metadat2,
                       add = TRUE) # 'add = TRUE' handles negative eigenvalues from PCoA

anova.cca(cap_result, by="terms")

smry <- summary(cap_result)
smry
sc_si <- scores(cap_result, display="sites", choices=c(1,2), scaling=1)
sc_si[,2]<-sc_si[,2]*-1 # flip y axis 

# Extract the model's adjusted R2
RsquareAdj(cap_result)$adj.r.squared

# percent varience of total varience on RDA 1 and RDA 2
perc <- round(100*(summary(cap_result)$cont$importance[2, 1:2]), 2)
perc

### 4. plot 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/")
svg("cap.total.svg", width = 6 , height = 6)

#windows(6,6)
par(cex.lab = 1.8) # make all fonts in graphs little bigger
ordiplot(cap_result, choices=c(1,2), scaling =1, type="none",
         xlab=paste("CAP 1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("CAP 2 (",round(perc[2],2),"% variance explained)"),
         ylim = c(-.3, .3))
par(adj = 0)
#title(main= "E")
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
            cex=2)
 legend("topright", legend=c("L", "G", "B", "GB", "LB", "LG", "LGB"),
        fill= IBM,
        cex=1,
        bty = "n")
# legend("topright", legend=c("measured", "predicted"  ),
#        pch=c(16,8 ),
#        cex=1.8,
#        bty = "o")

dev.off()





######STATS ###########
# seperatated out by treatment
#GB
ps1 <-subset_samples(ps,  Treatment=="GB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 1745 taxa
# subset metadata
metadat2<-filter(metadat, Treatment=="GB")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2)


#LB
ps1 <-subset_samples(ps, Treatment=="LB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Treatment=="LB")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2)




#LG
ps1 <-subset_samples(ps, Treatment=="LG")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Treatment=="LG")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2)


#LGB
ps1 <-subset_samples(ps, Treatment=="LGB")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Treatment=="LGB")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# permanova
adonis2(otus.bray ~ Measurement, data = metadat2)









# compare pairwise distance ##############
ps1 <-subset_samples(ps, Fraction=="Total" & Treatment!="Soil" & n_species!="1" &Rep!="2" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil", n_species!="1", Rep!="2"  ) %>%
  select(Treatment, Measurement)

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
dist_mat<-as.matrix(otus.bray) 
dist_mat
key <-cbind(metadat2, mat)



# Define the Mapping
# We create a lookup table that links the Measured Sample ID to its Prediction ID
# Note: Ensure the order in 'measured_ids' matches the N1, N3, N4... order
mapping <- data.frame(
  Treatment = c(rep("GB", 5), rep("LB", 5), rep("LG", 5), rep("LGB", 5)),
  Prediction_ID = c(
    paste0("predict_GB_N", c(1,3,4,5,6)),
    paste0("predict_LB_N", c(1,3,4,5,6)),
    paste0("predict_LG_N", c(1,3,4,5,6)),
    paste0("predict_LGB_N", c(1,3,4,5,6))
  ),
  Measured_ID = c(
    "T_DNA_13_S130", "T_DNA_41_S90", "T_DNA_56_S168", "T_DNA_71_S159", "T_DNA_85_S139", # GB
    "T_DNA_12_S119", "T_DNA_42_S101", "T_DNA_57_S92", "T_DNA_72_S170", "T_DNA_86_S150", # LB
    "T_DNA_11_S108", "T_DNA_43_S112", "T_DNA_58_S103", "T_DNA_73_S94", "T_DNA_87_S161", # LG
    "T_DNA_14_S141",  "T_DNA_44_S123", "T_DNA_59_S114", "T_DNA_74_S105", "T_DNA_88_S172"  # LGB
  )
)

mapping
dist_mat
# Extract the Distances
# We loop through the mapping and pull the specific intersection from the matrix
mapping$BC_Distance <- mapply(function(p, m) dist_mat[p, m], 
                              mapping$Prediction_ID, 
                              mapping$Measured_ID)


# Statistical Comparison
# Test if the prediction error (distance) differs by treatment
fit <- aov(BC_Distance ~ Treatment, data = mapping)
summary(fit)

kruskal.test(BC_Distance ~ Treatment, data = mapping)

# 5. Visualization
ggplot(mapping, aes(x = Treatment, y = BC_Distance, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1)) +
  labs(title = "Within-Pair Prediction Accuracy",
       y = "Bray-Curtis Distance (Measured vs. Predicted)",
       x = "Treatment") +
  theme_bw()


######### compare pairwise distance -1 to 1 ##############
ps1 <-subset_samples(ps, Fraction=="Total" & Treatment!="Soil" &Rep!="2" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil",  Rep!="2"  ) %>%
  select(Treatment, Measurement)

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps1), method = "bray")
dist_mat<-as.matrix(otus.bray) 
dist_mat
key <-cbind(metadat2, dist_mat)
key


# Define the Mapping
# We create a lookup table that links the Measured Sample ID to its Prediction ID
# Note: Ensure the order in 'measured_ids' matches the N1, N3, N4... order
mapping <- data.frame(
  Treatment = c(rep("GB", 5), rep("LB", 5), rep("LG", 5), rep("LGB", 5)),
  Prediction_ID = c(
    paste0("predict_GB_N", c(1,3,4,5,6)),
    paste0("predict_LB_N", c(1,3,4,5,6)),
    paste0("predict_LG_N", c(1,3,4,5,6)),
    paste0("predict_LGB_N", c(1,3,4,5,6))
  ),
  Measured_ID = c(
    "T_DNA_13_S130", "T_DNA_41_S90", "T_DNA_56_S168", "T_DNA_71_S159", "T_DNA_85_S139", # GB
    "T_DNA_12_S119", "T_DNA_42_S101", "T_DNA_57_S92", "T_DNA_72_S170", "T_DNA_86_S150", # LB
    "T_DNA_11_S108", "T_DNA_43_S112", "T_DNA_58_S103", "T_DNA_73_S94", "T_DNA_87_S161", # LG
    "T_DNA_14_S141",  "T_DNA_44_S123", "T_DNA_59_S114", "T_DNA_74_S105", "T_DNA_88_S172"  # LGB
  )
)

mapping
dist_mat
# Extract the Distances
# We loop through the mapping and pull the specific intersection from the matrix
mapping$BC_Distance <- mapply(function(p, m) dist_mat[p, m], 
                              mapping$Prediction_ID, 
                              mapping$Measured_ID)


# Statistical Comparison
# Test if the prediction error (distance) differs by treatment
fit <- aov(BC_Distance ~ Treatment, data = mapping)
summary(fit)

kruskal.test(BC_Distance ~ Treatment, data = mapping)

# 5. Visualization
ggplot(mapping, aes(x = Treatment, y = BC_Distance, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1)) +
  labs(title = "Within-Pair Prediction Accuracy",
       y = "Bray-Curtis Distance (Measured vs. Predicted)",
       x = "Treatment") +
  theme_bw()

