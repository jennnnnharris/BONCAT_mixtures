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
# Install and load
#install.packages("compositions")
library(compositions)
#BiocManager::install
#library(multcompView)
#library(BiodiversityR)


#import data ###############
# Set the working directory 
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
#import data
taxon <- read.csv("16S_taxonomy.csv", header=T)
asvs <- read.table("16S_feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("16S_metadata.csv", header = T, row.names = 1)
# Transpose ASVS table #
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

# order metadata #
metadat
metadat<-as.data.frame(metadat[order(row.names(metadat)),])

metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
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
biomass<-read.csv("biomass_percent.csv")
# filter for only n+ and required columns
biomass<-biomass %>% filter(N==1) %>% select(-Total.Root.g, -Total.withbulk, -Stem.Biomass.g, -Root.Biomass.g, -Bulk.Root.g)

# Note: interpolate the mean is for missing biomass observation LB mixture rep 5 #




# predict for LG 

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

# predict for GB 

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


# predict for LB 
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



# predict for LBG
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
write.csv(otus, "total.feature.table_predicted16S.csv")

#taxon
tax_table<-as.data.frame(tax_table(ps))
write.csv(tax_table, "total.taxonomy_predicted16S.csv.csv")



##### import predicted #####
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/Data_for_upload")
taxon <- read.csv("total.taxonomy_predicted16S.csv", row.names = 1)
asvs <- read.csv("total.feature.table_predicted16S.csv", row.names = 1)
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
# make col
metadat$Trt_measure <- paste0(metadat$Treatment, metadat$Measurement)
metadat$Trt_measure<-gsub("measured", "", metadat$Trt_measure)
head(metadat)

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps


##### remove soil
ps <-subset_samples(ps, Fraction=="Total" & Treatment!="Soil" )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps

# subset metadata
metadat<-filter(metadat, Fraction=="Total" & Treatment!="Soil" )

#factor
metadat$Treatment   <- factor(metadat$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))


# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps), method = "bray")




########## Aitchison distance PLOT  #####
ps1 <-subset_samples(ps, Fraction=="Total" &  Rep!=2 & Treatment!="Soil"  )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" &  Rep!=2 & Treatment!="Soil") %>%
  select(Treatment, Measurement, Rep)
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

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
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_predicted_total")
#svg("total_predicted_PCA.svg",  width=6, height=6)

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
 
#dev.off() 
 
# 
# 
# 
# 
# ######pairwise adonis test
# # seperatated out by treatment
# #GB
# ps1 <-subset_samples(ps,  Treatment=="GB")
# ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
# ps1
# # 1745 taxa
# # subset metadata
# metadat2<-filter(metadat, Treatment=="GB")
# # Calculate Bray-Curtis distance between samples
# otus.bray<-vegdist(otu_table(ps1), method = "bray")
# # Perform PCoA analysis of BC distances #
# otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# # permanova
# adonis2(otus.bray ~ Measurement, data = metadat2)
# 
# 
# #LB
# ps1 <-subset_samples(ps, Treatment=="LB")
# ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
# ps1
# # subset metadata
# metadat2<-filter(metadat, Treatment=="LB")
# # Calculate Bray-Curtis distance between samples
# otus.bray<-vegdist(otu_table(ps1), method = "bray")
# # Perform PCoA analysis of BC distances #
# otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# # permanova
# adonis2(otus.bray ~ Measurement, data = metadat2)
# 
# 
# 
# 
# #LG
# ps1 <-subset_samples(ps, Treatment=="LG")
# ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
# ps1
# # subset metadata
# metadat2<-filter(metadat, Treatment=="LG")
# # Calculate Bray-Curtis distance between samples
# otus.bray<-vegdist(otu_table(ps1), method = "bray")
# # Perform PCoA analysis of BC distances #
# otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# # permanova
# adonis2(otus.bray ~ Measurement, data = metadat2)
# 
# 
# #LGB
# ps1 <-subset_samples(ps, Treatment=="LGB")
# ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
# ps1
# # subset metadata
# metadat2<-filter(metadat, Treatment=="LGB")
# # Calculate Bray-Curtis distance between samples
# otus.bray<-vegdist(otu_table(ps1), method = "bray")
# # Perform PCoA analysis of BC distances #
# otus.pcoa <- cmdscale(otus.bray, k=(3-1), eig=TRUE)
# # permanova
# adonis2(otus.bray ~ Measurement, data = metadat2)
# 


 ###scalar projection GB #####
 
 # clr matrix
 ps1 <-subset_samples(ps, Fraction=="Total" &  Rep!=2 & Treatment!="Soil"  )
 ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
 ps1
 
 # subset metadata
 metadat2<-filter(metadat, Fraction=="Total" &  Rep!=2 & Treatment!="Soil") %>%
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
mix_data$rep <- rep(c("1","3", "4", "5", "6"), 2)
 

# # plot raw -1 to 1
 ggplot(mix_data, aes(x = Measurement, y = GB_index, fill = Measurement)) +
   geom_boxplot(alpha = 0.6, outlier.shape = NA) +
   scale_fill_manual(values = c("#DC267F", "grey70"))+
   geom_jitter(width = 0.1, size = 2) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
   annotate("text", x = 0.5, y = -1, label = "Pure G", hjust = 0) +
   annotate("text", x = 0.5, y = 1, label = "Pure B", hjust = 0) +
   labs(title = "Mixture Composition Along G-to-B Axis",
        y = "Projection Index (G → B)",
        x = "Treatment") +
   theme_minimal()+
   coord_flip()
### plot with predicted as zero ###
 # subtract each rep 
 mix_data<-mix_data %>% select(GB_index, Treatment, Measurement)
 P<-rep(mix_data$GB_index[which(mix_data$Measurement == "predicted")],2)
 mix_data$diff_predict <- mix_data$GB_index-P 
 
# plot
#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_index")
#svg("GB_total.svg",  width=4, height=1.5)
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
#dev.off()  

 
#stats #############
mix_data


mix_data$rep <- rep(c("1","3", "4", "5", "6"), 2)
# run anova
a1<-aov(GB_index ~ Measurement+rep, data = mix_data)
summary(a1)

# paired test 
p<-mix_data$GB_index[which(mix_data$Measurement=="predicted")]
m<-mix_data$GB_index[which(mix_data$Measurement!="predicted")]
paired_test<-t.test(p,m, paired = TRUE, data = mix_data)



 
 ###scalar projection LB #####
 
 # clr matrix
 ps1 <-subset_samples(ps, Fraction=="Total" &  Rep!=2 & Treatment!="Soil"  )
 ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
 ps1

 # subset metadata
 metadat2<-filter(metadat, Fraction=="Total" &  Rep!=2 & Treatment!="Soil") %>%
   select(Treatment, Measurement, Rep)
 metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

 # CLR transform
 clr_data <- clr(otu_table(ps1)+1)
 clr_data<-as.matrix(clr_data)
 clr_data[1:5, 1:5]

 # 1. Calculate the centroids (mean vector) for your baseline monocultures
 # (Assuming 'clr_matrix' contains only your numeric CLR-transformed columns)
 centroid_L <- colMeans(clr_data[ which(metadat2$Treatment=="L"), ])
 centroid_B <- colMeans(clr_data[which(metadat2$Treatment=="B"), ])
 
 # 2. Define the axis vector from L to B
 v <- centroid_B - centroid_L
 v_length_sq <- sum(v^2)
 

 # 4. Apply this to your target samples (LB observed and LB predicted)
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
 
   
 #### plot with predicted as zero ###
  # subtract each rep #
 mix_data<-mix_data %>% select(LB_index, Treatment, Measurement)
 P<-rep(mix_data$LB_index[which(mix_data$Measurement == "predicted")],2)
 P
 mix_data$diff_predict <- mix_data$LB_index-P 
 mix_data
 
 # plot
 #setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_predicted_total")
 #svg("LB_total.svg",  width=4, height=1.5)
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
   
   
   # stats #############
   
   mix_data     
   mix_data$rep <- rep(c("1","3", "4", "5", "6"), 2)
   
   # run anova
   a1<-aov(LB_index ~ Measurement+rep, data = mix_data)
   summary(a1)
   
   
   # paired test 
   p<-mix_data$LB_index[which(mix_data$Measurement=="predicted")]
   m<-mix_data$LB_index[which(mix_data$Measurement!="predicted")]
   paired_test<-t.test(p,m, paired = TRUE, data = mix_data)
   paired_test
   
   
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
 mix_data<-mix_data %>% select(index, Treatment, Measurement)
 P<-rep(mix_data$index[which(mix_data$Measurement == "predicted")],2)
 mix_data$diff_predict <- mix_data$index-P 
 mix_data
 

 
 #setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_pt")
 #svg("LG_total.svg", width=4, height=1.5)
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
 
 
 ###stats #############
 
 mix_data
 mix_data$rep <- rep(c("1","3", "4", "5", "6"), 2)
 # run anova
 a1<-aov(index ~ Measurement+rep, data = mix_data)
 summary(a1)
 
  # paired test 
 p<-mix_data$index[which(mix_data$Measurement=="predicted")]
 m<-mix_data$index[which(mix_data$Measurement!="predicted")]
 paired_test<-t.test(p,m, paired = TRUE, data = mix_data)
 paired_test
 
 
 
 ######simple ordination three species ####
 
 # LGB, G, B, G example
 ps1 <-subset_samples(ps, Fraction=="Total" &  Rep!=2) %>% subset_samples(Treatment=="L" | Treatment=="G" | Treatment=="B" | Treatment=="LGB")    
 ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
 ps1
 
 # subset metadata
 metadat2<-filter(metadat, Fraction=="Total" &  Rep!=2) %>% filter(Treatment=="L" | Treatment=="G" | Treatment=="B" | Treatment=="LGB") %>%
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

 
 ##### scalar projection LGB ########
 

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
 
 # 4. Apply this to your target samples (LGB observed and LGB predicted)
 # Let's say you target a data frame of your mixtures
 mix_data <- clr_data[which(metadat2$Treatment=="LGB"), ]
 mix_data[1:5, 1:5]
 dim(mix_data)
  
 # Calculate the index for each row
 
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
 
 
 #### plot with predicted as zero ###
 # subtract each rep #
 mix_data<-mix_data %>% select(LG_index, GB_index, BL_index, Treatment, Measurement)
 P<-rep(mix_data$GB_index[which(mix_data$Measurement == "predicted")],2)
 mix_data$adj_GBindex <- mix_data$GB_index-P 
 mix_data
 P<-rep(mix_data$LG_index[which(mix_data$Measurement == "predicted")],2)
 mix_data$adj_LGindex <- mix_data$LG_index-P 
 mix_data
 P<-rep(mix_data$BL_index[which(mix_data$Measurement == "predicted")],2)
 mix_data$adj_BLindex <- mix_data$BL_index-P 
 mix_data
 
 
 mix_data %>% filter(Measurement=="measured") %>%
   ggplot(aes(x = Treatment, y = adj_LGindex, fill = Measurement)) +
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
 
 
 
 mix_data %>% filter(Measurement=="measured") %>%
   ggplot(aes(x = Treatment, y = adj_BLindex, fill = Measurement)) +
   geom_boxplot(alpha = 0.6, outlier.shape = NA) +
   geom_jitter(width = 0.1, size = 2) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
   annotate("text", x = 0.7, y = -1, label = "B", hjust = 0) +
   annotate("text", x = 0.7, y = 1, label = "L", hjust = 0) +
   labs(title = "",
        y = "Projection Index (L → G)",
        x = "Mixture") +
   scale_fill_manual(values = c("#FFB000", "grey70"))+
   coord_flip()+
   theme_minimal()+
   theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
 
 mix_data %>% filter(Measurement=="measured") %>%
   ggplot(aes(x = Treatment, y = adj_GBindex, fill = Measurement)) +
   geom_boxplot(alpha = 0.6, outlier.shape = NA) +
   geom_jitter(width = 0.1, size = 2) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
   annotate("text", x = 0.7, y = -1, label = "G", hjust = 0) +
   annotate("text", x = 0.7, y = 1, label = "B", hjust = 0) +
   labs(title = "",
        y = "Projection Index (L → G)",
        x = "Mixture") +
   scale_fill_manual(values = c("#FFB000", "grey70"))+
   coord_flip()+
   theme_minimal()+
   theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

 
 
 ### triplot ######
 
 
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
 
 #### plot with predicted as zero ###
 # # subtract each rep #
 # P<-rep(mix_data$GB_index[which(mix_data$Measurement == "predicted")],2)
 # mix_data$adj_GBindex <- mix_data$GB_index-P 
 # mix_data
 # P<-rep(mix_data$LG_index[which(mix_data$Measurement == "predicted")],2)
 # mix_data$adj_LGindex <- mix_data$LG_index-P 
 # mix_data
 # P<-rep(mix_data$BL_index[which(mix_data$Measurement == "predicted")],2)
 # mix_data$adj_BLindex <- mix_data$BL_index-P 
 mix_data
  
 #stats #
 
 mix_data
 
 # Run an independent t-test (Welch's t-test is default, which handles unequal variance safely)
 t_test_result <- t.test(GB_index ~ Measurement, data = mix_data, alternative = "two.sided")
  print(t_test_result)
 
  t_test_result <- t.test(LB_index ~ Measurement, data = mix_data, alternative = "two.sided")
  print(t_test_result)
 
 t_test_result <- t.test(LG_index ~ Measurement, data = mix_data, alternative = "two.sided")
 print(t_test_result)
 
 
 library(ggplot2)
 library(ggtern)
 
 # View one of the triads
 # Use ggtern which is a specialized package using ggplot2 tools.
 # You will define the data source (demodata1) and then provide the three column names where ggtern
 # will find the x, y, and z coordinate data (T1A, T1B, and T1C)
 # Then you tell ggtern what kind of geometry to use to display the data, in this case "geom_point"


 p4<-ggtern(data=mix_data, aes(x=LG_index, y=GB_index, z=BL_index, colour = Measurement, shape = factor(Rep))) +
   geom_point(size=2)+
   theme_minimal()+
   scale_color_manual(values = c( "#865338", "grey70"), labels=c("measured", "expectation"))+
   scale_shape_manual(values=c(8, 15, 17, 19, 9), name="Rep")+
   xlab("L")  +                  
   ylab("G") +
   zlab("B")   

#setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig_index")
svg("LGB_total.svg", width=4, height=4) 
 p4  
 dev.off()
 
 
 
 