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


# Load required libraries #
library(tidyverse)
library(vegan)
library(readxl)
library(lubridate)
library(phyloseq)
library(MicEco)
library(multcompView)
library(paletteer)


# colors
mycols7<-c( "#715b8a", "#4D8F8BFF", "#CDD6ADFF", "#365C83FF", "#AD5A6BFF", "#E3C1CBFF",  "#384351FF")
mycols7v<-c( "#715b8a", "#4D8F8BFF", "#b1de64", "#365C83FF", "#bd0262", "#c77597",  "#384351FF")
mycols4 <- c("#715b8a",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )



# set shapes
myshapes <- c(1, 12,15 ,21, 22, 23 , 24)
myshapes2 <- c(21 , 12, 24,1, 15 , 22, 23 )



  
#import data#
# Set the working directory 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/16S_sequencing/")

taxon <- read.csv("all/taxonomy.csv", header=T)
asvs <- read.table("all/feature.table.tsv", sep="\t", header=T, row.names = 1)
metadat<-read.csv("metadat.csv", header = T)

## Transpose ASVS table ##
asvs[1:5,1:5]#taxa are columns
asvs<-t(asvs)

#T_DNA_23_S153 has really few reads so I am omitting it.
asvs<-asvs[which(row.names(asvs)!= "T_DNA_23_S153"),]

#get min number of reads in a sample
min.s<-min(rowSums(asvs))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs<-rrarefy(asvs, min.s)

# order metadata
metadat<-as.data.frame(metadat[order(metadat$SampleID),])
metadat$Treatment   <- factor(metadat$Treatment, levels= c("Soil", "L", "G", "B", "GB", "LB", "LG", "LGB"))
row.names(metadat) <- metadat$SampleID
# make composition column


# make composition var
metadat$composition <- metadat$Treatment
metadat$composition<-gsub("L", "Legume", metadat$composition)
metadat$composition<-gsub("G", "Grass", metadat$composition)
metadat$composition<-gsub("B", "Brassica", metadat$composition)
metadat$composition<-gsub("LegumeBrassica", "Legume_Brassica", metadat$composition)
metadat$composition<-gsub("LegumeGrass", "Legume_Grass", metadat$composition)
metadat$composition<-gsub("GrassBrassica", "Grass_Brassica", metadat$composition)
metadat$composition 
metadat$composition<- factor(metadat$composition, levels = c("Soil", "Legume", "Grass", "Brassica", "Grass_Brassica", "Legume_Brassica", "Legume_Grass", "Legume_Grass_Brassica"))


#make taxon matrix row names OTUs
row.names(taxon) <- taxon$Feature.ID

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
# 218 k asvs



####### 2. DIVERSITY ####

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

#BCAT_73_S35 is kind of a weird outlier
#rich<-rich[which(rich$SampleID!="BCAT_73_S35"),]


##plots##

# nitrogen effect
rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=as.factor(N), y=Shannon))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  geom_jitter(size=1.5)+
  theme_classic(base_size = 16)
  

# Treatment effect?
# shannon diversity letters
rich<-rich%>%  filter(Fraction=="Total") %>% filter(Treatment!="Soil")
  
lab<-as.character(rich$Treatment)
lab<- gsub("LGB", "X", lab)
lab <- gsub("LG", "B", lab)
lab <- gsub("LB", "B", lab)
lab <- gsub("GB", "B", lab)
lab
lab <- gsub("L", "Y", lab) 
lab <- gsub("G", "B", lab) 
lab <- gsub("Y", "A", lab) 
lab <- gsub("X", "AB", lab) 
lab

p1<-rich%>%  
#filter(N=="0") %>%
  ggplot(aes(x=Treatment, y=Shannon,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  #geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  geom_jitter(size=1.5, width=.1)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Total DNA Shannon Diversity")+
 geom_text(y=7.9, label = lab, size=5)+
  xlab("")
  #scale_shape_discrete() 
p1

#label for observed

lab<-as.character(rich$Treatment)
lab<- gsub("LGB", "X", lab)
lab<- gsub("LB", "Y", lab)
lab<- gsub("GB", "Y", lab)
lab<- gsub("LG", "Y", lab)
lab<- gsub("B", "Y", lab)
lab<- gsub("G", "X", lab)
lab<- gsub("L", "Z", lab)
#sub letters
lab<- gsub("X", "AB", lab)
lab<- gsub("Y", "B", lab)
lab<- gsub("Z", "A", lab)


p2<- rich%>% filter(Fraction=="Total") %>% filter(Treatment!="Soil") %>%
  ggplot(aes(x=composition, y=Observed,  fill=Treatment))+
  geom_boxplot(alpha=.5, outlier.shape = NA) +
  scale_color_manual(values=mycols7) +
  scale_fill_manual(values = mycols7)+
  #geom_jitter(aes(shape = as.factor(Rep) ), width = .1, size=2,  )+
  geom_jitter(width = .1,size=1.5)+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5),legend.position="none")+
  ylab("Total DNA Number of Asvs")+
  xlab("")+
  geom_text(y= 5100, label = lab, size=5)
  
#scale_shape_discrete() 
p2
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig_diversity")
svg(file="diversity.total.svg",width = 10, height=6)
require(gridExtra)
#windows(10,4)
grid.arrange(p1, p2, ncol=2)
dev.off()



#STATS#


# remove soil
rich <- rich %>% filter(Fraction=="Total") %>% filter(Treatment!="Soil")
rich$Treatment   <- factor(rich$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))
rich$N <- factor(rich$N, levels= c("0", "1"))


# Shannon
rich1<-rich%>%  
  filter(N=="1")
m1<-lm(Shannon ~ N,  data = rich1)
summary(m1)

m1<-lm(Shannon ~ Treatment*N,  data = rich)
summary(m1)
# no interaction with nitrogen, so I'm just gonna bin +N =N together

# Analysis of variance 
anova1<- aov(Shannon ~ Treatment, data = rich)
summary(anova1)
library(multcompView)
tukey.a1 <- TukeyHSD(anova1)

print(tukey.a1) # all difference except LG-LB and LGB-LG
plot(anova1) #homoscedasticity looks fine

##observed ANOVA

# #chao1 Analysis of variance 
anova1<- aov(Observed ~ Treatment, data = rich)
summary(anova1)
# 3. Run the Tukey HSD test on the ANOVA model
tukey_output <- TukeyHSD(anova1)

# The 'tukey_output' shows the pairwise comparisons and p-values
print(tukey_output)

# 4. Extract the p-values for the factor of interest
p_values <- tukey_output$Treatment[, 4]

# 5. Generate the grouping letters using multcompLetters()
# The 'multcompLetters()' function takes a named vector of p-values
cld <- multcompLetters(p_values)

# 6. Print the results
print(cld)



##chao1 Analysis of variance 
anova1<- aov(Chao1 ~ Treatment, data = rich)
summary(anova1)
# 3. Run the Tukey HSD test on the ANOVA model
tukey_output <- TukeyHSD(anova1)

# The 'tukey_output' shows the pairwise comparisons and p-values
print(tukey_output)

# 4. Extract the p-values for the factor of interest
p_values <- tukey_output$Treatment[, 4]

# 5. Generate the grouping letters using multcompLetters()
# The 'multcompLetters()' function takes a named vector of p-values
cld <- multcompLetters(p_values)

# 6. Print the results
print(cld)

######## 3. Filtering rare taxa ##################

# remove true singletons 
ps<-prune_taxa(taxa_sums(ps) > 1, ps)
ps
# 183095 

#remove asvs with a mean of less than 5
mean.reads <- rowSums(t(otu_table(ps)))/nsamples(ps)
keep<-row.names(t(otu_table(ps))[ mean.reads > 5, ])
ps<-prune_taxa(keep, ps)
ps
#1696 asvs

# remove asvs that are in less than 5 samples
#ps <-ps_prune(ps, min.samples = 3) # no features to group!
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
#1696 asvs

## plot
#plot(sort(taxa_sums(total), TRUE), type="h", ylim=c(0, 8000))


######## 4. CAP ##################

# Constrained ordination
# Perform vegdist analysis of BC distances #
ps1 <-subset_samples(ps, Treatment !="Soil" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
otus.bray<-vegdist(otu_table(ps1), method = "bray")
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Treatment!="Soil")
#set factors 
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

#overall L+B+G
#p1.cap <- ordinate(ps1, method='CAP',distance='bray',formula=~Grass*Legume*Brassicae)
p1.cap <- ordinate(ps1, method='CAP',distance='bray',formula=~N*Treatment)
anova.cca(p1.cap, by="terms")



p1.cap


#col
# colors
mycols7<-c( "#715b8a", "#4D8F8BFF", "#a5c76f", "#365C83FF", "#bd0262", "#c77597",  "#384351FF")
mycols7<-c( "#715b8a", "#4D8F8BFF", "#a5c76f", "#365C83FF", "#bd0262", "#c77597",  "white")

#mycols7<-c( "#715b8a", "#4D8F8BFF", "#a5c76f", "white", "white", "white",  "#384351FF")
mycols4 <- c("#715b8a",  "#AD5A6BFF", "#E3C1CBFF", "#384351FF" )


# cap plot total
p1 <-plot_ordination(ps1,  p1.cap, color="composition")+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
        legend.position="left")+
  scale_color_manual(values =  mycols7, name="composition")+
  scale_shape_discrete(name= "Nitrogen")
#ggtitle("A  All cover crops")
p1

pathfig4 <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/Fig4_CAPtotal"
svg("cap.total.treatment3.svg", width = 8 , height = 4)
p1
dev.off()

# Legume present# #########

# subset data
ps1 <-subset_samples(ps, Legume == "1")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Legume == "1" )
#set factors 
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c("L", "LB", "LG", "LGB"))

#cap
p1.cap <- ordinate(ps1, method='CAP',distance='bray',formula=~Treatment)
anova.cca(p1.cap, by="terms")

p2<-plot_ordination(ps1,  p1.cap, color="Treatment", )+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  
  #stat_ellipse(    geom = "polygon",aes(fill = Treatment),alpha = 0.25,level = 0.95 )+  
  theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
        legend.position="none")+
  scale_color_manual(values =  mycols4, name="Treatment")+
  scale_fill_manual(values =  mycols4, name="Treatment")+
    
  scale_shape_discrete(name= "Nitrogen")+
  
  ggtitle("C Legume")
p2





# Grass# #####
# subset  data
ps1 <-subset_samples(ps, Grass=="1")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
metadat2<-filter(metadat, Fraction=="Total" & Grass=="1"  )
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "G", "GB",  "LG", "LGB"))

#cap
p1.cap <- ordinate(ps1, method='CAP',distance='bray',formula=~Treatment)
anova.cca(p1.cap, by="terms")


mycols4<- c( "#443983", "#21918c", "#90d743",  "#f7b307")

p3<-plot_ordination(ps1,  p1.cap, color="Treatment", )+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  
  #stat_ellipse(    geom = "polygon",aes(fill = Treatment),alpha = 0.25,level = 0.95 )+  
  theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
        legend.position="none")+
  scale_color_manual(values =  mycols4, name="Treatment")+
  scale_fill_manual(values =  mycols4, name="Treatment")+
  
  scale_shape_discrete(name= "Nitrogen") +
  ggtitle("C Grass")
p3


# Brass ###########

# subset data
ps1 <-subset_samples(ps, Brassicae=="1")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
metadat2<-filter(metadat, Fraction=="Total" & Brassicae=="1" )
metadat2$Treatment   <- factor(metadat2$Treatment, levels= c( "B", "GB", "LB",  "LGB"))

# cap
p1.cap <- ordinate(ps1, method='CAP',distance='bray',formula=~Treatment)
anova.cca(p1.cap, by="terms")

#Brass colors

mycols4<- c( "#31688e", "#21918c", "#35b779",  "#f7b307")


p4<-plot_ordination(ps1,  p1.cap, color="Treatment", )+
  theme_bw()+
  geom_point(aes(shape = as.factor(N) ), size=2.5)+
  stat_ellipse(aes(group=Treatment), linetype=1)+
  
  #stat_ellipse(    geom = "polygon",aes(fill = Treatment),alpha = 0.25,level = 0.95 )+  
  theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
        legend.position="none")+
  scale_color_manual(values =  mycols4, name="Treatment")+
  scale_fill_manual(values =  mycols4, name="Treatment")+
  
  scale_shape_discrete(name= "Nitrogen")  +
  ggtitle("D Brassicae")
p4




setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures/fig4_CAP")
svg(file="CAP.total.svg",width = 6, height=2.5)
require(gridExtra)
#windows(10,10)
grid.arrange( p2, p3, p4, ncol=3)
dev.off()







######## 5. PCOA ##################

#set shapes and cols
  myshapes <- c(1, 12,15 ,21, 22, 23 , 24)
  myshapes2 <- c(21 , 12, 24,1, 15 , 22, 23 )
  
  
######PCOA between fractions##
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
  
  
  
   
  
######PCOA of total ###

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


######PCOA nitrogen + no nitrogen###

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



######PCOA PLOTS  facet by number of species##

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


######PCOA PLOTS by group#

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



######PCOA STATS: BETA DISPERSION#
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

######PERMANOVA ######

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

######PERMANOVA by Treatment##

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



######PERMANOVA by Fraction###
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


######## 6. PERCENT abundance figure ##################
#get data
#otu_table(total)[1:5, 1:5]

# Transform ASV counts to proportions / relative abundance

total <-subset_samples(total, Treatment !="Soil" )
total<-prune_taxa(taxa_sums(total) > 0, total)
total

total.prop <- transform_sample_counts(subset_samples(total),function(ASV) ASV/sum(ASV))
otu_table(total.prop)[1:5, 1:5]

# remove taxa groups finer than Phyla
total.prop <- phyloseq::tax_glom(total.prop, "Phyla")
# check number of Phylas and names
length(unique(as.data.frame(tax_table(total.prop))$Phyla))
#
#head(as.data.frame(tax_table(total.prop)))
# 23 Phyla

# Melt data frame with phyloseq function for plotting
df <- psmelt(total.prop)
head(df)
#
#aggregate by Phyla
df <- df %>% group_by(across(c(-OTU, -Feature.ID))) %>%
summarise(Abundance= sum(Abundance))
#
#rename low abundance taxa as other
df.mean<-df %>%
group_by(Phyla) %>%
summarise(mean=mean(Abundance)) %>%
arrange(.,mean)
df.mean

remove<-filter(df.mean, mean<.0005 )$Phyla
#for anything in this list, remove it and call it other

for (i in remove) {
    print(i)
    df$Phyla<-sub(i, "other", df$Phyla)
  
    }

#df$Phy <- "other"
# check Phyla names
unique(df$Phyla)
# check number of Phylas
length(unique(df$Phyla))
#18
# remove the o__ in front of Phyla
df$Phyla<-sub(" p__", "", df$Phyla)
unique(df$Phyla)
# Reorder levels to put other at the end - otherwise taxa are in alphabetical Phyla
#install.packages("forcats")
library(forcats)
df$Phyla<-fct_relevel(as.factor(df$Phyla), "other", after = Inf)
# make really low abundance taxa other
unique(df$Phyla)
colnames(df)

df$Treatment <- factor(df$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))



#set working directory for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Figures")

# Visualize

svg("percent.abund.total.svg", width = 12, height = 7)
windows(12,7)
ggplot(df, aes(x=Sample, y=Abundance, fill=Phyla))+
  geom_bar(stat = "identity")+
  facet_wrap(~Treatment, scales="free_x", nrow = 1)+
  #           labeller=as_labeller(c(BC='Bush clover', YSC='Yellow sweet clover',
  #                                  FP='Field pea', CC='Crimson clover',
  #                                  WC='White clover', CP='Cowpea')))+
  #labs(y="Relative abundance", title="Nodule")+
  theme_classic()+
  theme(legend.position='right', axis.title.x=element_blank(),
        axis.text.x=element_blank(), text=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=16), legend.text=element_text
        (size=14))+
  ggtitle("Total DNA")+ 
  scale_fill_manual(values=c("#855C75", "#dba7c7",
                             "#996e29", "#D9AF6B",
                             "#AF6458", "#eba69b" , 
                             "#736F4C", "#bdb993", 
                             "#455669", "#8095ab",
                             "#625377", "#baafc9",
                             "#68855C", "#c7dbbf", 
                             #"#9C9C5E", "#e3e3bc",
                             "#467378", "#b9c8c9",
                             "#8C785D", "#ded6cc",
                             "#7C7C7C", "#C3C3C3"))
dev.off()



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




