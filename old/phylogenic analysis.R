# phylogenic analysis for BONCAT gradients
# cut from manuscript Nov 2024
# Jennifer Harris

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



