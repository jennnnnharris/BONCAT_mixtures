# kayla coode
# Color palette
mycolors <- c("#855C75FF", "#D9AF6BFF", "#AF6458FF", "#736F4CFF", "#526A83FF", "#625377FF", "#68855CFF", "#9C9C5EFF", "#A06177FF", "#C38961FF", "#467378FF", "#7C7C7CFF", "#586028FF", "#755028FF", "#006475FF", "#8C785DFF)

# Across-sample thresholding: 3x25 (throw out "non-reproducible" ASVs)
threshold <- kOverA(3,A=25) # set threshold values (require k samples with A reads)
phy.nobadASVs.thresholded <- filter_taxa(phy.nobadASVs,threshold,TRUE)

# What proportion of original reads remain?
sum(taxa_sums(phy.nobadASVs.thresholded))/(sum(taxa_sums(phy.nobadASVs)))

# Remove samples with <300 usable reads
phy.nobadASVs.thresholded.highcoverage <- subset_samples(phy.nobadASVs.thresholded,UsableReads>=300) %>% prune_taxa(taxa_sums(.)>0,.)

# Constrained ordination
p1.cap <- ordinate(p1.ccmAim3.asv.clr,method='CAP',distance='bray',formula=~CommonCoverCropCode*InoculumType+Condition(log(UsableReads)))

# Visualize ordination
plot_ordination(p1.ccmAim3.asv.clr,p1.cap,color="CommonCoverCropCode")+
  facet_wrap(~InoculumType)+
  theme_bw()+
  geom_point(size=2.5)+
  stat_ellipse(aes(group=CommonCoverCropCode), linetype=2)+
  theme(text=element_text(size=15), strip.text.x=element_text(size=15.5),
        legend.position="top")+
  scale_color_manual(name="Cover crop", labels=c(BC='Bush clover', CP='Cowpea', 
                                                 WC='White clover', CC='Crimson clover', 
                                                 FP='Field pea', YSC='Yellow sweet clover'))



                                                 # Transform ASV counts to proportions / relative abundance
p1.ccmAim3.nod.prop <- transform_sample_counts(subset_samples(p1.ccmAim3.nod.asv),function(ASV) ASV/sum(ASV))

# Subset by family
p1.ccmAim3.nod.prop.family <- phyloseq::tax_glom(p1.ccmAim3.nod.prop, "Family") # 196 taxa

# Take x number of taxa per taxonomic rank based on relative abundance
# Taxa > n will be added to as "other" label
p1.ccmAim3.nod.prop.family.top10 <- fantaxtic::get_top_taxa(physeq_obj = p1.ccmAim3.nod.prop.family, n=10, relative = TRUE, discard_other = FALSE, other_label = "Other")

# Melt data frame with phyloseq function for plotting
p1.family.nod.top10 <- psmelt(p1.ccmAim3.nod.prop.family.top10)

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