#Pairwise adonis tables

#Note, Combined_16S = my phyloseq object

library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pairwiseAdonis)

# 1) Subset phyloseq object
# subset data
ps1 <-subset_samples(ps, Treatment !="Soil" & N=="0" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# subset metadata
meta<-filter(metadat, Fraction=="Total" & Treatment!="Soil"  & N=="0")
meta$Treatment   <- factor(meta$Treatment, levels= c("L", "G", "B", "GB", "LB", "LG", "LGB"))

# bray curtis 
bray<-vegdist(otu_table(ps1), method = "bray")

# Quick checks
unique(meta$Treatment)

# 2) Global PERMANOVA
global_adonis <- adonis2(
  bray ~ Treatment,
  data = meta,
  strata = meta$block,
  by = "terms"
)

global_adonis

# 3) Pairwise PERMANOVA and heatmap (Treatment term only)
pw <- pairwise.adonis2(
  bray ~ Treatment,
  data = meta,
  p.adjust.m = "fdr",
  by = "terms",
  perm = 999
)

treatment_pw <- bind_rows(lapply(names(pw), function(x) {
  as.data.frame(pw[[x]]) %>%
    rownames_to_column("Term") %>%
    mutate(Comparison = x)
})) %>%
  filter(Term == "Treatment") %>%
  separate(Comparison, into = c("Group1", "Group2"), sep = "_vs_") %>%
  mutate(
    `P-value` = as.numeric(`Pr(>F)`),
    Label = sprintf("R2 = %.2f", R2)
  ) %>%
  filter(Group1 != Group2) %>%
  mutate(
    # keep one triangle only (avoids duplicate mirrored pairs)
    Group1 = as.character(Group1),
    Group2 = as.character(Group2)
  ) %>%
  rowwise() %>%
  mutate(
    g_min = min(c(Group1, Group2)),
    g_max = max(c(Group1, Group2))
  ) %>%
  ungroup() %>%
  transmute(Group1 = g_min, Group2 = g_max, `P-value`, Label)

treatment_pw

#treatment_pw$Group1<-factor(treatment_pw$Group1,levels= c("L", "G", "B", "GB", "LB", "LG", "LGB") )
#treatment_pw$Group2<-factor(treatment_pw$Group2,levels= c("L", "G", "B", "GB", "LB", "LG", "LGB") )

# 4) Plot heatmap
p_recol_covercrop <- ggplot(treatment_pw, aes(x = Group2, y = Group1, fill = `P-value`)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 2) +
  scale_fill_gradient(
    low = "darkgreen",
    high = "white",
    na.value = "grey90",
    limits = c(0, 0.1),
    name = "P-value"
  ) +
  theme_minimal() +
  #labs(x = "", y = "", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_recol_covercrop





