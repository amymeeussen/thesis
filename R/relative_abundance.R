# This Rscript uses the package 'mia' to find the dominant bacteria in each of the different
# sample types. 

# Input: phyloseq object

# Output: Barplots of dominant phyla in three different sample types/table with
#         relative frequency of top 6 Families/Phyla


library(phyloseq)
library(tidyverse)
library(DT)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(pheatmap)


# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")


# filter out negative controls
ps.no = subset_samples(ps, Area %in% c("ML", "SF"))
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.no)

# Define ordering of labels in plots
custom_order_type <- c("M", "C", "F")
custom_order_area <- c("SF", "ML")

# Cloaca, Mouth and Foot separate - other code path
# Code from: https://www.yanh.org/2021/01/01/microbiome-r/#abundance-bar-plot
ps.rel = transform_sample_counts(ps.no, function(x) x/sum(x)*100)
glom = tax_glom(ps.rel, taxrank = "Phylum", NArm = FALSE)
ps.melt = psmelt(glom)
ps.melt$Phylum = as.character(ps.melt$Phylum)
ps.melt = ps.melt %>% group_by(type, Phylum) %>% mutate(median=median(Abundance))
keep = unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[! (ps.melt$Phylum %in% keep)] <- "< 1%"

# Make sure sum adds up to 100 for each category
ps.melt_sum = ps.melt %>% group_by(Area, type, Phylum) %>% summarise(Abundance=sum(Abundance))
sum_by_category <- ps.melt_sum %>%
  group_by(type, Area) %>%
  summarise(sum_abundance = sum(Abundance))
ps.melt_sum = merge(ps.melt_sum, sum_by_category, by=c("Area", "type"))
ps.melt_sum$Abundance = ps.melt_sum$Abundance * 100 / ps.melt_sum$sum_abundance


# Order the phylums by size
sum_phylum <- ps.melt_sum %>%
  group_by(Phylum) %>%
  summarise(sum_abundance = sum(Abundance))
sum_phylum$sum_abundance[sum_phylum$Phylum == "< 1%"] = 0
sum_phylum = sum_phylum[order(sum_phylum$sum_abundance, decreasing=TRUE), ]
custom_order_phylum = as.list(sum_phylum$Phylum)

# Set custom order for plot
ps.melt_sum$type <- factor(ps.melt_sum$type, levels = custom_order_type)
ps.melt_sum$Area <- factor(ps.melt_sum$Area, levels = custom_order_area)
ps.melt_sum$Phylum <- factor(ps.melt_sum$Phylum, levels = custom_order_phylum)

ggplot(ps.melt_sum, aes(x=Area, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity", aes(fill=Phylum)) +
  labs(x="Region", y="Relative Abundance") +
  facet_wrap(~type, scales="free_x", nrow=1, labeller = labeller(type = c("M" = "Mouth", "C" = "Cloaca","F" = "Feet"))) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  theme(text=element_text(size=15))

