# This Rscript uses the phyloseq packages to rarefy the samples, it creates
# shannon, observed richness, and Faith's PD diversity plots + unweighted,
# weighted unifrac, jaccard, and bray ordination plots along with a PERMANOVA tests.

# Input: phyloseq object: dada2-table.qza, SILVAtree.qza, taxonomy.qza, metadata_phyloseq.tsv
# Output: Rarification plot, venn diagrams of shared/unique taxa, alpha diversity plots by sample type, 
#         and ordination plots + PERMANOVA results

library(phyloseq)
library(tidyverse)
library(DT)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ANCOMBC)
library(mia)
library(ade4)
library(dplyr)
library(MicEco)
library(breakaway)
library(microbiome)
library(ggplot2)
library(vegan)
library(ggsci)

blue = rgb(60/255,191/255,196/255)
red = rgb(247/255,123/255,114/255)
green = rgb(100/255,180/255,136/255)
custom_order_type <- c("M", "C", "F")
custom_order_area <- c("SF", "ML")

# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

# Remove negative controls

ps.no = subset_samples(ps, Area %in% c("ML", "SF"))



#-------------------------------Rarefy-----------------------------------

# rarefy samples without controls

set.seed(999) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps.no, rngseed=1, replace=F)


# Check that sample sizes are all even
barplot(sample_sums(ps.rarefied), las =2)

print(barplot)

#---------------------------Shared and unique taxa between groups--------

# unique and shared ASVs per area
ps_venn(ps.rarefied, fraction = 0.01, group = "Area")

# unique and shared ASVs in mouth samples
ps.mouth = subset_samples(ps.rarefied, type %in% "M")
ps_venn(ps.mouth, fraction = 0.10, group = "Area") 
 
# unique and shared ASVs in foot samples
ps.foot = subset_samples(ps.rarefied, type %in% "F")
ps_venn(ps.foot, fraction = 0.10, group = "Area")

# unique and shared ASVs in cloaca samples
ps.cloaca = subset_samples(ps.rarefied, type %in% "C")
ps_venn(ps.cloaca, fraction = 0.10, group = "Area")



#---------------------------------ALL SAMPLES ALPHA DIVERSITY------------------------

# Create dataframe with Observed, Shannon's, Faith's PD, and Pielou's evenness diversity metrics

adiv = data.frame(
    "Observed" = estimate_richness(ps.rarefied, measures = "Observed"),
    "Shannon" = phyloseq::estimate_richness(ps.rarefied, measures = "Shannon"),
    "PD" = picante::pd(samp = data.frame(t(data.frame(otu_table(ps.rarefied)))), tree = phyloseq::phy_tree(ps.rarefied))[, 1],
    "Pielou" = microbiome::evenness(ps.rarefied, index = "pielou", zeroes = TRUE, detection = 0),
    "Area" = sample_data(ps.rarefied)$Area, 
    "type" = sample_data(ps.rarefied)$type, 
    "condition" = sample_data(ps.rarefied)$body_condition)
head(adiv)
adiv$type = factor(adiv$type, levels = custom_order_type)
adiv$Area = factor(adiv$Area, levels = custom_order_area)

# ------------------------------------Diversity by Area-------------------------------

# Shannon Diversity, by area

adiv %>%
  gather(key = metric, value = value, c("Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Shannon"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "Population", y = "Shannon diversity index") +
  theme(legend.position="none") +
  theme_q2r() + 
  labs(color = "Population") +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 0))+
  theme(text = element_text(size = 20)) +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) 

# Faith's PD, by area

adiv %>%
  gather(key = metric, value = value, c("PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  adiv %>% theme(legend.position="none") +
  labs(x = "Location", y = "Faith\'s PD index") +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Faith\'s Phylogenetic Diversity Index of CAGU samples in SF and ML, by type") 

# Peilou's Evenness, by area

adiv %>%
  gather(key = metric, value = value, c("pielou")) %>%
  mutate(metric = factor(metric, levels = c("pielou"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Region), height = 0, width = .2) +
  labs(x = "Location", y = "Pielou/'s evenness index") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Feet", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Pielou/'s evenness index for each body site")


# ------------------------------------Diversity by body site------------------------------------

# Shannon's by type
adiv %>%
  gather(key = metric, value = value, c("Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Shannon"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "Population", y = "Shannon diversity index") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Feet", "M" = "Mouth"))) +
  theme(legend.position="none") +
  theme_q2r() + 
  labs(color = "Population") +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 0))+
  theme(text = element_text(size = 20)) +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) 


# Faith's PD, by type
adiv %>%
  gather(key = metric, value = value, c("PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  adiv %>% theme(legend.position="none") +
  labs(x = "Location", y = "Faith\'s PD index") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Faith\'s Phylogenetic Diversity Index of CAGU samples in SF and ML, by type") 
  


# Observed richness, by type
adiv %>%
  gather(key = metric, value = value, c("Observed")) %>%
  mutate(metric = factor(metric, levels = c("Observed"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "Location", y = "Number of observed taxa") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Feet", "M" = "Mouth"))) +
  theme(legend.position="none") +
  theme_classic() + 
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))+
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Observed richness of CAGU samples in SF and ML, by type") 



# Pielou's evenness across areas, by type
adiv$Region = adiv$Area
adiv %>%
  gather(key = metric, value = value, c("pielou")) %>%
  mutate(metric = factor(metric, levels = c("pielou"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Region), height = 0, width = .2) +
  labs(x = "Location", y = "Pielou/'s evenness index") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Feet", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Pielou/'s evenness index for each body site")
  


#----------------------------------Mouth Samples Alpha and BETA DIVERSITY---------------------

# filter out mouth samples
ps.rarefied.mouth = subset_samples(ps.rarefied, type %in% c("M"))
ps.rarefied.foot = subset_samples(ps.rarefied, type %in% c("F"))
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

# Shannon alpha diversity
q = plot_richness(ps.rarefied.mouth, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Mouth Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
q + stat_compare_means(method = "wilcoxon.test", label.x = 1.5, label.y = 6)

# unweighted unifrac PCoA
unifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: Mouth samples")

# PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist ~ sample_data(ps.rarefied.mouth)$Area)


# weighted unifrac PCoA
wunifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_classic() + 
  theme(text = element_text(size = 20))+
  ggtitle("Weighted Unifrac: Mouth samples")

# PERMANOVA test for weighted unifrac
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.mouth)$Area)

# Bray curtis PCoA
bray_curtis = phyloseq::distance(ps.rarefied.mouth, method = "bray")
ordination = ordinate(ps.rarefied.mouth, method="PCoA")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: Mouth samples")

# PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis ~ sample_data(ps.rarefied.mouth)$Area)

# Jaccard PCoA
Jaccard = phyloseq::distance(ps.rarefied.mouth, method = "jaccard")
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Mouth samples")

# PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard ~ sample_data(ps.rarefied.mouth)$Area)


#----------------------------------Cloaca Samples Alpha and BETA DIVERSITY---------------------

# filter out cloaca samples
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

# Shannon alpha diversity
p = plot_richness(ps.rarefied.cloaca, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 6)

# unweighted unifrac PCoA
unifrac_dist.c = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: Cloaca samples")

# PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist.c ~ sample_data(ps.rarefied.cloaca)$Area)

# weighted unifrac test (graph not working)
wunifrac_dist.c = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) +
  theme_classic() + 
  ggtitle("Weighted Unifrac PCoA: Cloaca")

# PERMANOVA test
vegan::adonis2(wunifrac_dist.c ~ sample_data(ps.rarefied.cloaca)$Area)

# Bray curtis PCoA
bray_curtis.c = phyloseq::distance(ps.rarefied.cloaca, method = "bray")
ordination = ordinate(ps.rarefied.cloaca, method="PCoA")
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: Cloaca samples")

# PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis.c ~ sample_data(ps.rarefied.cloaca)$Area)

# Jaccard PCoA
Jaccard.c = phyloseq::distance(ps.rarefied.cloaca, method = "jaccard")
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Cloaca samples")

# PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard.c ~ sample_data(ps.rarefied.cloaca)$Area)


#----------------------------------Foot Samples Alpha and BETA DIVERSITY---------------------

# filter out foot samples
ps.rarefied.foot = subset_samples(ps.rarefied.new, type %in% c("F"))

# Shannon alpha diversity
p = plot_richness(ps.rarefied.foot, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Foot Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 6)

# unweighted unifrac PCoA
unifrac_dist.f = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied.foot, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: Foot samples")

# PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist.f ~ sample_data(ps.rarefied.foot)$Area)


# weighted unifrac test
wunifrac_dist.f = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="Area") + 
  theme(aspect.ratio=1) + 
  ggtitle("Weighted Unifrac: Foot")

# PERMANOVA test
vegan::adonis2(wunifrac_dist.f ~ sample_data(ps.rarefied.foot)$Area)

# Bray curtis PCoA
bray_curtis.f = phyloseq::distance(ps.rarefied.foot, method = "bray")
ordination = ordinate(ps.rarefied.foot, method="PCoA")
plot_ordination(ps.rarefied.foot, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: Foot samples")

# PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis.f ~ sample_data(ps.rarefied.foot)$Area)

# Jaccard PCoA
Jaccard.f = phyloseq::distance(ps.rarefied.foot, method = "jaccard")
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied.foot, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Mouth samples")

# PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard.f ~ sample_data(ps.rarefied.foot)$Area)




#--------------------------------Between type Alpha and Beta diversity----------------------------

# Shannon alpha diversity
q = plot_richness(ps.rarefied, x="type", measures="Shannon", color = "type") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 180)) +
  ggtitle("Shannon Diversity of Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
q + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)
                       

# unweighted unifrac PCoA
unifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied, ordination, color="type") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: body site")

# PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist ~ sample_data(ps.rarefied)$type)

# weighted unifrac PCoA by body site
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination, color="type") + theme(aspect.ratio=1) + 
  ggtitle("Weighted Unifrac: body sites (both region included)") +
  scale_color_aaas() +
  theme(text = element_text(size = 15)) +
  theme_classic() + 
  labs(color = "Body Site") +
  stat_ellipse(type = "norm", level=0.90, linetype = 2)

# weighted unifrac PCoA for mouth by area
wunifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  ggtitle("Weighted Unifrac for Mouth samples by Area") +
  labs(color = "Population") +
  stat_ellipse(type = "norm", level=0.90, linetype = 2)

# weighted unifrac PCoA for foot by area
wunifrac_dist = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  ggtitle("Weighted Unifrac for Foot samples by Area") +
  labs(color = "Population") +
  stat_ellipse(type = "norm", level=0.9, linetype = 2)

# weighted unifrac PCoA for cloaca by area
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  ggtitle("Weighted Unifrac for Cloaca samples by Area") +
  labs(color = "Population") +
  stat_ellipse(type = "norm", level=0.9, linetype = 2)

# PERMANOVA test for weighted unifrac
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$Area)

# Bray curtis PCoA
bray_curtis = phyloseq::distance(ps.rarefied, method = "bray")
ordination = ordinate(ps.rarefied, method="PCoA")
plot_ordination(ps.rarefied, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: body site")

# PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis ~ sample_data(ps.rarefied.mouth)$type)

# Jaccard PCoA
Jaccard = phyloseq::distance(ps.rarefied, method = "jaccard")
ordination = ordinate(ps.rarefied, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Mouth samples")

# PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard ~ sample_data(ps.rarefied.mouth)$type)






#----------------------------------Cloaca Samples (M vs. F) Alpha and BETA DIVERSITY---------------------

# filter out cloaca samples
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

# Shannon alpha diversity
p = plot_richness(ps.rarefied.cloaca, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="sex") + theme(aspect.ratio=1)

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$sex)

#----------------------------------Mouth Samples (M vs. F) Alpha and BETA DIVERSITY---------------------

# filter out mouth samples
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("M"))

# Shannon alpha diversity
p = plot_richness(ps.rarefied.cloaca, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="sex") + theme(aspect.ratio=1)

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$sex)

#----------------------------------Foot Samples (M vs. F) Alpha and BETA DIVERSITY---------------------

# filter out mouth samples
ps.rarefied.foot = subset_samples(ps.rarefied, type %in% c("F"))

# Shannon alpha diversity
p = plot_richness(ps.rarefied.foot, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="sex") + theme(aspect.ratio=1)

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.foot)$sex)

#----------------------------------Body Condition by area and sample type Alpha and BETA DIVERSITY---------------------

# filter out foot
ps.rarefied.foot = subset_samples(ps.rarefied, type %in% c("F"))
ps.rarefied.F.no = subset_samples(ps.rarefied.foot, !Bird %in% c("19", "41"))

# Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.F.no, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 6)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.F.no, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="condition_score") + theme(aspect.ratio=1) 

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.F.no)$condition_score)

# filter out mouth
ps.rarefied.mouth = subset_samples(ps.rarefied, type %in% c("M"))
ps.rarefied.M.no = subset_samples(ps.rarefied.mouth, !Bird %in% c("19", "41"))

# Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.M.no, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.mouth, ordination, color="condition_score") + theme(aspect.ratio=1) 

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.mouth)$condition_score)

# filter out cloaca
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

# Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.cloaca, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="condition_score") + theme(aspect.ratio=1) 

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$condition_score)


#------------------------Sex-----------------------------
# Do not use until sex data comes in. 

# Shannon alpha diversity for sexes (sexes )
p = plot_richness(ps.rarefied, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 2)

# weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.new, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.new, ordination, color="sex") + theme(aspect.ratio=1) 

# PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied)$sex)






