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


# This Rscript uses the phyloseq packages to rarefy the samples, it creates shannon, observed richness, and 
#Faith's PD diversity plots + unweighted, weighted unifrac, jaccard, and bray ordination plots along
# with a PERMANOVA tests. 
# Input: phyloseq object
# Output: Rarification plot, venn diagrams of shared/unique taxa, alpha diversity plots by sample type, 
# and ordination plots + PERMANOVA results


# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

# Remove negative controls

ps.no = subset_samples(ps, Area %in% c("ML", "SF"))



#-------------------------------Rarefy-----------------------------------

#rarefy samples without controls

set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps.no, rngseed=1, sample.size=11000, replace=F)

print(ps.rarefied)

#Check that sample sizes are all even
barplot(sample_sums(ps.rarefied), las =2)

#---------------------------Shared and unique taxa between groups--------

#unique and shared ASVs per area
ps_venn(ps.rarefied, fraction = 0.01, group = "Area")

#unique and shared ASVs in mouth samples
ps.mouth = subset_samples(p.rarefied, type %in% "M")
ps_venn(ps.mouth, fraction = 0.10, group = "Area") 
 
#unique and shared ASVs in foot samples
ps.foot = subset_samples(ps.rarefied, type %in% "F")
ps_venn(ps.foot, fraction = 0.10, group = "Area")

#unique and shared ASVs in cloaca samples
ps.cloaca = subset_samples(ps.rarefied, type %in% "C")
ps_venn(ps.cloaca, fraction = 0.10, group = "Area")



#---------------------------------ALL SAMPLES ALPHA and BETA DIVERSITY------------------------


#Shannon alpha diversity for two areas
p = plot_richness(ps.rarefied, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_q2r()+
  theme(text = element_text(size = 15))+
#  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")+
  facet_wrap(~type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth")))

p + stat_compare_means(method = "wilcox.test", label.x = 1.1, label.y = 6)

p = plot_richness(ps.rarefied, x="Area", measures="ACE", color = "Area") +
  geom_boxplot() +
  theme_q2r()+
  theme(text = element_text(size = 15))+
  ggtitle("ACE Diversity in CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("ACE index")+
  facet_wrap(~type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth")))

p + stat_compare_means(method = "wilcox.test", label.x = 1.1, label.y = 6)

p = plot_richness(ps.rarefied, x="Area", measures="Observed", color = "Area") +
  geom_boxplot() +
  theme_q2r()+
  theme(text = element_text(size = 15))+
  ggtitle("Observed Diversity in CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Observed")+
  facet_wrap(~type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth")))

p + stat_compare_means(method = "wilcox.test", label.x = 1.1, label.y = 6)


----------------Faith's Phylogenetic Diversity across areas, by type------------------

adiv <- data.frame(
    "Observed" = estimate_richness(ps.rarefied, measures = "Observed"),
    "Shannon" = phyloseq::estimate_richness(ps.rarefied, measures = "Shannon"),
    "PD" = picante::pd(samp = data.frame(t(data.frame(otu_table(ps.rarefied)))), tree = phyloseq::phy_tree(ps.rarefied))[, 1],
    "Pielou" = microbiome::evenness(ps.rarefied, index = "pielou", zeroes = TRUE, detection = 0),
    "Area" = sample_data(ps.rarefied)$Area, 
    "type" = sample_data(ps.rarefied)$type)
head(adiv)

#Faith's PD, by type

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
  

#Shannon's by type

adiv %>%
  gather(key = metric, value = value, c("Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Shannon"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "Location", y = "Shannon/s diversity index") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme(legend.position="none") +
  theme_bw()
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Shannon Diversity Index for each body site")
 

#Observed richness, by type

adiv %>%
  gather(key = metric, value = value, c("Observed")) %>%
  mutate(metric = factor(metric, levels = c("Observed"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "Location", y = "Number of observed taxa") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Observed richness of CAGU samples in SF and ML, by type") 

#Pielou's evenness across areas, by type

adiv %>%
  gather(key = metric, value = value, c("pielou")) %>%
  mutate(metric = factor(metric, levels = c("pielou"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "Location", y = "Pielou/'s evenness index") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Pielou/'s evenness index for each body site")
  


#----------------------------------Mouth Samples Alpha and BETA DIVERSITY---------------------

#filter out mouth samples
ps.rarefied.mouth = subset_samples(ps.rarefied, type %in% c("M"))

#Shannon alpha diversity
q = plot_richness(ps.rarefied.mouth, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Mouth Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity") +
q + stat_compare_means(method = "wilcoxon.test", label.x = 1.5, label.y = 6)

#unweighted unifrac PCoA
unifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: Mouth samples")

#PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist ~ sample_data(ps.rarefied.mouth)$Area)


#weighted unifrac PCoA
wunifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Weighted Unifrac: Mouth samples")

#PERMANOVA test for weighted unifrac
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.mouth)$Area)

# Bray curtis PCoA
bray_curtis = phyloseq::distance(ps.rarefied.mouth, method = "bray")
ordination = ordinate(ps.rarefied.mouth, method="PCoA")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: Mouth samples")

#PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis ~ sample_data(ps.rarefied.mouth)$Area)

#Jaccard PCoA
Jaccard = phyloseq::distance(ps.rarefied.mouth, method = "jaccard")
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Mouth samples")

#PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard ~ sample_data(ps.rarefied.mouth)$Area)


#----------------------------------Cloaca Samples Alpha and BETA DIVERSITY---------------------

#filter out cloaca samples
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.cloaca, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 6)

#unweighted unifrac PCoA
unifrac_dist.c = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: Cloaca samples")

#PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist.c ~ sample_data(ps.rarefied.cloaca)$Area)

#weighted unifrac test
wunifrac_dist.c = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) +
  theme_classic() + 
  ggtitle("Weighted Unifrac PCoA: Cloaca")

#PERMANOVA test
vegan::adonis2(wunifrac_dist.c ~ sample_data(ps.rarefied.cloaca)$Area)

# Bray curtis PCoA
bray_curtis.c = phyloseq::distance(ps.rarefied.mouth, method = "bray")
ordination = ordinate(ps.rarefied.mouth, method="PCoA")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: Cloaca samples")

#PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis.c ~ sample_data(ps.rarefied.mouth)$Area)

#Jaccard PCoA
Jaccard.c = phyloseq::distance(ps.rarefied.mouth, method = "jaccard")
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Mouth samples")

#PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard.c ~ sample_data(ps.rarefied.mouth)$Area)


#----------------------------------Foot Samples Alpha and BETA DIVERSITY---------------------

#filter out foot samples
ps.rarefied.foot = subset_samples(ps.rarefied.new, type %in% c("F"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.foot, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Foot Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#unweighted unifrac PCoA
unifrac_dist.f = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance="unifrac")
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Unifrac: Cloaca samples")

#PERMANOVA test for unweighted unifrac
vegan::adonis2(unifrac_dist.f ~ sample_data(ps.rarefied.cloaca)$Area)


#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="Area") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.foot)$Area)

# Bray curtis PCoA
bray_curtis.f = phyloseq::distance(ps.rarefied.mouth, method = "bray")
ordination = ordinate(ps.rarefied.mouth, method="PCoA")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Bray-curtis: Cloaca samples")

#PERMANOVA test for Bray curtis 
vegan::adonis2(bray_curtis.f ~ sample_data(ps.rarefied.mouth)$Area)

#Jaccard PCoA
Jaccard.c = phyloseq::distance(ps.rarefied.mouth, method = "jaccard")
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance = "jaccard")
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) + 
  theme_bw() + 
  ggtitle("Jaccard: Mouth samples")

#PERMANOVA test for Jaccard 
vegan::adonis2(Jaccard.c ~ sample_data(ps.rarefied.mouth)$Area)


#----------------------------------Cloaca Samples (M vs. F) Alpha and BETA DIVERSITY---------------------

#filter out cloaca samples
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.cloaca, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="sex") + theme(aspect.ratio=1)

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$sex)

#----------------------------------Mouth Samples (M vs. F) Alpha and BETA DIVERSITY---------------------

#filter out mouth samples
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("M"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.cloaca, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="sex") + theme(aspect.ratio=1)

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$sex)

#----------------------------------Foot Samples (M vs. F) Alpha and BETA DIVERSITY---------------------

#filter out mouth samples
ps.rarefied.foot = subset_samples(ps.rarefied, type %in% c("F"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.foot, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Cloaca Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="sex") + theme(aspect.ratio=1)

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.foot)$sex)

#----------------------------------Body Condition by area and sample type Alpha and BETA DIVERSITY---------------------

#filter out foot
ps.rarefied.foot = subset_samples(ps.rarefied, type %in% c("F"))
ps.rarefied.foot.ML = subset_samples(ps.rarefied.foot, Area %in% "ML")

#Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.foot.ML, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.foot.SF, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="condition_score") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.foot.SF)$condition_score)

#filter out mouth
ps.rarefied.mouth = subset_samples(ps.rarefied, type %in% c("M"))

#Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.mouth, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.mouth, ordination, color="condition_score") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.mouth)$condition_score)

#filter out cloaca
ps.rarefied.cloaca = subset_samples(ps.rarefied, type %in% c("C"))

#Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.cloaca, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="condition_score") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$condition_score)


#------------------------Sex-----------------------------
# Do not use until sex data comes in. 

#Shannon alpha diversity for sexes (sexes )
p = plot_richness(ps.rarefied, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 2)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.new, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.new, ordination, color="sex") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied)$sex)



