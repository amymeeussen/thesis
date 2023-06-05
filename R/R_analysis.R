library(phyloseq)
library(tidyverse)
library(DT)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ANCOMBC)
library(mia)

# This Rscript uses the phyloseq packages to create relative abundance tables per sample; it rarefies
# the samples, and it creates shannon alpha diversity plots + weighted unifrac ordination plots along
# with a PERMANOVA test. 
# Input: phyloseq object
# Output: Relative abundance bar chart, rarification plot, shannon diversity plots by sample type, 
# and weighted unifrac ordination plot + PERMANOVA results


# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#---------------------------------Unrarefied Data Exploration----------------------------------


#filter out negative controls
ps.no = subset_samples(ps, Area %in% c("ML", "SF"))

#filter out empty rows for body condition
ps.bc.one = subset_samples(ps.no, Bird != 19)
ps.bc = subset_samples(ps.bc.one, Bird != 41)


#Relative Abundance of Phylums between sample types
ps.rel = transform_sample_counts(ps.no, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(type, Phylum) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Bird,type,Phylum) %>%
  summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Bird, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="%") +
  facet_wrap(~type, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))


#Relative Abundance of Phylums between Areas
ps.melt <- ps.melt %>%
  group_by(sex, Phylum) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Bird,sex,Phylum) %>%
  summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Bird, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="%") +
  facet_wrap(~sex, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))


#-------------------------------Rarefy-----------------------------------

#rarefy samples without controls

set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps.no, rngseed=1, sample.size=15000, replace=F)

#rarefy samples for body condition (without samples 19 and 41)
set.seed(111) # keep result reproductive
ps.rarefied.bc = rarefy_even_depth(ps.bc, rngseed=1, sample.size=15000, replace=F)


#---------------------------------ALL SAMPLES ALPHA and BETA DIVERSITY------------------------

#Shannon alpha diversity for two areas
p = plot_richness(ps.rarefied, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.new, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.new, ordination, color="Area") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied)$Area)

#Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.bc, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.bc, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.new, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.bc, ordination, color="condition_score") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.bc)$condition_score)

#Shannon alpha diversity for sexes
p = plot_richness(ps.rarefied, x="sex", measures="Shannon", color = "sex") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.new, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.new, ordination, color="sex") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied)$sex)









#----------------------------------Mouth Samples Alpha and BETA DIVERSITY---------------------

#filter out mouth samples
ps.rarefied.mouth = subset_samples(ps.rarefied, type %in% c("M"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.mouth, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Mouth Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.mouth, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.mouth, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.mouth, ordination, color="Area") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.mouth)$Area)

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
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.cloaca, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.cloaca, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.cloaca, ordination, color="Area") + theme(aspect.ratio=1)

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.cloaca)$Area)

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

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="Area") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.foot)$Area)


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

#----------------------------------Body Condition (by sample type) Alpha and BETA DIVERSITY---------------------

#filter out foot
ps.rarefied.foot = subset_samples(ps.rarefied, type %in% c("F"))

#Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.foot, x="condition_score", measures="Shannon", color = "condition_score") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Shannon Diversity")
p + stat_compare_means(method = "kruskal.test", label.x = 1.5, label.y = 6)

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.foot, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.foot, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied.foot, ordination, color="condition_score") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.foot)$condition_score)

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



