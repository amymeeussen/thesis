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


#Relative Abundance of Phylums between sample types
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)



#Relative Abundance of Phylums between Areas
ps.melt <- ps.melt %>%
  group_by(type, Phylum) %>%
  mutate(median=median(Abundance))
# select group median > 1%
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Area,type,Phylum) %>%
  summarise(Abundance=sum(Abundance))

ps.area_sum <- ps.melt %>%
  group_by(Area) %>%
  summarise(Abundance=sum(Abundance))



ggplot(ps.melt_sum, aes(x = Area, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="%") +
#  facet_wrap(~sex, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))








#-------------------------------Rarefy-----------------------------------

#rarefy samples without controls

set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps.no, rngseed=1, sample.size=15000, replace=F)

print(ps.rarefied)

#Check that sample sizes are all even
barplot(sample_sums(ps.rarefied), las =2)

#rarefy samples for body condition (without samples 19 and 41)
set.seed(111) # keep result reproductive
ps.rarefied.bc = rarefy_even_depth(ps.bc, rngseed=1, sample.size=15000, replace=F)

#---------------------------Shared and unique taxa between groups--------

#unique and shared ASVs per area
ps_venn(ps.no, fraction = 0.01, group = "Area")

#unique and shared ASVs in mouth samples
ps.mouth = subset_samples(ps.no, type %in% "M")
ps_venn(ps.mouth, fraction = 0.10, group = "Area") 
 
#unique and shared ASVs in foot samples
ps.foot = subset_samples(ps.no, type %in% "F")
ps_venn(ps.foot, fraction = 0.10, group = "Area")

#unique and shared ASVs in cloaca samples
ps.cloaca = subset_samples(ps.no, type %in% "C")
ps_venn(ps.cloaca, fraction = 0.10, group = "Area")



#---------------------------------ALL SAMPLES ALPHA and BETA DIVERSITY------------------------

#Test for Shannon diversity with phyloseq package
sd = estimate_richness(ps.no, measures = "Shannon")
sd

#Test to see if shannon diversity has an equal distribution
shapiro.test(sd$Shannon)


#create a dataframe with diversity metrics and groups you want to test for diversity. You will need this
#test variance of diversity metrics between groups. 

diversity = as.data.frame(sd)
meta = meta(ps.no)
diversity$Area = meta$Area
diversity$bc = meta$condition_score
diversity$sex = meta$sex
diversity$condition = meta$body_condition

ML = filter(diversity, Area == "ML")
SF = filter(diversity, Area == "SF")

hist(ML$condition)
hist(SF$condition)




#test for normality
shapiro.test(diversity$condition)
#use levene's test to test for equal variance, since data is not normal
result = leveneTest(Shannon ~ interaction(Area), data = diversity)

#Test for equal variance (2nd assumption of the wilcoxon rank sum test)
var.test(diversity$Shannon ~ diversity$Area, alternative = "two.sided")

#Remove birds that don't have bc or sex data
diversity$Bird = meta$Bird
diversity_zbc = diversity %>% filter(!Bird %in% c("41", "37", "2", "19", "57"))

#test for equal variance between body conditions and sexes
var.test(diversity_zbc$Shannon ~ diversity_zbc$bc, alternative = "two.sided")
var.test(diversity_zbc$Shannon ~ diversity_zbc$sex, alternative = "two.sided")



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
  labs(x = "", y = "") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme_pubr() +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Faith\'s Phylogenetic Diversity Index of CAGU samples in SF and ML, by type") +
  ylab("Faith\'s PD index") +
  xlab("Area") 

#Shannon's by type

adiv %>%
  gather(key = metric, value = value, c("Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Shannon"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Shannon Diversity Index of CAGU samples in SF and ML, by type") +
  ylab("Shannon diversity index") +
  xlab("Area")

#Observed richness, by type

adiv %>%
  gather(key = metric, value = value, c("Observed")) %>%
  mutate(metric = factor(metric, levels = c("Observed"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Observed richness of CAGU samples in SF and ML, by type") +
  ylab("Number of taxa observed") +
  xlab("Area")

#Pielou's evenness across areas, by type

adiv %>%
  gather(key = metric, value = value, c("pielou")) %>%
  mutate(metric = factor(metric, levels = c("pielou"))) %>%
  ggplot(aes(x = Area, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Area), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ type, labeller = labeller(type = c("C" = "Cloaca","F" = "Foot", "M" = "Mouth"))) +
  theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1, label.y = 1) + 
  ggtitle("Pielou/'s evenness index for CAGU samples in SF and ML, by type") +
  ylab("Pielou/'s evenness index") +
  xlab("Area")





------------------------Body condition-----------------------------

#Retrieve metadata from phyloseq object
meta = meta(ps)

#Count how many samples with low and high body condition score in each colony
ML_low = filter(meta, Area == "ML", condition_score == "low")
nrow(ML_low)
#20
ML_high = filter(meta, Area == "ML", condition_score == "high")
nrow(ML_high)
#2
SF_low = filter(meta, Area == "SF", condition_score == "low")
nrow(SF_low)
#15
SF_high = filter(meta, Area == "SF", condition_score == "high")
nrow(SF_high)
#22
  
  
#Observed, ACE, and Shannon alpha diversity for two body conditions
p = plot_richness(ps.rarefied.bc, x="condition_score", measures= c("Shannon","Observed", "ACE"), color = "condition_score") +
  geom_boxplot() +
  theme_bw() +
 # facet_wrap(~Area) +
  ggtitle("Diversity in CAGU Samples Taken in Mono Lake and San Francisco Bay: High and Low Body Conditions") +
  ylab("Diversity")+
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 6)
p

#weighted unifrac test
wunifrac_dist = phyloseq::distance(ps.rarefied.bc, method="unifrac", weighted=T)
ordination = ordinate(ps.rarefied.new, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.bc, ordination, color="condition_score") + theme(aspect.ratio=1) 

#PERMANOVA test
vegan::adonis2(wunifrac_dist ~ sample_data(ps.rarefied.bc)$condition_score)



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



#----------------------------------Mouth Samples Alpha and BETA DIVERSITY---------------------

#filter out mouth samples
ps.rarefied.mouth = subset_samples(ps.rarefied, type %in% c("M"))

#Shannon alpha diversity
p = plot_richness(ps.rarefied.mouth, x="Area", measures="Shannon", color = "Area") +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  ggtitle("Shannon Diversity of CAGU Mouth Samples Taken in Mono Lake and San Francisco Bay") +
  ylab("Shannon Diversity") +
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


  


