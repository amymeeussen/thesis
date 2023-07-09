library(phyloseq)
library(tidyverse)
library(DT)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ANCOMBC)
library(mia)
library(pheatmap)

#This Rscript uses the package 'mia' to find the dominant bacteria in each of the different
#sample types. 
#Input: phyloseq object
#Output: Barplots of dominant phyla in three different sample types/table with raw count and
#relative frequency of top 6 Families/Phyla

# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#filter out negative controls
ps.no = subset_samples(ps, Area %in% c("ML", "SF"))

#----------------------------Explore Data-------------------------------
#Get raw count and the relative frequency of 6 dominant Families in ML and SF

# Mouth: dominant bacteria


tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.no)

tse_subset <- tse[ , tse$type %in% c("M")]
tse_ML_M = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_M = countDominantTaxa(tse_ML_M,rank = "Family")


tse_subset <- tse[ , tse$type %in% c("M")]
tse_SF_M = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_M = countDominantTaxa(tse_SF_M,rank = "Family")
dominant_taxa_SF_M

# Cloaca: dominant bacteria

tse_subset <- tse[ , tse$type %in% c("C")]
tse_ML_C = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_C = countDominantTaxa(tse_ML_C,rank = "Family")

tse_SF_C = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_C = countDominantTaxa(tse_SF_C,rank = "Family")


# Foot: dominant bacteria

tse_subset <- tse[ , tse$type %in% c("F")]
tse_ML_F = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_F = countDominantTaxa(tse_ML_F,rank = "Family")

tse_SF_F = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_F = countDominantTaxa(tse_SF_F,rank = "Family")

# Stacked barplot mouth

dominant_mouth = rbind(dominant_taxa_ML_M, dominant_taxa_SF_M)

dominant_mouth$Area <- c(rep("ML", 6), rep("SF", 6))

ggplot(dominant_mouth, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Families in Mono Lake and SF Bay Mouth Samples")

# Stacked barplot cloaca

dominant_foot = rbind(dominant_taxa_ML_F, dominant_taxa_SF_F)

dominant_foot$Area <- c(rep("ML", 10), rep("SF", 18))

ggplot(dominant_foot, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Families in Mono Lake and SF Bay Foot Samples")

# Stacked barplot cloaca

dominant_cloaca = rbind(dominant_taxa_ML_C, dominant_taxa_SF_C)

dominant_cloaca$Area <- c(rep("ML", 10), rep("SF", 12))

ggplot(dominant_cloaca, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Families in Mono Lake and SF Bay Cloaca Samples")


#Mouth: Barplot level phylum

tse_subset <- tse[ , tse$type %in% c("M")]
tse_ML_M = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_M = countDominantTaxa(tse_ML_M,rank = "Phylum")

tse_subset <- tse[ , tse$type %in% c("M")]
tse_SF_M = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_M = countDominantTaxa(tse_SF_M,rank = "Phylum")

dominant_mouth = rbind(dominant_taxa_ML_M, dominant_taxa_SF_M)

dominant_mouth$Area <- c(rep("ML", 3), rep("SF", 3))

ggplot(dominant_mouth, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Phylum in Mono Lake and SF Bay Mouth Samples")

# Cloaca: rel.abundance (phylum)

tse_subset <- tse[ , tse$type %in% c("C")]
tse_ML_C = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_C = countDominantTaxa(tse_ML_C,rank = "Phylum")

tse_subset <- tse[ , tse$type %in% c("C")]
tse_SF_C = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_C = countDominantTaxa(tse_SF_C,rank = "Phylum")

dominant_cloaca = rbind(dominant_taxa_ML_C, dominant_taxa_SF_C)

dominant_mouth$Area <- c(rep("ML", 5), rep("SF", 5))

ggplot(dominant_mouth, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Phylum in Mono Lake and SF Bay Cloacal Samples")


# Foot: rel.abundance (phylum)

tse_subset <- tse[ , tse$type %in% c("F")]
tse_ML_F = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_F = countDominantTaxa(tse_ML_F,rank = "Phylum")

tse_subset <- tse[ , tse$type %in% c("F")]
tse_SF_F = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_F = countDominantTaxa(tse_SF_F,rank = "Phylum")

dominant_foot = rbind(dominant_taxa_ML_F, dominant_taxa_SF_F)

dominant_foot$Area <- c(rep("ML", 5), rep("SF", 4))

ggplot(dominant_foot, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Phylum in Mono Lake and SF Bay Foot Samples")

