# Note: this code does not work yet, I ran ANCOM in QIIME instead

# This Rscript creates heatmaps

# Input: phyloseq object: dada2-table.qza, SILVAtree.qza, taxonomy.qza, metadata_phyloseq.tsv

# Output: Heatmap plots


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
library(ggrepel)
library(MicEco)

# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#filter out negative controls
ps.no = subset_samples(ps, Area %in% c("ML", "SF"))
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.no)

#-------------------------Heatmap of counts (not differential abundance)---------------------------------
# Agglomerate tse by phylum
tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_phylum <- transformCounts(tse_phylum, method = "clr", assay_name = "counts", psuedocount = 1)

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, MARGIN = "features", assay.type = "clr", 
                              method = "z", name = "clr_z")

# Take subset: only samples from mouth
tse_phylum_subset <- tse_phylum[ , tse_phylum$type %in% c("M") ]

# Add clr-transformation
tse_phylum_subset <- transformCounts(tse_phylum_subset, method = "clr",
                                     MARGIN= "samples",
                                     assay_name = "counts", pseudocount=1)
# Does z-transformation
#tse_phylum_subset <- transformCounts(tse_phylum_subset, assay_name = "clr",
#                                     MARGIN = "features", 
#                                    method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_phylum_subset)
tse_phylum <- tse_phylum_subset[top_taxa, ]



# Gets the assay table
mat <- assay(tse_phylum)

#separate colomns by Area
SF = as.data.frame(mat[,c(1:33,44)])
ML = as.data.frame(mat[, c(34:43, 45:50 )])

#make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

#make a dataframe out of the means for each area and add rownames
heatmap = data.frame(SF$Colony, ML$Colony)
rownames(heatmap) = c("Proteobacteria", "Bacteroidota", "Firmicutes", "Actinobacteriota", "Fusobacteriota")

mat <- assay(heatmap_F)

# Creates the heatmap
pheatmap(mat, main = "Heatmap of Five Most Abundant Phylums in Mouth Samples from San Francisco and Mono Lake")

#-----------------------------Heatmap: foot-------------------------------------------

# Agglomerate tse by phylum
tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_phylum <- transformCounts(tse_phylum, method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, MARGIN = "features", assay.type = "clr", 
                              method = "z", name = "clr_z")

# Take subset: only samples from foot
tse_phylum_subset_f <- tse_phylum[ , tse_phylum$type %in% c("F") ]

# Add clr-transformation
tse_phylum_subset_f <- transformCounts(tse_phylum_subset_f, method = "clr",
                                       MARGIN="samples",
                                       assay_name = "counts", pseudocount=1)
# Does z-transformation
tse_phylum_subset_f <- transformCounts(tse_phylum_subset_f, assay_name = "clr",
                                       MARGIN = "features", 
                                       method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_phylum_subset_f)
tse_phylum <- tse_phylum_subset_f[top_taxa, ]



# Gets the assay table
mat_f <- assay(tse_phylum)

#separate colomns by Area
SF = as.data.frame(mat_f[,c(1:33,44)])
ML = as.data.frame(mat_f[, c(34:43, 45:50 )])

#make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

heatmap_F = data.frame(SF$Colony, ML$Colony)
rownames(heatmap_F) = c("Firmicutes", "Proteobacteria", "Bacteroidota", "Actinobacteriota", "Halobacterota")

mat <- assay(heatmap_F)

# Creates the heatmap
pheatmap(mat, main = "Heatmap of Five Most Abundant Phylums in Foot Samples from San Francisco and Mono Lake")

#------------------------------------Heatmap: Cloaca--------------------------------------

# Agglomerate tse by phylum
tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_phylum <- transformCounts(tse_phylum, method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, MARGIN = "features", assay.type = "clr", 
                              method = "z", name = "clr_z")

# Take subset: only samples from cloaca
tse_phylum_subset_c <- tse_phylum[ , tse_phylum$type %in% c("C") ]

# Add clr-transformation
tse_phylum_subset_c <- transformCounts(tse_phylum_subset_c, method = "clr",
                                       MARGIN="samples",
                                       assay_name = "counts", pseudocount=1)
# Does z-transformation
tse_phylum_subset_c <- transformCounts(tse_phylum_subset_c, assay_name = "clr",
                                       MARGIN = "features", 
                                       method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_phylum_subset_c)
tse_phylum <- tse_phylum_subset_c[top_taxa, ]



# Gets the assay table
mat_c <- assay(tse_phylum)

#separate colomns by Area
SF = as.data.frame(mat_c[,c(1:33,44)])
ML = as.data.frame(mat_c[, c(34:43, 45:50 )])

#make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

heatmap_c = data.frame(SF$Colony, ML$Colony)
rownames(heatmap_c) = c("Firmicutes", "Actinobacteriota", "Bacteroidota", "Proteobacteria", "Fusobacteriota")

mat <- assay(heatmap_c)

# Creates the heatmap
pheatmap(mat, main = "Heatmap of Five Most Abundant Phylums in Cloaca Samples from San Francisco and Mono Lake")


#-------------------------------Family Heatmap: mouth-------------------------------------

# Agglomerate tse by phylum
tse_family <- agglomerateByRank(tse,
                                rank = "Family",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_family <- transformCounts(tse_family, method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_family <- transformCounts(tse_family, MARGIN = "features", assay.type = "clr", 
                              method = "z", name = "clr_z")

# Take subset: only samples from mouth
tse_family_subset <- tse_family[ , tse_family$type %in% c("M") ]

# Add clr-transformation
tse_family_subset <- transformCounts(tse_family_subset, method = "clr",
                                     MARGIN="samples",
                                     assay_name = "counts", pseudocount=1)
# Does z-transformation
tse_family_subset <- transformCounts(tse_family_subset, assay_name = "clr",
                                     MARGIN = "features", 
                                     method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_family_subset, top = 10)
tse_family <- tse_family_subset[top_taxa, ]



# Gets the assay table
mat <- assay(tse_family)

#separate colomns by Area
SF = as.data.frame(mat[,c(1:33,44)])
ML = as.data.frame(mat[, c(34:43, 45:50 )])

#make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

#make a dataframe out of the means for each area and add rownames
heatmap = data.frame(SF$Colony, ML$Colony)
rownames(heatmap) = c("Pasteurellaceae", "Cardiobacteriaceae", "Leptotrichiaceae
", "Micrococcaceae", "Veillonellaceae", "Weeksellaceae", "Tannerellaceae", "Moraxellaceae", "Corynebacteriaceae
", "Vibrionaceae", "Flavobacteriaceae", "Streptococcaceae", "Actinomycetaceae", "Rikenellaceae", "Porphyromonadaceae", "Fusobacteriaceae", "Prevotellaceae", "Peptostreptococcaceae", "Aerococcaceae", "Nitrincolaceae
")

mat <- assay(heatmap)

# Creates the heatmap
pheatmap(mat, main = "Heatmap of Twenty Most Abundant Families in Mouth Samples from San Francisco and Mono Lake")


#-------------------------------Family Heatmap: foot-------------------------------------

# Agglomerate tse by phylum
tse_family <- agglomerateByRank(tse,
                                rank = "Family",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_family <- transformCounts(tse_family, method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_family <- transformCounts(tse_family, method = "z", MARGIN = "features", assay.type = "clr",  name = "clr_z")

# Take subset: only samples from foot
tse_family_subset <- tse_family[ , tse_family$type %in% c("F") ]

# Add clr-transformation
tse_family_subset <- transformCounts(tse_family_subset, method = "clr",
                                     MARGIN="samples",
                                     assay_name = "counts", pseudocount=1)
# Does z-transformation
tse_family_subset <- transformCounts(tse_family_subset, assay_name = "clr",
                                     MARGIN = "features", 
                                     method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_family_subset, top = 7)
tse_family <- tse_family_subset[top_taxa, ]



# Gets the assay table
mat <- assay(tse_family)

#separate colomns by Area
SF = as.data.frame(mat[,c(1:33,44)])
ML = as.data.frame(mat[, c(34:43, 45:50 )])

# make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

# make a dataframe out of the means for each area and add rownames
heatmap = data.frame(SF$Colony, ML$Colony)
rownames(heatmap) = c("Planococcaceae", "Halomonadaceae", "Flavobacteriaceae", "Bacillaceae", "Catellicoccaceae", "Haloferacaceae
", "Rhodobacteraceae", "Fusobacteriaceae
", "Micrococcaceae
", "Nitrococcaceae
", "Nitrincolaceae
", "Streptococcaceae", "Lactobacillaceae", "Halomicrobiaceae", "Salisediminibacteriaceae", "Moraxellaceae", "ML635J-40_aquatic_group", "Clostridiaceae", "Vibrionaceae", "Leuconostocaceae
")

mat <- assay(heatmap)

# Creates the heatmap
pheatmap(mat, main = "Heatmap of Twenty Most Abundant Families in Foot Samples from San Francisco and Mono Lake")

#-----------------------------Heatmap: Cloaca----------------------------------------------

# Agglomerate tse by phylum
tse_family <- agglomerateByRank(tse,
                                rank = "Family",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_family <- transformCounts(tse_family, method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_family <- transformCounts(tse_family, MARGIN = "features", assay.type = "clr", 
                              method = "z", name = "clr_z")

# Take subset: only samples from cloaca
tse_family_subset <- tse_family[ , tse_family$type %in% c("C") ]

# Add clr-transformation
tse_family_subset <- transformCounts(tse_family_subset, method = "clr",
                                     MARGIN="samples",
                                     assay_name = "counts", pseudocount=1)
# Does z-transformation
tse_family_subset <- transformCounts(tse_family_subset, assay_name = "clr",
                                     MARGIN = "features", 
                                     method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_family_subset, top = 20)
tse_family <- tse_family_subset[top_taxa, ]



# Gets the assay table
mat <- assay(tse_family)

#separate colomns by Area
SF = as.data.frame(mat[,c(1:33,44)])
ML = as.data.frame(mat[, c(34:43, 45:50 )])

#make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

#make a dataframe out of the means for each area and add rownames
heatmap = data.frame(SF$Colony, ML$Colony)
rownames(heatmap) = c("Corynebacteriaceae", "Actinomycetaceae
", "Catellicoccaceae", "Peptostreptococcales-Tissierellales", "Hungateiclostridiaceae", "Fusobacteriaceae
", "Porphyromonadaceae
", "Lachnospiraceae
", "Veillonellaceae
", "Peptococcaceae
", "Planococcaceae
", "Enterobacteriaceae", "Synergistaceae", "Tannerellaceae", "Halomonadaceae
", "Enterococcaceae", "Lactobacillaceae", "Streptococcaceae", "uncultured
", "Micrococcaceae

")

mat <- assay(heatmap)

# Creates the heatmap
pheatmap(mat, main = "Heatmap of Twenty Most Abundant Families in Cloaca Samples from San Francisco and Mono Lake")

