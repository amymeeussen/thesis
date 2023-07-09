# Note: this code does not work yet, I ran ANCOM in QIIME instead

# This Rscript runs the phyloseq and mia package find relative abundances, run ANCOMBC
# and create a volcano plot

# Input: phyloseq object: dada2-table.qza, SILVAtree.qza, taxonomy.qza, metadata_phyloseq.tsv

# Output: None yet


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

#----------------------------Explore Data-------------------------------

# Mouth: dominant bacteria

tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.no)

tse_subset <- tse[ , tse$type %in% c("M")]
tse_ML_M = tse_subset[, tse_subset$Area %in% "ML"]

dominant_taxa_ML_M = countDominantTaxa(tse_ML_M,rank = "Family")


tse_subset <- tse[ , tse$type %in% c("M")]
tse_SF_M = tse_subset[, tse_subset$Area %in% "SF"]

dominant_taxa_SF_M = countDominantTaxa(tse_SF_M,rank = "Family")


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

dominant_foot$Area <- c(rep("ML", 4), rep("SF", 4))

ggplot(dominant_foot, aes(x = Area, y = rel.freq, fill = dominant_taxa)) + 
  geom_col(aes(fill = dominant_taxa)) +
  ggtitle("Dominant Phylum in Mono Lake and SF Bay Foot Samples")
------------------------DESeq2-----------------------
  
# pairwise comparison between Areas
ps.taxa.sub <- subset_samples(ps, Area %in% c("ML", "SF"))
# filter sparse features, with > 90% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)
ps_ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ Area)
# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(ps_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
# Only keep gut and tongue samples
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.pse.sub)), ps.taxa.rel.sig)


matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
# Define the annotation color for columns and rows
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
  Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
)
rownames(annotation_row) = rownames(matrix)

# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  `Body site` = c(gut = "purple", tongue = "yellow"),
  Phylum = phylum_col
)

ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors)
  
  



#-----------------------------ANCOMBC------------------------------
out <- ancombc(phyloseq = ps.no, formula = "Area", 
               p_adj_method = "holm", lib_cut = 1000, 
               group = "Area", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
#out <- ancombc(
#  phyloseq = ps,
 # p_adj_method = "fdr",
  #formula = "group",
  #lib_cut = 0, 
  #group = "Area", 
  #struc_zero = TRUE, 
  #neg_lb = TRUE,
  #alpha = 0.05, 
  #global = TRUE # multi group comparison will be deactivated automatically 
#)
res = out$res

res.or_p <- rownames(res$q_val["AreaSF"])[base::order(res$q_val[,"AreaSF"])]
taxa_sig <- res.or_p[1:20]
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
# Only keep gut and tongue samples 
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.sub)), ps.taxa.rel.sig)





---------------
tab_lfc = res$lfc

tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result")


tab_se = res$se
tab_se %>% 
  datatable(caption = "SEs from the Primary Result")

tab_w = res$W
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result")

tab_p = res$p_val
tab_p %>% 
  datatable(caption = "P-values from the Primary Result")

tab_q = res$q
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result")

tab_diff = res$diff_abn
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

#Bias corrected abundances
samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")

res = as.data.frame(res)

a = ggplot(res, aes(tab_lfc, color = Area ))+
  geom_point() +theme_bw() +
  geom_segment(aes(x=0, xend=lfc, y=name, yend=name, color = Area)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +ylab(NULL)











#Graph not working....

df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_area = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, AreaSF, SE) %>%
  dplyr::arrange(desc(AreaSF)) %>%
  dplyr::mutate(direct = ifelse(Area == "ML", "ML", "SF"))
df_fig_area$taxon_id = factor(df_fig_area$taxon_id, levels = df_fig_area$taxon_id)
df_fig_area$direct = factor(df_fig_area$direct, 
                           levels = c("ML", "SF"))

p_age = ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = age, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = age - ageSE, ymax = age + ageSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_age

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

#make a new column with the means of each row
SF$Colony = rowMeans((SF))
ML$Colony = rowMeans((ML))

#make a dataframe out of the means for each area and add rownames
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


#----------------------Volcano Plot--------------------------------------------

# https://www.youtube.com/watch?v=sIRnaKo1aKE

data = read_tsv("/home/amyparsons/data.tsv")

plot(data[,2:3])

sig.w = data$W >125 #arbitrarily setting sig W value
data$sig.w = sig.w
volc = ggplot(data) +
geom_point(aes(x=clr, y=W, colour = sig.w)) +
ggtitle("Differentially abundant Genus in Mouth Samples of ML and SF CAGU") +
 xlab("clr") + ylab("W") +
  theme_minimal()
  
# add a column of NAs
data$diffexpressed = "NO"
# if clr > 2.5 and W < 500, set as "SF Bay" 
data$diffexpressed[data$clr > 2.25 & data$W > 100] <- "SF Bay"
# if clr < -2.5 and pvalue < 500, set as "ML"
data$diffexpressed[data$clr < -2.5 & data$W > 100] <- "Mono Lake"

# Re-plot but this time color the points with "diffexpressed"
p = ggplot(data=data, aes(x=clr, y= W, colour=diffexpressed)) + geom_point() + theme_minimal()
p3 = p + scale_color_manual(values=c("blue", "black", "red"))
p3

#Create an column full of NAs
data$label = NA

#Create labels
data$label[data$id == "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Varibaculum"] = "(f) Actinomycetaceae"

data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Cardiobacteriales;f__Cardiobacteriaceae;g__Cardiobacterium"] = "(f) Cardiobacteriaceae"

data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio"] = "(f) Vibrionaceae"

data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Oceanospirillales;f__Nitrincolaceae;g__Nitrincola"] = "(f) Nitrincolaceae"

ggplot(data=data, aes(x= clr, y = W, col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
#  geom_label_repel(force_pull = 100, point.padding = 5, grid::arrow()) +
  ggtitle("Volcano plot of differentially abundant taxa (genus) in Mono Lake and SF Bay")

# Make a table with significantly differentially abundant taxa

#Make a column called "area" with NAs
Data$Area = "NO"
a
#Add labels to all significantly differentially abundant bacteria
data$area[data$clr > 0 & data$sig.w == "TRUE" ] <- "SF Bay"
data$area[data$clr < 0 & data$sig.w == "TRUE" ] <- "Mono Lake"



-----------------
  



