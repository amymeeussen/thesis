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

# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#This Rscript creates manipulate-able volcano plots showing the differential abundance of 
#microbes in each area
#input: .tsv file created from qiime2 view website after uploading .qzv file from qiime2 ANCOM results
#output: volcano plots!! hooray!



#----------------------Volcano Plot: Cloaca--------------------------------------------


data = read_tsv("/home/amyparsons/thesis/R/cloaca_ancom_data.tsv")

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
data$diffexpressed[data$clr < -2.25 & data$W > 100] <- "Mono Lake"

# Re-plot but this time color the points with "diffexpressed"
p = ggplot(data=data, aes(x=clr, y= W, colour=diffexpressed)) + geom_point() + theme_minimal()
p3 = p + scale_color_manual(values=c("blue", "black", "red"))
p3

#Create an column full of NAs
data$label = NA

#Create labels
data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio	"] = "(g) Vibrio"

data$label[data$id == "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus	"] = "(g) Streptococcus"

data$label[data$id == "d_Bacteria;p_Fusobacteriota;c_Fusobacteriia;o_Fusobacteriales;f_Fusobacteriaceae;g_Cetobacterium "] = "(g) Cetobacterium"


ggplot(data=data, aes(x= clr, y = W, col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
    geom_label_repel(force_pull = 100, point.padding = 5, grid::arrow()) +
  ggtitle("Volcano plot of differentially abundant taxa (genus) in CAGU cloaca at Mono Lake and SF Bay")


#----------------------Volcano Plot: Mouth--------------------------------------------


data = read_tsv("/home/amyparsons/thesis/R/mouth_ancom_data.tsv")

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
data$diffexpressed[data$clr < -2.25 & data$W > 100] <- "Mono Lake"

# Re-plot but this time color the points with "diffexpressed"
p = ggplot(data=data, aes(x=clr, y= W, colour=diffexpressed)) + geom_point() + theme_minimal()
p3 = p + scale_color_manual(values=c("blue", "black", "red"))
p3

#Create an column full of NAs
data$label = NA

#Create labels
data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Cardiobacteriales;f__Cardiobacteriaceae;g__Cardiobacterium"] = "(g) Cardiobacterium"

data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio"] = "(g) Vibrio"

data$label[data$id == "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Alloprevotella"] = "(g) Alloprevotella"

data$label[data$id == "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Aerococcaceae;__"] = "(f) Aerococcaceae"

ggplot(data=data, aes(x= clr, y = W, col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  geom_label_repel(force_pull = 100, point.padding = 5, grid::arrow()) +
  ggtitle("Volcano plot of differentially abundant taxa (genus) in CAGU mouth at Mono Lake and SF Bay")

#----------------------Volcano Plot: Foot--------------------------------------------


data = read_tsv("/home/amyparsons/thesis/R/foot_ancom_data.tsv")

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
data$diffexpressed[data$clr < -2.25 & data$W > 100] <- "Mono Lake"

# Re-plot but this time color the points with "diffexpressed"
p = ggplot(data=data, aes(x=clr, y= W, colour=diffexpressed)) + geom_point() + theme_minimal()
p3 = p + scale_color_manual(values=c("blue", "black", "red"))
p3

#Create an column full of NAs
data$label = NA

#Create labels
data$label[data$id == "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__Salinimicrobium	"] = "(g) Salinimicrobium"

data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Ectothiorhodospirales;f__Ectothiorhodospiraceae;g__Thioalkalivibrio	"] = "(g) Thioalkalivibrio"

data$label[data$id == "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Oceanospirillales;f__Halomonadaceae;g__Marinospirillum	"] = "(g) Marinospirillum"

data$label[data$id == "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Virgibacillus	"] = "(g) Virgibacillus"

ggplot(data=data, aes(x= clr, y = W, col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  geom_label_repel(force_pull = 100, point.padding = 5, grid::arrow()) +
  ggtitle("Volcano plot of differentially abundant taxa (genus) in CAGU mouth at Mono Lake and SF Bay")



