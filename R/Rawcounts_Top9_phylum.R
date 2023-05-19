library(phyloseq)
library(tidyverse)
library(DT)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ANCOMBC)
library(mia)



# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#filter out negative controls
ps.no = subset_samples(ps, Area %in% c("ML", "SF"))


#Get count of phyla
table(phyloseq::tax_table(ps.no)[, "Phylum"])


ps_phylum <- phyloseq::tax_glom(ps.no, "Phylum")
phyloseq::taxa_names(ps_phylum) = phyloseq::tax_table(ps_phylum)[, "Phylum"]



filter = phyloseq::genefilter_sample(ps_phylum, filterfun_sample(function(x) x >= 250), 
                                      A = 0.2*nsamples(ps_phylum))
ps_filtered <- prune_taxa(filter, ps_phylum)

phyloseq::otu_table(ps_filtered)[1:9, 1:9]


phyloseq::psmelt(ps_filtered) %>%
  ggplot(data = ., aes(x = Area, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") + 
  ggtitle("Raw Count data of top 9 Phyla in both Mono Lake and SF Bay")



