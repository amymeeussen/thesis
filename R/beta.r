library('qiime2R')
library('ggplot2')

ps = qza_to_phyloseq(features = "/home/aparsons/table-filtered.qza", 
                     tree = "/home/aparsons/PhylogeneticTree/SILVAtree.qza",
                     taxonomy = "/home/aparsons/taxonomy.qza",
                     metadata = "/home/aparsons/metadata.tsv")
ps


library('phyloseq')

#explore phyloseq object

nsamples(ps)
sample_names(ps)
sample_variables(ps)
head(sample_data(ps))
table(sample_data(ps)$colony)
metadata <- data.frame(sample_data(ps))
head(metadata)
sample_sums(ps)
sort(sample_sums(ps))
metadata$total_reads <- sample_sums(ps)

#Examine ASV table

ntaxa(ps)
head(taxa_names(ps))
head(taxa_sums(ps))
(asv_tab <- data.frame(otu_table(ps)[1:5, 1:5]))

#Examine taxonomy

rank_names(ps)
rank_names(ps)
head(tax_table(ps)[, 2])
table(tax_table(ps)[, 2])

#Subsetting taxa

(ps_phylum <- tax_glom(ps, "Phylum"))
taxa_names(ps_phylum)
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, 2]
taxa_names(ps_phylum)
otu_table(ps_phylum)[1:5, c(1:3, 5, 7)]

#Relative abundance counts

ps_relabund <- transform_sample_counts(ps, function(x) x / sum(x))
otu_table(ps_relabund)[1:5, 1:5]

#Rarify

(ps_rare <- rarefy_even_depth(ps, sample.size = 4000, rngseed = 123, replace = FALSE))
sample_sums(ps_rare)

#Alpha diversity

png("richness.png", 490, 350)
head(estimate_richness(ps))
(p <- plot_richness(ps, x = "Area", color = "Area", measures = c("Observed", "Shannon")))
dev.off()

#Beta diversity

png("bray.png", 490, 350)
ps_rare_bray <- ordinate(ps_rare, "NMDS", "bray")
plot_ordination(ps_rare, ps_rare_bray, type="samples", color="colony") + geom_point(size = 3) 
dev.off()

png("bar.png", 1200, 1600)
plot_bar(ps, fill="Phylum")
plot_bar = plot_bar(ps_relabund, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank())
dev.off()

#Heatmap

png("heatmap.png", 690, 550)
(ps_fam_rare <- rarefy_even_depth(ps_rare, sample.size = 4000, rngseed = 123, replace = FALSE))
plot_heatmap(ps_fam_rare, sample.label = "Area", taxa.label = "Family")
dev.off()

