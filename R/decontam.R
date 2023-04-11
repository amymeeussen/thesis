library(phyloseq)
library(decontam)
library(tidyverse)
library(DT)
library(qiime2R)

# Use the line below to debug reading the metadata file
# meta = read_q2metadata("~/thesis/metadata/metadata_phyloseq.tsv")

ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
#  tree = "~/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

datatable(sample_data(ps))

# Compute library size
df = as.data.frame(sample_data(ps))
df$LibrarySize = sample_sums(ps)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

# Get the PCoA of the data
a.ord = ordinate(ps, "PCoA", "bray")
plot_ordination(ps, a.ord, type="samples", color="Sample_or_control")

# 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.01 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev.01$contaminant)
View(contamdf.prev.01)

# result = FALSE:21248  TRUE:326 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.05 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev.05$contaminant)
View(contamdf.prev.05)
# result = FALSE: 20222 TRUE:1352

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa = transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg = prune_samples(sample_data(ps.pa)$Sample_or_control == "control", ps.pa)
ps.pa.pos = prune_samples(sample_data(ps.pa)$Sample_or_control == "sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa = data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev.05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove controls
ps_no_controls = subset_samples(ps, Sample_or_control != "control")

# Identify and extract blank contaminants using the negative controls. The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants
sample_data(ps_no_controls)$is.neg <- sample_data(ps_no_controls)$Sample_or_control == "control"
contamdf.prev.05 = isContaminant(ps_no_controls, method="prevalence", neg="is.neg", threshold=0.5)

# Extract contaminants
ps_no_control_no_contaminants = prune_taxa(!contamdf.prev.05$contaminant, ps_no_controls)

# Extract negative control samples 
ps_no_control_no_contaminants_new = subset_samples(ps_no_control_no_contaminants, Sample_or_control=="sample")
write.csv(otu_table(ps_no_control_no_contaminants_new),file = "~/thesis/metadata/ps_no_control_no_contaminants.csv")

# Pull taxonomic table out of phyloseq object
tax = as(tax_table(ps), "matrix")
View(tax)

# Merge taxonomy table and contaminant table by sequence
tax = merge(tax, contamdf.prev.01, by="row.names", all=TRUE)
View(tax)

# Make a df of contaminants only
contaminants = filter(tax, contaminant == "TRUE")
View(contaminants)

#histogram of the prevelance of each contaminant in all of the samples
hist(contaminants$prev)

#histogram of the frequency of each contaminant in all of the samples
hist(contaminants$freq)

#remove all sequences that were classified as contaminants from the feature table
final_biom <- prune_taxa(!contamdf.prev.01$contaminant, ps)
final_biom

#--------------phyloseq to Qiime2---------------------


# Export taxonomy table as "tax.txt"
tax = as(tax_table(final_biom),"matrix")
tax_cols = colnames(tax)
tax = as.data.frame(tax)
tax$taxonomy = do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "tax.txt", quote=FALSE, col.names=FALSE, sep="\t")

# Export feature/OTU table as a biom file

library(biomformat);packageVersion("biomformat")

otu = as(otu_table(final_biom),"matrix") # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
#otu<-as(otu_table(globalpatterns),"matrix"))
otu_biom = make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")

#import re-seq file from qiime

seqs = read_qza("~/qiime/denoising/rep-seqs.qza")

