library(phyloseq)
library(decontam)
library(tidyverse)
library(DT)
library(qiime2R)

# Use the line below to debug reading the metadata file
# meta = read_q2metadata("~/thesis/metadata/metadata_phyloseq.tsv")

ps = qza_to_phyloseq(
  features = "~/denoising/dada2-table.qza",
#  tree = "~/PhylogeneticTree/SILVAtree.qza",
#  taxonomy = "~/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

datatable(sample_data(ps))

df = as.data.frame(sample_data(ps))
df$LibrarySize = sample_sums(ps)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

a.ord = ordinate(ps, "PCoA", "bray")
plot_ordination(ps, a.ord, type="samples", color="Sample_or_control")

sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.01 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev.01$contaminant)

# result = FALSE:21248  TRUE:326 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.05 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev.05$contaminant)
View(contamdf.prev.05)

#result = FALSE: 20222 TRUE:1352

# View list of contaminents
View(contamdf.prev.05)

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
ps_no_control_no_contaminents = prune_taxa(!contamdf.prev.05$contaminant, ps_no_controls)

# Extract negative control samples 
ps_no_control_no_contaminents_new = subset_samples(ps_no_control_no_contaminents, Sample_or_control=="sample")

write.csv(otu_table(ps_no_control_no_contaminents_new),file = "ps_no_control_no_contaminents.csv")
