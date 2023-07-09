# This Rscript runs 3 methods for decontaminating samples (prevelance, frequency, and combined)

# Input: phyloseq object, DNA concentration list from Argonne labs, and a list of negative
#        lab controls

# Output: list of contaminants for each method with .01 and .05 threshold, a new feature table to be
#         reread back into QIIME2

library(phyloseq)
library(decontam)
library(tidyverse)
library(DT)
library(qiime2R)

# Use the line below to debug reading the metadata file
# meta = read_q2metadata("~/thesis/metadata/metadata_phyloseq.tsv")


ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
# tree = "~/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

datatable(sample_data(ps))

# Compute library size
df = as.data.frame(sample_data(ps))
df$LibrarySize = sample_sums(ps)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

#-----------------Frequency method----------------------

#Creates dataframe with features, contaminant TRUE/FALSE, and p value
contamdf.freq = isContaminant(ps, method="frequency", conc="Concentration")
head(contamdf.freq)

#Table tells how many contaminants (98 contaminants out of 22608)
table(contamdf.freq$contaminant)

#lists the ranking of each contaminant (shows where contaminant falls in the list of most abundant to least)
head(which(contamdf.freq$contaminant))

#plots frequency against DNA concentration, dotted line shows non-contaminent distribution, red line shows contaminent 
#distribution which would be inversley proportional to DNA concentration.
plot_frequency(ps, taxa_names(ps)[c(1,119)], conc="Concentration") + 
  xlab("DNA Concentration")

#try a few more known contaminants to see if they show inversely proportional distribution
plot_frequency(ps, taxa_names(ps)[c(141,174)], conc="Concentration") + 
  xlab("DNA Concentration")

# Get the PCoA of the data
a.ord = ordinate(ps, "PCoA", "bray")
plot_ordination(ps, a.ord, type="samples", color="Sample_or_control")

#-------------------------Explore the data-------------------------------

# Prevelance method: test for contaminants at p=0.01 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.05 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.01)
table(contamdf.prev.05$contaminant)
View(contamdf.prev.05)
# result = FALSE:21248  TRUE:140 

#Prevelance method: test for contaminants at p=0.05 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.05 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.05)
table(contamdf.prev.05$contaminant)
View(contamdf.prev.05)
# result = FALSE: 20222 TRUE:1352

# Prevelance and Freq. method = "both", threshold=0.05
#Using method ="both" only identifies contaminants that are found by both prevelance and frequency methods
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.both.05 = isContaminant(ps, method="both", conc = "Concentration", neg="is.neg")
table(contamdf.both.05$contaminant)
View(contamdf.both.05)
# result = FALSE:22606 TRUE:2

# Prevelance and Freq. method = "combined", threshold=0.05
#Using method ="combined" uses the fisher method to give a contamination score using both freq. and prev.
#Used this method, since it was shown to be most robust in Davis et al.2018
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.comb.05 = isContaminant(ps, method="combined", conc = "Concentration", neg="is.neg")
table(contamdf.comb.05$contaminant)
View(contamdf.comb.05)
# result = FALSE:22506 TRUE:102

#------------------------process data--------------------------------

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa = transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg = prune_samples(sample_data(ps.pa)$Sample_or_control == "control", ps.pa)
ps.pa.pos = prune_samples(sample_data(ps.pa)$Sample_or_control == "sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa = data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.comb.05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove controls
ps_no_controls = subset_samples(ps, Sample_or_control != "control")

# Identify and extract contaminants using the negative controls, combined method, threshold = 0.05
sample_data(ps_no_controls)$is.neg = sample_data(ps_no_controls)$Sample_or_control == "control"
contamdf.comb.05 = isContaminant(ps_no_controls, method="combined", neg="is.neg", conc = "Concentration")

# Extract contaminants
ps_no_control_no_contaminants = prune_taxa(!contamdf.comb.05$contaminant, ps_no_controls)

# Extract negative control samples and create phyloseq object
ps_no_control_no_contaminants_new = subset_samples(ps_no_control_no_contaminants, Sample_or_control=="sample")
write.csv(otu_table(ps_no_control_no_contaminants_new),file = "~/thesis/metadata/ps_no_control_no_contaminants.csv")

# Pull taxonomic table out of phyloseq object
tax = as(tax_table(ps), "matrix")
View(tax)

# Merge taxonomy table and contaminant table by sequence (this is where it gets weird)
tax2 = merge(tax, contamdf.comb.05, by="row.names", all=TRUE)
View(tax2)

# Make a df of contaminants only
contaminants = filter(tax2, contaminant == "TRUE")
View(contaminants)

#write the list of contaminants to a csv
write.csv(contaminants, "~/thesis/metadata/contaminants_05.csv", row.names=FALSE)

#histogram of the prevelance of each contaminant in all of the samples
#hist(contaminants$prev)

#histogram of the frequency of each contaminant in all of the samples
#hist(contaminants$freq)

#remove all sequences that were classified as contaminants from the feature table
final_biom = prune_taxa(!contamdf.comb.05$contaminant, ps)
#final_biom_nocontrols = subset_samples(final_biom, Sample_or_control=="sample")
#View(final_biom_no controls)
View(final_biom)


#--------------phyloseq to Qiime2---------------------


# Export taxonomy table as "tax.txt"
#tax = as(tax_table(final_biom),"matrix")
#tax_cols = colnames(tax)
#tax = as.data.frame(tax)
#tax$taxonomy = do.call(paste, c(tax[tax_cols], sep=";"))
#for(co in tax_cols) tax[co]<-NULL
#write.table(tax, "tax.txt", quote=FALSE, col.names=FALSE, sep="\t")

# Export feature/OTU table as a biom file

library(biomformat);packageVersion("biomformat")

#otu = as(otu_table(final_biom_nocontrols),"matrix") # 't' to transform if taxa_are_rows=FALSE
otu = as(otu_table(final_biom),"matrix") # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
#otu<-as(otu_table(globalpatterns),"matrix"))
otu_biom = make_biom(data=otu)
write_biom(otu_biom,"otu_biom_test.biom")


