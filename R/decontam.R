library(phyloseq)
library(decontam)
library(tidyverse)
library(DT)
library(qiime2R)

#This file gives the code to decontaminate sequences using the fisher method. This method uses both 
#frequency and prevelance to give each sequence a score. The default threshold is p = 0.01. 
#Input: phylsoseq object that includes: dada2-table.qza, taxonomy.qza, metadata.tsv
#Outputs: .tsv file that lists all of the contaminants and a .biom file of the decontaminated sOTU 
#table. Import this file to the HPC and then use biom.slurm to convert this file to a .qza file for 
#downstream analysis in qiime2.


# Use the line below to debug reading the metadata file
# meta = read_q2metadata("~/thesis/metadata/metadata_phyloseq.tsv")


#Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  # tree = "~/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")


#Check the data
datatable(sample_data(ps))

# Graph the distribution of samples and controls by library size
df = as.data.frame(sample_data(ps))
df$LibrarySize = sample_sums(ps)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()



#-----------------Explore decontam (frequency method)----------------------

# Creates dataframe with features, contaminant TRUE/FALSE, and p value
contamdf.freq = isContaminant(ps, method="frequency", conc="Concentration")
head(contamdf.freq)

# Table tells how many contaminants (98 contaminants out of 22608)
table(contamdf.freq$contaminant)

# Lists the ranking of each contaminant (shows where contaminant falls in the list of most abundant to least)
head(which(contamdf.freq$contaminant))

# Plots frequency against DNA concentration, dotted line shows non-contaminant distribution, red line shows contaminant 
#distribution which would be inversely proportional to DNA concentration.
plot_frequency(ps, taxa_names(ps)[c(1,119)], conc="Concentration") + 
  xlab("DNA Concentration")

# Try a few more known contaminants to see if they show inversely proportional distribution
plot_frequency(ps, taxa_names(ps)[c(141,174)], conc="Concentration") + 
  xlab("DNA Concentration")

# Get the PCoA of the data
a.ord = ordinate(ps, "PCoA", "bray")
plot_ordination(ps, a.ord, type="samples", color="Sample_or_control")


#-------------------------Explore different methods-------------------------------

# Prevelance method: test for contaminants at threshold = 0.1 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.01 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev.01$contaminant)
View(contamdf.prev.01)
# result = FALSE:21248  TRUE:140 

# Prevelance method: test for contaminants at threshold=0.5 
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.prev.05 = isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev.05$contaminant)
View(contamdf.prev.05)
# result = FALSE: 20222 TRUE:1352

# Prevelance and Freq. method = "both", threshold=0.5
#Using method ="both" only identifies contaminants that are found by both prevelance and frequency methods
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.both.05 = isContaminant(ps, method="both", conc = "Concentration", neg="is.neg")
table(contamdf.both.05$contaminant)
View(contamdf.both.05)
# result = FALSE:22606 TRUE:2

# Prevelance and Freq. method = "combined", threshold=0.5
#Using method ="combined" uses the fisher method to give a contamination score using both freq. and prev.
#Used this method, since it was shown to be most robust in Davis et al.2018
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.comb.05 = isContaminant(ps, method="combined", conc = "Concentration", neg="is.neg")
table(contamdf.comb.05$contaminant)
View(contamdf.comb.05)
# result = FALSE:22506 TRUE:102


#------------------------Decontam Run---------------------------------------

# Prevelance and Freq. method = "combined", threshold=0.05
#Using method ="combined" uses the fisher method to give a contamination score using both freq. and prev.
#Used this method, since it was shown to be most robust in Davis et al.2018
sample_data(ps)$is.neg = sample_data(ps)$Sample_or_control == "control"
contamdf.comb.01 = isContaminant(ps, method="combined", conc = "Concentration", neg="is.neg")
table(contamdf.comb.01$contaminant)
View(contamdf.comb.01)
# result = FALSE:22506 TRUE:102

# Pull taxonomic table out of phyloseq object
tax = as(tax_table(ps), "matrix")
View(tax)

# Merge taxonomy table and contaminant table by sequence
tax2 = merge(tax, contamdf.comb.05, by="row.names", all=TRUE)
View(tax2)

# Make a df of contaminants only
contaminants = filter(tax2, contaminant == "TRUE")
View(contaminants)

#write the list of contaminants to a csv
write.csv(contaminants, "~/thesis/metadata/contaminants_05.csv", row.names=FALSE)

#Prune contaminant taxa out of OTU table
final_biom <- prune_taxa(!contamdf.comb.01$contaminant, ps)
final_biom

#-----------------------Phyloseq object to Qiime2---------------------------------

library(biomformat);packageVersion("biomformat")

#if taxa_are_rows=TRUE: otu = as(otu_table(globalpatterns),"matrix"))
#if taxa_are_rows=FALSE: t(otu_table(globalpatterns),"matrix"))
#Upload otu_biom.biom to HPC and run biom.slurm to convert into .qza file

otu = as(otu_table(final_biom),"matrix") 
otu_biom = make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")





