

ps = qza_to_phyloseq(
  features = "~/qiime/denoising/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#subset variables

sample_variables(ps)

ps_mouth = subset_samples(ps, type == "M")

ps_foot = subset_samples(ps, type == "F")

ps_cloaca = subset_samples(ps, type == "C")

ps_samples_only = subset_samples(ps, type == c("C", "F", "M"))

# check to make sure subsetting worked

with(sample_data(ps), table(type == "M"))
ps_mouth

with(sample_data(ps), table(type == "F"))
ps_foot

with(sample_data(ps), table(type == "C"))
ps_cloaca

with(sample_data(ps), table(type == c("C", "F", "M")))
ps_samples_only

#---------------------rarefy----------------------------

rarecurve(t(otu_table(ps)), step=50, cex=0.5)

#check if taxa area rows. 
taxa_are_rows(otu_table(ps))

#vegan package requires samples in rows and taxa in columns. Transpose data:
trans = t(otu_table(ps))

#turn otu table into a matrix
trans_mat = as(t(otu_table(ps)), "matrix")
class(trans_mat)

raremax = min(rowSums(trans_mat))

system.time(rarecurve(trans_mat, step = 100, sample = raremax, col = "blue", label = FALSE))

ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)



#plot abundances

#ps_fam_mouth = tax_glom(ps_mouth, "Phylum")

#ps_fam_rel_mouth = transform_sample_counts(ps_fam_mouth, function(x) x/sum(x))

#couldn't get ps_prune to work
                             
#ps_df = psmelt(ps_fam_rel_mouth)

#p = ggplot(ps_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
#  theme_bw() +
#  geom_bar(stat = "identity") +
#  facet_grid(~ Area, space = "free", scales = "free") +
#  theme(axis.text.x = element_text(angle=90, size=6))
#p

