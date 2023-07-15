# This Rscript makes a linear regression of body condition and diversity. 

# Input: phyloseq object

# Output: Linear regression test

library(phyloseq)
library(tidyverse)
library(DT)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ANCOMBC)
library(mia)
library(ade4)
library(dplyr)
library(MicEco)
library(breakaway)
library(microbiome)

# Create a phyloseq object
ps = qza_to_phyloseq(
  features = "~/qiime/denoising_ss/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree_ss/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy_ss/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

# filter out negative controls
ps.no = subset_samples(ps, Area %in% c("ML", "SF"))

# subset individual areas
ps_ML = subset_samples(ps, Area %in% "ML")
ps_SF = subset_samples(ps, Area %in% "SF")

# Test for Shannon diversity with phyloseq package
sd = estimate_richness(ps.no, measures = "Shannon")
sd

Ob = estimate_richness(ps.no, measures = "Observed")
Ob
Ob.df = as.data.frame(Ch)

# FPD = picante::pd(samp, tree, include.root=TRUE)

# create a dataframe with diversity metrics and groups you want to test for diversity. You will need this
# test variance of diversity metrics between groups. 

diversity = as.data.frame(sd)
meta = meta(ps.no)
diversity$Area = meta$Area
diversity$sex = meta$sex
diversity$condition = meta$body_condition
diversity$Choa1 = Chao.df$Chao1
diversity$type = meta$type

#----------------------------------------Linear Regression--------------------------------

#Assumptions of the linear regression test

# 1. Linear relationship
# 2. Normality
# 3. No multicolinearity
# 4. Equal variances


#test for normality
shapiro.test(diversity$condition)

diversity = na.omit(diversity)
hist(diversity$condition)

#use levene's test to test for equal variance, since data is not normal
leveneTest(Shannon ~ interaction(Area), data = diversity)

#Test for equal variance (2nd assumption of the wilcoxon rank sum test)
var.test(diversity$Shannon ~ diversity$Area, alternative = "two.sided")

# Remove birds that don't have bc or sex data
diversity$Bird = meta$Bird
diversity_zbc = diversity %>% filter(!Bird %in% c("41", "37", "2", "19", "57"))

#test for equal variance between body conditions and sexes

#var.test(diversity_zbc$Shannon ~ diversity_zbc$bc, alternative = "two.sided")
#var.test(diversity_zbc$Shannon ~ diversity_zbc$sex, alternative = "two.sided")
 
# Check to see if there's a linear relationship
scatter.smooth(x=diversity$condition, y=diversity$Shannon)
#yes
scatter.smooth(x=diversity$condition, y=diversity$Chao1)
#no

# Check for outliers
par(mfrow=c(1, 2))  # divide graph area in 2 columns
boxplot(diversity$Shannon, sub=paste("Outlier rows: ", boxplot.stats(diversity$Shannon)$out))  
boxplot(diversity$condition, sub=paste("Outlier rows: ", boxplot.stats(diversity$condition)$out))

# Look at histograms to check for normality

hist(ML$condition)
hist(SF$condition)

# Shapiro test for normality
shapiro.test(sd$Shannon) # normal
shapiro.test(Ch$Chao1) # not normal

# Check for correlation 
cor.test(diversity$condition, diversity$Shannon) 

# Build linear model
L_model = lm(Shannon ~ condition, data=diversity)
summary(L_model)


result <- t.test(sample1, sample2, alternative = "two.sided", var.equal = TRUE)

#-------------------------------Equal variance in Mouth samples------------------------

# Assumptions of levenes test: 
# 1. independent observations
# 2. variables are quantitative


diversity_M= filter(diversity, type == "M")

leveneTest(Shannon ~ interaction(Area), data = diversity_M)


