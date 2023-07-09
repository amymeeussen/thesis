# This script checks for assumption and runs the two way anova

# Input: phyloseq object: dada2-table.qza, SILVAtree.qza, taxonomy.qza, metadata_phyloseq.tsv

# Outputs: This will tell me if the effect of Area on diversity level is affected by body condition

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

ps = qza_to_phyloseq(
  features = "~/qiime/denoising_ss/dada2-table.qza",
  tree = "~/qiime/PhylogeneticTree_ss/SILVAtree.qza",
  taxonomy = "~/qiime/taxonomy_ss/taxonomy.qza",
  metadata = "~/thesis/metadata/metadata_phyloseq.tsv")

#filter out negative controls
ps_ML = subset_samples(ps, Area %in% c("ML"))

#------------------------Test for 2 way ANOVA assumptions--------------------

# This will tell me if the effect of Area on diversity level is affected by body condition

#Assumption 1:factor A and B are normally distributed (Area and body condition)
#Assumption 2: Factor A and B have equal variances 
#Assumption 3: There is no interaction between factors A and B
#Assumption 4: No outliers


#Test for Shannon diversity with phyloseq package 
sd = estimate_richness(ps.no_ss, measures = c("Shannon"))
sd


#Make a dataframe with variables you'll use in 2 way ANOVA

diversity = as.data.frame(sd)
meta = meta(ps_ML)
diversity$Area = meta$Area
diversity$bc = meta$condition_score
diversity$sex = meta$sex
diversity$condition = meta$body_condition
diversity$egg = meta$egg

#Replace empty values in the dataset with NAs so that you can then remove them
diversity[diversity == ""] = NA

#Remove all rows with NA in the bc column
diversity = diversity[complete.cases(diversity$bc), ]

#Check for equal means in 


#Check structure of the dataframe
str(diversity)

#Look at the balance of the design. It's unbalanced, so you have to use the sum of squares for the ANOVA and it's recommended to use Type III
table(diversity$bc, diversity$Area)

library(car)
my_anova <- aov(Shannon ~  bc + , data = diversity)
Anova(my_anova, type = "III")


# Result: Area has a significant effect on Shannon diversity (p <0.01), while body condition alone does not. However, the 
#interaction between body condition and Shannon diversity was significant (p<0.001)








