#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=biom
#SBATCH --error=/home/aparsons/qiime/log/biom.log
#SBATCH --output=/home/aparsons/qiime/log/biom.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8

#Feature table import
#qiime tools import --input-path ~/thesis/R/otu_biom.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ~/qiime/denoising/taxa-table-decontam.qza && \ 
\
#Taxonomy table import
qiime tools import --type 'FeatureData[Taxonomy]’ --source-format HeaderlessTSVTaxonomyFormat --input-path ~/thesis/R/metadata/taxa.txt --output-path ~/qiime/taxonomy/taxa-table-filtered.qza && \
date \
