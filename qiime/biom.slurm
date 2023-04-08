#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=dem_den_tax
#SBATCH --error=/home/aparsons/qiime/biom.log
#SBATCH --output=/home/aparsons/qiime/biom.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8
qiime tools import --input-path ~/qiime/otu_biom.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ~/qiime/denoising/feature-table2-decontam.qza 
