#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=diversity
#SBATCH --error=diversity.log
#SBATCH --output=diversity.log

date \
&& \
echo "Activating qiime environemnt" \
&& \
source activate qiime2-2022.8 \
&& \
echo "Starting measuring" \
&& \
qiime diversity core-metrics-phylogenetic \
--i-phylogeny SILVAtree.qza \
--i-table PhylogeneticTree/SILVA-removed-table.qza \
--p-sampling-depth 21000 \
--m-metadata-file Meta_HPC.tsv \
--output-dir core-metric-results \
&& \
date

