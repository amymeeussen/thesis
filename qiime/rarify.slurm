#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=rarify
#SBATCH --error=rarify.log
#SBATCH --output=rarify.log

date \
&& \
echo "Activating qiime environemnt" \
&& \
source activate qiime2-2022.8 \
&& \
echo "Starting rarifying" \
&& \
qiime diversity alpha-rarefaction --i-table table-filtered.qza --i-phylogeny PhylogeneticTree/SILVAtree.qza --m-metadata-file metadata.tsv --p-max-depth 28964 --p-steps 30 --p-iterations 20 --o-visualization alpha-rarefaction-rooted.qzv && \
date \

