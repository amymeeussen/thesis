#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=48:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=collapse
#SBATCH --error=/home/aparsons/qiime/log/collapse.log
#SBATCH --output=/home/aparsons/qiime/log/collapse.log

echo "Activating qiime environemnt"
source activate qiime2-2022.8

echo "Start collapsing" && \
qiime taxa collapse  --i-table ~/qiime/denoising/decontam_table_filtered.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza   --o-collapsed-table ~/qiime/taxonomy/deconotam_table-level6.qza --p-level 7 && \
qiime tools export --input-path /home/aparsons/qiime/taxonomy/deconotam_table-level6.qza --output-path /home/aparsons/qiime/
biom convert -i ~/qiime/feature-table.biom -o feature-table.tsv --to-tsv  && \
date \



