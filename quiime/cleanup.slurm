#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=clean
#SBATCH --error=clean.log
#SBATCH --output=clean.log

date && \
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "starting clean" && \
qiime taxa filter-table --i-table DADA2/dada2-table.qza --i-taxonomy taxonomy.qza --p-include p_ --p-exclude 'p_; ,Mitochondria,Chloroplast,Eukaryota' --o-filtered-table table-filtered.qza && \
echo "finish clean" \
&& \
echo "start reclassification" && \
qiime feature-table filter-seqs --i-data DADA2/rep-seqs-dada2.qza --i-table table-filtered.qza --o-filtered-data filtered-sequences.qza && \
echo "finish reclassification" \
&& \ 
echo "start taxa barplot" && \
qiime taxa barplot --i-table table-filtered.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.tsv --o-visualization taxa-bar-plot-filtered.qzv && \
echo "finished taxa barplot" \
