#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=clean
#SBATCH --error=cleanup.log
#SBATCH --output=cleanup.log

date && \
echo "Activating qiime environment" && \
source activate qiime2-2022.8 && \
echo "starting sample filtering" && \
qiime feature-table filter-samples --i-table ~/DADA2/dada2-table.qza --m-metadata-file ~/thesis/metadata/Samples_to_keep.tsv --o-filtered-table ~/filtered/id-filtered-table.qza
echo "finish sample filtering" \
&& \
echo "starting clean" && \
qiime taxa filter-table --i-table ~/filtered/id-filtered-table.qza/ --i-taxonomy ~/taxonomy/taxonomy.qza --p-include p_ --p-exclude 'p_; ,Mitochondria,Chloroplast,Eukaryota' --o-filtered-table ~/filtered/table-filtered.qza && \
echo "finish clean" \
&& \
echo "start reclassification" && \
qiime feature-table filter-seqs --i-data ~/DADA2/rep-seqs-dada2.qza --i-table ~/filtered/table-filtered.qza --o-filtered-data ~/filtered/filtered-sequences.qza && \
echo "finish reclassification" \
&& \ 
echo "start taxa barplot" && \
qiime taxa barplot --i-table ~/filtered/table-filtered.qza --i-taxonomy ~/taxonomy/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/filtered/taxa-bar-plot-filtered.qzv && \
echo "finished taxa barplot" \
