#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=demultiplex
#SBATCH --error=wholerun.log
#SBATCH --output=wholerun.log


echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "Starting demultiplexing" && \
date && \
qiime demux emp-paired --m-barcodes-file metadata.tsv --m-barcodes-column barcodes --i-seqs emp-paired-end-sequences.qza --o-per-sample-sequences demux.qza --o-error-correction-details demux-details.qza --p-no-golay-error-correction && \
echo "Finished demultiplexing" && \
echo "Starting visualization" && \
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv &&\
date && \
&& \
echo "Starting denoising" && \
date && \
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-n-threads 24 --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 232 --p-trunc-len-r 209 --o-table DADA2/dada2-table.qza --o-representative-sequences DADA2/rep-seqs-dada2.qza --o$
echo "Finished denoising" && \
echo "start visualization" && \
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv && \
date && \
&&\
echo "Starting denoisingstats" && \
date && \
qiime feature-table summarize --i-table DADA2/dada2-table.qza --m-sample-metadata-file metadata.tsv --o-visualization DADA2/dada2-table.qzv && \
qiime feature-table tabulate-seqs --i-data DADA2/rep-seqs-dada2.qza --o-visualization DADA2/rep-seqs.qzv && \
qiime metadata tabulate --m-input-file DADA2/denoising-stats.qza --o-visualization DADA2/denoising-stats.qzv && \
echo "Finished denoisingstats" && \
date && \
&& \
echo "starting taxonomy assignment" && \
qiime feature-classifier classify-sklearn --i-classifier TrainClassifier/classifier.qza --i-reads DADA2/rep-seqs-dada2.qza --verbose --o-classification taxonomy.qza && \
echo "finish taxonomy assignment" && \
\
echo "start visualization" && \
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv && \
echo "finished visualization" && \
\
echo "start barplot" && \
qiime taxa barplot --i-table DADA2/dada2-table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.tsv --o-visualization TaxonomicAssignment/taxa-bar-plots.qzv && \
echo "finished barplot" && \
date\

echo "Starting removing unwanted" && \
time qiime fragment-insertion sepp --i-representative-sequences filtered-sequences.qza --i-reference-database sepp-refs-silva-128.qza --o-tree PhylogeneticTree/SILVAtree.qza --o-placements SILVAtree_placements.qza --p-threads 10 && \
qiime fragment-insertion filter-features --i-table table-filtered.qza --i-tree PhylogeneticTree/SILVAtree.qza --o-filtered-table PhylogeneticTree/SILVA-filtered-table.qza --o-removed-table PhylogeneticTree/removed_table.qza && \
date \
