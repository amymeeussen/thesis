B#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=demultiplex
#SBATCH --error=~/qiime/dem_den_tax.log
#SBATCH --output=~/qiime/dem_den_tax.log

echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \


echo "Starting demultiplexing" && \
rm -rf ~/qiime/demux/ && \
mkdir -p ~/qiime/demux && \
date && \
qiime demux emp-paired --m-barcodes-file ~/thesis/metadata/metadata_noenv.tsv --m-barcodes-column barcodes --i-seqs ~/qiime/lab_data/emp-paired-end-sequences.qza --o-per-sample-sequences ~/qiime/demux/demux.qza --o-error-correction-details ~/qiime/demux/demux-details.qza --p-no-golay-error-correction && \
qiime demux summarize --i-data ~/qiime/demux/demux.qza --o-visualization ~/qiime/demux/demux.qzv && \
echo "finished demultiplexing" && \
\
\
\
echo "Starting denoising" && \
rm -rf ~/qiime/denoising/ && \
mkdir -p ~/qiime/denoising && \
date && \
qiime dada2 denoise-paired --i-demultiplexed-seqs ~/qiime/demux/demux.qza --p-n-threads 24 --p-trim-left-f 13 --p-trim-left-r 7 --p-trunc-len-f 170 --p-trunc-len-r 169 --o-table ~/qiime/denoising/dada2-table.qza --o-representative-sequences ~/qiime/denoising/rep-seqs-dada2.qza --o-denoising-stats ~/qiime/denoising/denoising-stats.qza && \
echo "Finished denoising" && \
echo "Starting denoisingstats" && \
date && \
qiime feature-table summarize --i-table ~/qiime/denoising/dada2-table.qza --m-sample-metadata-file ~/thesis/metadata/metadata_noenv.tsv --o-visualization ~/qiime/denoising/dada2-table.qzv && \
qiime feature-table tabulate-seqs --i-data ~/qiime/denoising/rep-seqs-dada2.qza --o-visualization ~/qiime/denoising/rep-seqs.qzv && \
qiime metadata tabulate --m-input-file ~/qiime/denoising/denoising-stats.qza --o-visualization ~/qiime/denoising/denoising-stats.qzv && \
echo "Finished denoisingstats" && \
\
\
\
echo "starting taxonomy assignment" && \
rm -rf ~/qiime/taxonomy/ && \
mkdir -p ~/qiime/taxonomy && \
date && \
qiime feature-classifier classify-sklearn --i-classifier ~/qiime/train_classifier/classifier.qza --i-reads ~/qiime/denoising/rep-seqs-dada2.qza --verbose --o-classification ~/qiime/taxonomy/taxonomy.qza && \
echo "finish taxonomy assignment" && \
date && \
echo "start visualization" && \
qiime metadata tabulate --m-input-file ~/qiime/taxonomy/taxonomy.qza --o-visualization ~/qiime/taxonomy/taxonomy.qzv && \
echo "finished visualization" && \
date && \
echo "start barplot" && \
qiime taxa barplot --i-table ~/qiime/denoising/dada2-table.qza --i-taxonomy ~/qiime/taxonomy/taxonomy.qza --m-metadata-file ~/thesis/metadata/metadata_noenv.tsv --o-visualization ~/qiime/taxonomy/taxa-bar-plots.qzv && \
echo "finished barplot" && \
date\


