#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=denoising_sue
#SBATCH --error=denoising.log	
#SBATCH --output=denoising.log

mkdir -p DADA2 &&  \
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "Starting denoising" && \
date && \
qiime dada2 denoise-paired --i-demultiplexed-seqs ~/demux/demux.qza --p-n-threads 24 --p-trim-left-f 13 --p-trim-left-r 7 --p-trunc-len-f 170 --p-trunc-len-r 169 --o-table ~/DADA2/dada2-table.qza --o-representative-sequences ~/DADA2/rep-seqs-dada2.qza --o-denoising-stats ~/DADA2/denoising-stats.qza && \
echo "Finished denoising" && \
echo "Starting denoisingstats" && \
date && \
qiime feature-table summarize --i-table ~/DADA2/dada2-table.qza --m-sample-metadata-file ~/thesis/metadata/metadata.tsv --o-visualization ~/DADA2/dada2-table.qzv && \
qiime feature-table tabulate-seqs --i-data ~/DADA2/rep-seqs-dada2.qza --o-visualization ~/DADA2/rep-seqs.qzv && \
qiime metadata tabulate --m-input-file ~/DADA2/denoising-stats.qza --o-visualization ~/DADA2/denoising-stats.qzv && \
echo "Finished denoisingstats" && \
date
 
