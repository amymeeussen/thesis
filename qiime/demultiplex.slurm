#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=16G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=demultiplex
#SBATCH --error=demultiplex.log
#SBATCH --output=demultiplex.log


echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "Starting demultiplexing" && \
date && \
qiime demux emp-paired --m-barcodes-file ~/thesis/metadata/metadata.tsv --m-barcodes-column barcodes --i-seqs emp-paired-end-sequences.qza --o-per-sample-sequences demux.qza --o-error-correction-details demux-details.qza --p-no-golay-error-correction && \
echo "Finished demultiplexing" && \
echo "Starting visualization" && \
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv &&\
date \





