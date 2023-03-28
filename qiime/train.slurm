#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=train
#SBATCH --error=train.log
#SBATCH --output=train.log

date && \
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "start get silva data" && \
qiime rescript get-silva-data --p-version '138' --p-target 'SSURef_NR99' --p-include-species-labels --o-silva-sequences TrainClassifier/silva-rna-138-99-seqs.qza --o-silva-taxonomy TrainClassifier/silva-138-99-tax.qza && \
echo "finished get silva" && \
\
qiime rescript reverse-transcribe --i-rna-sequences TrainClassifier/silva-rna-138-99-seqs.qza --o-dna-sequences TrainClassifier/silva-138-99-seqs.qza && \
\
echo "start training" && \
qiime feature-classifier extract-reads --i-sequences TrainClassifier/silva-138-99-seqs.qza --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT --o-reads TrainClassifier/515-806-ref-seqs.qza && \
echo "finished training" && \
\
echo "starting this other thing" && \
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads TrainClassifier/515-806-ref-seqs.qza --i-reference-taxonomy TrainClassifier/silva-138-99-tax.qza --o-classifier TrainClassifier/classifier.qza && \
echo "finished this other thing" && \
date

