#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=train
#SBATCH --error=rescript.log
#SBATCH --output=rescript.log

date && \
echo "Activating qiime environemnt" && \
source activate qiime2-2022.8 && \
echo "start rescript" && \
qiime rescript get-silva-data --p-version '138' --p-target 'SSURef_NR99' --p-include-species-labels --o-silva-sequences silva-138-ssu-nr99-seqs.qza --o-silva-taxonomy silva-138-ssu-nr99-tax.qza && \
echo "finished rescript" && \
\
date
