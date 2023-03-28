#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --partition=nodes
#SBATCH --job-name=tree
#SBATCH --error=convert.log
#SBATCH --output=convert.log

date \
&& \
echo "Activating qiime environemnt" \
&& \
source activate qiime2-2022.8 \
&& \
echo "conversion" && \
for i in *.qza; do
qiime tools export --input-path $i --output-path .
done
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
date \
